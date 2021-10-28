function [Y,res]= SCH_NTF_ICV(X,J,ohm,iter,er,X0)
% Segmented Convex Hull NTF for Image Completion, Vertical Version
% [Y,res]= SCH_NTF_ICV(X,J,ohm,iter,er,X0)
%Inputs:
% X - incomplite data (image), 
% J - rank of factorization
% ohm - a binary tensor of indicators for the missing entries,  
% % maxiter - maximal number of iterations, 
% er - threshold for stagnation 
% X0 - intial undisterbed data (image), to calculate residual error
%
% Outputs:
% Y_hat - completed data (image)
% res - residual error
% Requires: MATLAB Tensor Toolbox, Sandia Corporation
% (https://www.tensortoolbox.org)
%
% Cite:
% Rafal Zdunek and Tomasz Sadowski,
% Image completion with approximate convex hull tensor decomposition. 
% Signal Processing : Image Communication, 2021, vol. 95, art. 116276, pp 1-13.
% DOI: doi.org/10.1016/j.image.2021.116276	
%
% Implemented by Rafal Zdunek and Tomasz Sadowski 2018-2021
ohm_bar = not(ohm);
X(ohm_bar)=255;
Y=X;

%Size
DimX = size(X);

 % update z
 for i=1:DimX(3)
     [Y(:,:,i)]=update_z_mean(Y(:,:,i),ohm_bar(:,:,i),1:size(X,1));
 end 
 Y(ohm) = X(ohm);

F=rand(DimX(2),J);
F = diag(1./sum(F,2))*F;

q=1; 
res=zeros(1,iter);
res(q)=norm(tensor(Y-X0))/norm(tensor(X0));

   while (q<iter)&&(res(end)<er)
       
        q=q+1;
        
        %Select extreme rays    
        if (q == 2) || ~mod(q,10)
            K_final = micro_model_means_n_way_slices(tensor(permute(Y,[2 1 3])),J,2);  
        end
        
        % Update Z
        Z = Y(:,K_final,:);
           
      %% updating F by (7) modified
        F = RuleF_proxi_coordinate(permute(Y,[2 1 3]),permute(Z,[2 1 3]),F,q); 
        
        %updating Y by (8) 
        Y_old = Y;
        Yk=double(ttm(tensor(Z),F,2));
        Y(ohm_bar) = Yk(ohm_bar);  
        res_ohm(q) = norm(X(ohm) - Yk(ohm))/norm(X(ohm));
        disp(['SCH_NTF_IC_V: iter: ' , num2str(q), ', res_ohm: ' , num2str(res_ohm(q))])
        
       % update z
        for i=1:DimX(3)
            [Y(:,:,i)]=update_z_mean(Y(:,:,i),ohm_bar(:,:,i),1:size(X,1));
        end        
        Y(ohm) = X(ohm);        
       
        res(q)=norm(tensor(Y-X0))/norm(tensor(X0));
    %    disp(['SCH_NTF_IC_V:: iter: ' , num2str(q), ', res: ' , num2str(res(q))])
        if q>2
            if res_ohm(q-1) < res_ohm(q)  
                Y_tmp = Y_old;
                if max(res_ohm(q-2),res_ohm(q-1))<res_ohm(q)
                   res(q-1) = res(q-2); 
                   res(q) = res(q-2);    
                   Y = Y_tmp;
                   break
                end
            end
        end        
   
   end 
