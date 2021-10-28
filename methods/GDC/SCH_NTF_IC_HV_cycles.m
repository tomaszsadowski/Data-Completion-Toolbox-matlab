function [Y,res]= SCH_NTF_IC_HV_cycles(X,J,ohm,iter,er,X0)
% Segmented Convex Hull NTF for Image Completion, Horizontal-Vertical
% version
% [Y,res]= SCH_NTF_IC_HV_cycles(X,J,ohm,iter,er,X0)
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
DimY = size(Y);

% update z
for i=1:DimY(3)
    [Y(:,:,i)]=update_z_mean(Y(:,:,i),ohm_bar(:,:,i),1:size(Y,1));
end 
Y(ohm) = X(ohm);

Fh=rand(DimY(1),J);
Fh = diag(1./sum(Fh,2))*Fh;

Fv=rand(DimY(2),J);
Fv = diag(1./sum(Fv,2))*Fv;

q=0; MaxIter = 5;
r = 0; iter_q = 1;
inx_h = 0; inx_v = 0;
res=zeros(1,iter);
res(iter_q)=norm(tensor(Y-X0))/norm(tensor(X0));

  % while (q<iter)&&(res(end)<er)
for q = 1:MaxIter
          
        % Select extreme rays along the first mode  
        if inx_h < 1
           Kh = micro_model_means_n_way_slices(tensor(Y),J,2);      
      %  end
       
      %% updating F        
        kf = 0;         
        while (kf < 10)
            
            kf = kf + 1;
            iter_q = iter_q + 1;
            Z = Y(Kh,:,:);    
            Fh = RuleF_proxi_coordinate(Y,Z,Fh,iter_q);  
              
            %updating Y
            Y_old = Y;
            Yk=double(ttm(tensor(Z),Fh,1));
            Y(ohm_bar) = Yk(ohm_bar);
            res_ohm(iter_q) = norm(X(ohm) - Yk(ohm))/norm(X(ohm));
        
            %% update z
            for i=1:DimY(3)
                [Y(:,:,i)]=update_z_mean(Y(:,:,i),ohm_bar(:,:,i),1:size(Y,1));
            end 
            Y(ohm) = X(ohm);  
            
            res(iter_q)=norm(tensor(Y-X0))/norm(tensor(X0));
            disp(['SCH_NTF_IC_HV_h: iter: ' , num2str(iter_q), ', res_ohm: ' , num2str(res_ohm(iter_q))])            
               
            if iter_q > 2
                if ((res_ohm(iter_q-1) - res_ohm(iter_q)) < 0) 
                    Y = Y_old; 
                    res(iter_q)=res(iter_q-1);
                    inx_h = inx_h + 1;
                    break
                end   
            end
        end % while kf

        end % if inx_h                
        
        % Select extreme rays along the second mode 
        if inx_v < 1
            Kv = micro_model_means_n_way_slices(tensor(permute(Y,[2 1 3])),J,2);  
       
        
        %% updating F        
        kv = 0;        
        while kv < 10
            
            kv = kv + 1;
            iter_q = iter_q + 1;            
            Zv = Y(:,Kv,:);
            Fv = RuleF_proxi_coordinate(permute(Y,[2 1 3]),permute(Zv,[2 1 3]),Fv,iter_q); 
            
            %updating Y 
            Y_old = Y;
            Yk=double(ttm(tensor(Zv),Fv,2));
            Y(ohm_bar) = Yk(ohm_bar);
            res_ohm(iter_q) = norm(X(ohm) - Yk(ohm))/norm(X(ohm));        
        
            %% update z
            for i=1:DimY(3)
                [Y(:,:,i)]=update_z_mean(Y(:,:,i),ohm_bar(:,:,i),1:size(Y,1));
            end     
            Y(ohm) = X(ohm);  
            res(iter_q)=norm(tensor(Y-X0))/norm(tensor(X0));
            disp(['SCH_NTF_IC_HV_v: iter: ' , num2str(iter_q), ', res_ohm: ' , num2str(res_ohm(iter_q))])           
         
            if res_ohm(iter_q) > res_ohm(iter_q-1)
               Y = Y_old; 
               res(iter_q)=res(iter_q-1);    
               inx_v = inx_v + 1;
               break
            end 
            
        end  % kv   
        end % inx_v
       
        if (iter_q > iter) | (res(iter_q) < er) | (min(inx_h,inx_v) > 1)            
            break
        end
        
      %  res(q)=norm(tensor(Y-X0))/norm(tensor(X0));
      %  disp(['SCH_NTF_IC:: iter: ' , num2str(q), ', res: ' , num2str(res(q))])
%         if q>5
%             if res(q-1)<res(q)
%                 break;
%             end
%         end        
   
   end % while
