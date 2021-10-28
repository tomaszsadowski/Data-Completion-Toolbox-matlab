function [ A ]=HALS_IC(M, ohm, r,maxiter,er)
%HALS for image completion method
% Y = HALS_IC(M, ohm, r,maxiter,er)
% Inputs:
% M - incomplite data (image),
% ohm - a binary matrix of indicators for the missing entries,
% R - rank of factorization,
% Outputs:
% Y - completed data (image)
%
% Cite:
% Tomasz Sadowski, Rafal Zdunek, Modified HALS algorithm for image completion and recommendation system,
% In: Information Systems Architecture and Technology : Proceedings of 38th International Conference on Information
% Systems Architecture and Technology, ISAT 2017. Pt. 2 / eds. Jerzy Świątek, Leszek Borzemski, Zofia Wilimowska.
% Cham : Springer, cop. 2018, pp 17-27 (Advances in Intelligent Systems and Computing, ISSN 2194-5357; vol. 656)
% DOI:	doi.org/10.1007/978-3-319-67229-8_2
%
% Implemented by Tomasz Sadowski 2017

I=size(M);

if length(I)>2
    a=rand(I(1,1),r);
    a1=repmat(a,1,1,I(3));
    x=rand(r,I(1,2));
    x1=repmat(x,1,1,I(3));
    
    Z=M;
    Y=zeros(I);
    A=Y;
    
    for i=1:I(3)
        %separatly for each color
        e1=0;
        e2=1e4+1;
        %er=1e3;
        k=0;
        
        while (abs(e1-e2) > er)&&(k~=maxiter)
            disp(['HALS_IC: Color: ', num2str(i), ' , iteration: ',num2str(k), ', error: ',num2str(abs(e1-e2))])
            
            k = k + 1;
            
            e1=norm(Z(:,:,i),'fro')^2;
            
            Y(:,:,i)=Z(:,:,i)+a1(:,:,i)*x1(:,:,i);
            tmp= Y(:,:,i);
            tmp1=M(:,:,i);
            tmp(ohm(:,:,i)) = tmp1(ohm(:,:,i));
            
            Y(:,:,i)=tmp;
            % choose method
            
            [a1(:,:,i),x1(:,:,i)]=nmf_fast_hals(Y(:,:,i),a1(:,:,i),x1(:,:,i),10e-12,10e-6,100,0);
            
            Z(:,:,i)=Y(:,:,i)-a1(:,:,i)*x1(:,:,i);
            tmp1=Z(:,:,i);
            tmp1(~ohm(:,:,i)) = 0;
            Z(:,:,i)=tmp1;
            
            
            e2=norm(Z(:,:,i),'fro')^2;
        end
        
        
    end
    for i=1:I(3)
        A(:,:,i)=a1(:,:,i)*x1(:,:,i);
    end
        
    
else
    
    a=rand(I(1,1),r);
    x=rand(r,I(1,2));
    
    Z=M;
    e1=0;
    e2=1e4+1;
    a1=a;
    x1=x;
    k=0;
    
    while (abs(e1-e2) > er)&&(k<maxiter)
                
        k = k + 1;
        disp(['HALS_IC: ' , 'iteration: ',num2str(k), ' error: ',num2str(abs(e1-e2))])
        e1=norm(Z,'fro')^2;
                
        Y=Z+a1*x1;
        Y(ohm) = M(ohm);
        
        [a1,x1]=nmf_fast_hals(Y,a1,x1,10e-12,10e-6,100,0);    
        
               
        Z=Y-a1*x1;
        Z(~ohm) = 0;

        e2=norm(Z,'fro')^2;
              
    end 
    
    %final result
    A=a1*x1;
    
end

end