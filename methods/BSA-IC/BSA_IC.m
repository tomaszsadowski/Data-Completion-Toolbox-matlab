function [ A , res ]=BSA_IC(M, ohm, r,maxiter,er,X0)
%B-Splines-based Algorithm for Image Completion
% Y = BSA_IC(M, ohm, r,maxiter,er)
% Inputs:
% M - incomplite data (image), 
% ohm - a binary matrix of indicators for the missing entries,, 
% r - rank of factorization, 
% maxiter - maximal number of iterations, 
% er - threshold for stagnation 
% X0 - intial undisterbed data (image), to calculate residual error
% Outputs:
% A - completed data (image)
% res - residual error
%
% Cite:
% Tomasz Sadowski, Rafal Zdunek, Image completion with smooth nonnegative matrix factorization, 
% In: Artificial Intelligence and Soft Computing : 17th International Conference, 
% ICAISC 2018, Zakopane, Poland, June 3-7, 2018 : proceedings. Pt. 1 / eds. Leszek Rutkowski [i in.]. Cham : Springer, cop. 2018, pp 62-72
% (Lecture Notes in Computer Science, Lecture Notes in Artificial Intelligence, ISSN 0302-9743; vol. 10841)
% DOI:	doi.org/10.1007/978-3-319-91262-2_6
%
% Implemented by Tomasz Sadowski 2018

I=size(M);
res=zeros(1,maxiter);

if length(I)>2
    a=rand(I(1,1),r);
    a1=repmat(a,1,1,I(3));
    x=rand(r,I(1,2));
    x1=repmat(x,1,1,I(3));
    
    Z=M;
    Y=zeros(I);
    A=Y;
    
    for i=1:I(3)
        k=2;
        res(1)=inf+1;
        res(2)=inf;
        while (res(k) > er)&&(k~=maxiter)&& (res(k) <= res(k-1))
            disp(['BSA_IC: Color: ', num2str(i), ' , iteration: ',num2str(k), ', error: ',num2str(res(k))])
            
            k = k + 1;
            
            e1=norm(Z(:,:,i),'fro')^2;
            
            Y(:,:,i)=Z(:,:,i)+a1(:,:,i)*x1(:,:,i);
            tmp= Y(:,:,i);
            tmp1=M(:,:,i);
            tmp(ohm(:,:,i)) = tmp1(ohm(:,:,i));
            
            Y(:,:,i)=tmp;
           
            [a1(:,:,i),x1(:,:,i)]=nmf_grbf_admm(Y(:,:,i),a1(:,:,i),x1(:,:,i),10e-36,200,10e-12,0);
            
            Z(:,:,i)=Y(:,:,i)-a1(:,:,i)*x1(:,:,i);
            tmp1=Z(:,:,i);
            tmp1(~ohm(:,:,i)) = 0;
            Z(:,:,i)=tmp1;
            
            
            e2=norm(Z(:,:,i),'fro')^2;
            res(k)=abs(e1-e2);

        end
        
        
    end
    for i=1:I(3)
        A(:,:,i)=a1(:,:,i)*x1(:,:,i);
    end
        
    
else
    
    a=rand(I(1,1),r);
    x=rand(r,I(1,2));
    
    Z=M;
    a1=a;
    x1=x;
    k=1;
    res(1)=inf+1;
    res(2)=inf;
    while (res(k) > er)&&(k<maxiter) &&(res(k) <= res(k-1))
                
        k = k + 1;
   
        
        disp(['BSA_IC: ' , 'iteration: ',num2str(k), ' error: ',num2str(res(k)),' SIR: ',num2str(CalcSIR(X0,a1*x1))])
        e1=norm(Z,'fro')^2;
                
        Y=Z+a1*x1;
        Y(ohm) = M(ohm);
        
   
        [a1,x1]=nmf_grbf_admm(Y,a1,x1,10e-36,200,10e-12,0);

        
               
        Z=Y-a1*x1;
        Z(~ohm) = 0;

        e2=norm(Z,'fro')^2;
   
        
        res(k)=abs(e1-e2);
        

    end 
    
    %final result
    A=a1*x1;
    
end

end