function  [Y_hat]  = fALS_tucker_IC(X, J,ohm,maxiter,er)
%filtered ALS algorithm
%[Y_hat]  = fALS_tucker_IC(X, J,ohm,maxiter,er)
%
% Inputs:
% X - incomplite data (image), 
% J - rank of factorization,
% ohm - a binary matrix of indicators for the missing entries,  
% maxiter - maximal number of iterations, 
% er - threshold for stagnation 
%
% Outputs:
% Y_hat - completed data (image)
%
% Requires: MATLAB Tensor Toolbox, Sandia Corporation
% (https://www.tensortoolbox.org)
%
% Cite:
% Tomasz Sadowski, Rafal Zdunek, Image completion with filtered alternating least squares Tucker decomposition,
% In: Signal processing, algorithms, architectures, arrangements, and applications, SPA 2019 : conference processing, Poznań, 
% 18th-20th September 2019, [Danvers, MA : IEEE, cop. 2019], pp 241-245.
% DOI: doi.org/10.23919/SPA.2019.8936733	
%
% Implemented by Tomasz Sadowski 2019

%Size
DimX = size(X);
N = length(DimX);
if N==2
    J=[J J];
    D=1;
else
    if DimX(3)>3
        J=[min(J,DimX(1)) min(J,DimX(2)) min(J,DimX(3))];
    else
        J=[J J 3];
    end
    D=N;
end
%Initialization
sigm2=0.5;
res=zeros(1,maxiter);res(1)=inf;
Y_hat=X;

for i = 1:D
    Y_hat(:,:,i) = imgaussfilt(squeeze(Y_hat(:,:,i)),sigm2);
end

k=1;

%Initial guess for factors of T
[Y_hat,U0] = tucker_als(tensor(Y_hat),J,'maxiters',10,'tol',er,'printitn',0);
 
while (res(k) > er)&&(k<maxiter)
    
disp(['fALS_tucker_IC: ' , 'iteration: ',num2str(k), ' error: ',num2str(res(k))])

k=k+1;

[Y_hat,U0] = tucker_als(Y_hat,J,'init',U0,'maxiters',10,'tol',er,'printitn',0);

Y_hat=double(Y_hat);

for i = 1:D
    Y_hat(:,:,i) = imgaussfilt(squeeze(Y_hat(:,:,i)),sigm2);
end


Yres=Y_hat; Yres(~ohm)=0;

Y_hat(ohm)= X(ohm); 

Y_hat=tensor(Y_hat);
        
res(k)=norm(tensor(Yres-X))^2/norm(tensor(X(ohm)))^2;

if k>5
    if res(k-1)<res(k)
        break;
    end
end

end

Y_hat=double(Y_hat);
    
%end