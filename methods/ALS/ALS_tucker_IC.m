function  [Y_hat]  = ALS_tucker_IC(X, J,ohm,maxiter,er)
%Image competion using tucker_als function
% [Y_hat]  = ALS_tucker_IC(X, J,ohm,maxiter,er)
% Inputs:
% X - incomplite data (image), 
% J - rank of factorization,
% ohm - a binary matrix of indicators for the missing entries,,  
% maxiter - maximal number of iterations, 
% er - residual error 
%
% Outputs:
% Y_hat - completed data (image)
%
% Requires: MATLAB Tensor Toolbox, Sandia Corporation
% (https://www.tensortoolbox.org)
%
% Cite:
% Tomasz Sadowski, Rafal Zdunek, Image completion with smooth nonnegative matrix factorization, 
% In: Artificial Intelligence and Soft Computing : 17th International Conference, 
% ICAISC 2018, Zakopane, Poland, June 3-7, 2018 : proceedings. Pt. 1 / eds. Leszek Rutkowski [i in.]. Cham : Springer, cop. 2018, pp 62-72
% (Lecture Notes in Computer Science. Lecture Notes in Artificial Intelligence, ISSN 0302-9743; vol. 10841)
% DOI: doi.org/10.23919/SPA.2019.8936733
%
% Implemented by Tomasz Sadowski 2019

%Size
DimX = size(X);
N = length(DimX);
if N==2
    J=[J J];
else
    if DimX(3)>3
        J=[J J J];
    else
        J=[J J 3];
    end
end
%Initialization

res=inf;
Y_hat=tensor(X);




k=1;

%Initial guess for factors of T
[Y_hat,U0] = tucker_als(Y_hat,J,'maxiters',10,'tol',er);
 
while (res > er)&&(k<maxiter)
    
disp(['ALS_tucker_IC: ' , 'iteration: ',num2str(k), ' error: ',num2str(res(end))])

k=k+1;

[Y_hat,U0] = tucker_als(Y_hat,J,'init',U0,'maxiters',10,'tol',er);

Y_hat=double(Y_hat);
Yres=Y_hat; Yres(~ohm)=0;

Y_hat(ohm)= X(ohm); 

Y_hat=tensor(Y_hat);
        
res=norm(tensor(Yres-X))^2/norm(tensor(X(ohm)))^2;

end

Y_hat=double(Y_hat);
    
%end