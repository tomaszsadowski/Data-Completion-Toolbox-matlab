function [A,X,res] = nmf_grbf_admm(Y,A,X,tol_obj,MaxIter,lambda,show_inx)

% INPUTS:
% Y - mixted spectra (I x T), whre: columns - spectra, rows - samples,
% A - initial factors (I x J),
% X - initial factor (R x T),
% tol_proj - tolerance for stagnation of the residual error (e.g. 1e-4),
% MaxIter - maximal number of iterations,
% lambda - initial value for regularization parameter (e.g. 1e2)
% show_inx, if show_inx = 1, the residual error is displayed

% OUTPUTS:
% A - source spectra,
% X - mixing operator,
% res - residual error

% ==========================================================================



[I,T] = size(Y);
J = size(A,2);
M = I;
i = (1:I)'; m = 1:M;
res(1) = norm(Y - A*X,'fro')/norm(Y,'fro');

% B-splines
k_order = 4; Ix = I;
tau = linspace(0,Ix,I);
knots = augknt(linspace(0,Ix,10),k_order);
Phi = spcol(knots,k_order,tau); 
R = size(Phi,2);

F_1 = Phi'*Phi;
LambdaA = zeros(size(A));
LambdaX = zeros(size(X));
k_inner = 10;
H = X;

for k = 1:MaxIter
    
    if show_inx && (k < 10 || ~mod(k,10))
        disp(['ADMM-GRBF iterations: ',num2str(k)]);
    end
    
    lambda = max(1e-12,lambda/2);
    
    % B-splines
    mx = max(2,min(100,k));
    k_order = 6; Ix = I/4;
    tau = linspace(0,Ix,I);
    knots = augknt(linspace(0,Ix,mx),k_order);
    Phi = spcol(knots,k_order,tau); 
   % R = size(Phi,2);
    
    F = Phi'*Phi;
    Finv = inv(F);
    U = Finv*Phi';
    if k == 1
       W = Finv*Phi'*A;
    end
    
    %% FC-NNLS
   %[X, Pset] = fcnnls(A,Y,1e-12);
    X = fast_hals_inner(A,Y,X,size(X,1),1e-12,k_inner);
      
% % Normalization
    dX = eps+(sum(X,2));
    X = bsxfun(@ldivide,dX,X);
    A = bsxfun(@times,dX',A);
   
   %% ADMM
    [A,W,LambdaA] = admm_grbf_W(Y,W,X,Phi,A,U,LambdaA,k_inner);
    
   % Normalization
     dA = eps+(sum(A,1));
     A = bsxfun(@rdivide,A,dA);
     X = bsxfun(@times,dA',X);

    res(k+1) = norm(Y - A*X,'fro')/norm(Y,'fro');
    
  % Stagnation 
    if k > 20
       cx(k) = abs(res(k) - res(k+1))/res(1); 
       if cx(k) < tol_obj
          disp(['Stagnation of ADMM-GRBF iterations: ',num2str(k+1), ', Residual: ',num2str(res(k+1)), ', Diff_res: ', num2str(res(k) - res(k+1))]);
          break;
       end
    end % if


end
