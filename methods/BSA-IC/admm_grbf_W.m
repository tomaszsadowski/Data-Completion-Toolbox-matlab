function [A,W,Lambda] = admm_grbf_W(Y,W,X,Phi,A,U,Lambda,k_inner)
%function [A,W,Lambda] = admm_grbf_W(Y,W,X,Phi,A,U,Lambda,lambda,k_inner)

   eta = 0.999; alpha = 1; Lambda_tilde = Lambda; c = 1;

   rho = 1/norm(X*X');
   Binv = inv(X*X' + rho*eye(size(X,1))); Yt = U*Y*X';
    
    for iter = 1:k_inner
       
       W = (Yt + U*(rho*A + Lambda))*Binv;
       A_tilde = Phi*W;
       A = max(eps,A_tilde - Lambda/rho);
       
       Lambda_old = Lambda;
       Lambda = Lambda + (A - A_tilde);
       
       A_old = A;

       c_old = c;
       c = (norm(Lambda - Lambda_tilde,'fro')^2)/rho + rho*norm(A - A_tilde,'fro')^2;
       if c < (eta*c_old)
           alpha_old = alpha;
           alpha = (1 + sqrt(1 + 4*alpha^2))/2;
           A_tilde = A + (A - A_old)*(alpha_old - 1)/alpha;
           Lambda_tilde = Lambda + (Lambda - Lambda_old)*(alpha_old - 1)/alpha;
       else
           alpha = 1;
           A_tilde = A_old;
           Lambda_tilde = Lambda_old;
           c = c_old/eta;
       end
                     
    end % for k

end % function spg_inner
