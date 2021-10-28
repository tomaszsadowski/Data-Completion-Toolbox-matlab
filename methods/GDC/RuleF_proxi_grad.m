function [C] = RuleF_proxi_grad(Q,Z,C,q)
%function [C] = RuleF_proxi_grad(Q,Z,C,q)

J = size(Z,1);
k_inner = 10;

Z = tensor(Z); 
ZZ = double(ttt(Z,Z,[2,3],[2,3]));
QZ = double(ttt(tensor(Q),Z,[2,3],[2,3]));


k = 0; W = C; tau = 1; 
while (k <= k_inner)
        
        k=k+1;
        lambda = min(1,10^(-2 + q*.15));
        X = (2*C - W + tau*QZ)/(ZZ + eye(J));
        W = W + lambda*(X - C);
        C = max(0,W);  
        C = C + (ones(size(C,1),1) - sum(C,2))*ones(1,J)/J;

        
end % k_inner
      
end
