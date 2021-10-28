function [C] = RuleF_proxi_coordinate(Q,Z,C,q)
%function [C] = RuleF_proxi_coordinate(Q,Z,C,q)

J = size(Z,1);
k_inner = 5;

Z = tensor(Z); 
ZZ = double(ttt(Z,Z,[2,3],[2,3]));
QZ = double(ttt(tensor(Q),Z,[2,3],[2,3]));


k = 0; W = C; tau = 1; 
while (k <= k_inner)
        
        k=k+1;
        lambda = min(1,10^(-2 + q*.15));
        X = 2*C - W;
        Xt = X;
        for kc = 1:1
            for i = 1:J-1            
                for j = i+1:J
               % j = i+1;
                    
                   ZZij = ZZ([i,j],[i,j]);
                   ZZij_det = (ZZij(1,1)+1)*(ZZij(2,2)+1) - ZZij(1,2)*ZZij(2,1);
                   ZZij_inv = [(ZZij(2,2)+1)/ZZij_det -ZZij(1,2)/ZZij_det; -ZZij(2,1)/ZZij_det (ZZij(1,1)+1)/ZZij_det];               
                   Xt(:,[i,j]) = (X(:,[i,j]) + Xt(:,[i,j])*ZZ([i,j],[i,j]) + tau*QZ(:,[i,j]) - Xt*ZZ(:,[i,j]))* ZZij_inv;
                
            
                end
            end
        end
        X = Xt;
        
        W = W + lambda*(X - C);
        C = max(0,W);  
        C = C + (ones(size(C,1),1) - sum(C,2))*ones(1,J)/J;
     %   Qm = ttm(Z,C,1);
     %   res(k) = norm(Q(:) - Qm(:))/norm(Q(:));
        
end % k_inner
      

end
