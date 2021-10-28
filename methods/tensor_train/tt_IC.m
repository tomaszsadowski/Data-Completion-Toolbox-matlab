function  [Y_hat]  = tt_IC(M, J, Omega, maxiter,er)
%Image competion using tensor train, mode permutation
%function  [Y_hat]  = tt_IC(M, J, Omega, maxiter,er)
%
% Inputs:
% X - incomplite data (image), 
% J - rank of factorization,
% Omega - a binary matrix of indicators for the missing entries,  
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
% Rafal Zdunek, Krzysztof Fonal and Tomasz Sadowski, Image completion with filtered low-rank tensor train approximations. 
% In: Advances in Computational Intelligence : 15th International Work-Conference on Artificial Neural Networks, IWANN 2019, 
% Gran Canaria, Spain, June 12-14, 2019 : proceedings. Pt. 2 / eds. Ignacio Rojas, Gonzali Joya, Andreu Catalia, 
% Cham : Springer, cop. 2019. pp 235-245 Lecture Notes in Computer Science, ISSN 0302-9743; vol. 11507)
% DOI: doi.org/10.1007/978-3-030-20518-8_20		
%
% Implemented by Krzysztof Fonal, Rafal Zdunek and Tomasz Sadowski 2019

Mp = permute(M,[1 3 2]); % original tensor
Omega_p = permute(Omega,[1 3 2]); % original tensor

ranks = [J J];
sigm2 = 1;

% setup test params
alg = 'svd';
tol = 1e-12;
k=1;

% Initialization
for i = 1:size(M,3)
    Ym(:,:,i) = imgaussfilt(M(:,:,i),10); % low-pass 2D filter
end
Ym = permute(Ym,[1 3 2]);

Ym(Omega_p) = Mp(Omega_p);
Mp_omega = Mp(Omega_p);
Mp_norm = norm(Mp_omega(:));
e(1) = 0;

while k < maxiter
    
k=k+1;

if k<2*J
    ranks = [10+k 10+k];
else
    ranks = [2*J 2*J];
end
Ym_tt = tt_decomposition(Ym, ranks, alg, 'tol', tol, 'maxiters', maxiter);
Ym = max(0,tt_to_tensor(Ym_tt));

for i = 1:size(M,3)
    Ym(:,i,:) = imgaussfilt(squeeze(Ym(:,i,:)),sigm2);
end

Rm = Mp(Omega_p) - Ym(Omega_p);
Ym(Omega_p) = Mp(Omega_p);

e(k) = norm(Rm(:))/Mp_norm;
disp(['tt_IC: ' , 'iteration: ',num2str(k), ' error: ',num2str(e(k))])

    % Stopping
    if (k > 2) && (e(k-1) - e(k)) < er
       break
    end

end

Y_hat=double(tt2img(permute(Ym,[1 3 2]),1));
    
%end
