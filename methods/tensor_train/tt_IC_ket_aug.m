function  [Y_hat]  = tt_IC_ket_aug(M, J, Omega, maxiter, er)
%Image competion using tensor train with Ket augmentation
%function  [Y_hat]  = tt_IC_ket_aug(M, J, Omega, maxiter, er)
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

if size(M,1)==256
    d = 16; % window size
elseif size(M,1)==512
    d=32;
else
    error('Define other window size ')
    
end
sigma2 = 1; % sigma

% setup test params
alg = 'svd';
tol = 1e-12;
k=1;

% Initialization
Mm=img2tt(M ,d); % multi-mode tensor M
for i = 1:3
    Ym(:,:,i) = imgaussfilt(M(:,:,i),10);
end
Ym=img2tt(Ym ,d); % multi-mode tensor M
Omega_m=logical(img2tt(double(Omega),d)); % multi-mode tensor Omega
Ym(Omega_m) = Mm(Omega_m);
Mm_omega = Mm(Omega_m);
Mm_norm = norm(Mm_omega(:));
e(1) = 0;

while k < maxiter
    
k=k+1;
j=k;
if k>2*J
    j=2*J;
end
ranks = [10 5+j 20 3]; % d = 32
%ranks = [10 30+j 3]; % d = 64
%ranks = [10 20+j 100 10 3]; % d = 16
%ranks = [10 100 47 2 3]; % 256

% TT decomposition
Ym_tt = tt_decomposition(Ym, ranks, alg, 'tol', tol, 'maxiters', maxiter);
Ym = max(0,tt_to_tensor(Ym_tt));

% Filtering
Lm = size(Ym);
Ym_12 = reshape(Ym,[Lm(1) Lm(2) prod(Lm)/(Lm(1)*Lm(2))]);
for i = 1:size(Ym_12,3)
    Ym_12(:,:,i) = imgaussfilt(squeeze(Ym_12(:,:,i)),sigma2);
end
Ym = reshape(Ym_12,Lm);

Rm = Mm(Omega_m) - Ym(Omega_m);
Ym(Omega_m) = Mm(Omega_m);

e(k) = norm(Rm(:))/Mm_norm;
disp(['tt_IC: ' , 'iteration: ',num2str(k), ' error: ',num2str(e(k))])

    % Stopping
    if (k > 2) && (e(k-1) - e(k)) < er
       break
       %disp('yyy')
    end

end

Y_hat=double(tt2img(Ym,1));
    
%end
