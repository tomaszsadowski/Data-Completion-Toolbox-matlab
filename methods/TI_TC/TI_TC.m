function [Y,res]= TI_TC(M,Omega,tau,z,func,H)
%Tensor Interpolation for Image Completion
%function [Y,res]= Biharmonic_fitting_block_param(M,Omega,tau,fun)
%Inputs:
% M - incomplite data (image), 
% Omega - a binary tensor of indicators for the missing entries,  
% tau - parameter =3 optimal for 'chebychev' i minkovsky P = 0.5
% z - interpolation function, switch cases 
%     6  :  mixed functions
%     5  :  only polynomial
%     4  :  only exp
% func - parameter, switch cases
%      'exp' - exponential finction
%      'cheb'- chebyshev function
%      'mink' - mikowski function
% H - parameter, H=3 : 3rd degree polynolmial, H=2 else : 2nd degree polynolmial
%
% Outputs:
% Y_hat - completed data (image)
%
% Requires: MATLAB Tensor Toolbox, Sandia Corporation
% (https://www.tensortoolbox.org)
%
% Cite:
% Rafal Zdunek and Tomasz Sadowski,
% Image completion with hybrid interpolation in tensor representation,
% Applied Sciences. 2020, vol. 10, num. 3, art. 797, pp 1-17.
% DOI: doi.org/10.3390/app10030797	
%
% Implemented by Rafal Zdunek and Tomasz Sadowski 2020


res = [];

W(:,:,1:3) = M;
W(:,:,4:6) = double(Omega);

% change depending on size of data, 
DimB = [26 26 3];

%tau = 3; % optimal for 'chebychev' i minkovsky P = 0.5

B = 1; C = 1;
F = cell([1 length(DimB)]);
for l = 1:length(DimB)
    inx = (1:DimB(l))';
    
    switch func
        
        case 'exp'
            F{l} = 10.^(-pdist2(inx,inx).^2/tau);
        case 'cheb'
            F{l} = exp(-(pdist2(inx,inx,'chebychev'))/tau);
        case 'mink'
            F{l} = exp(-(pdist2(inx,inx,'minkowski',.5))/tau);
        otherwise
            errror('Undefined function "func", check help')
    end
    if H==3
        P = [ones(length(inx),1) inx-mean(inx) (inx-mean(inx)).^2 (inx-mean(inx)).^3];
    else
        P = [ones(length(inx),1) inx inx.^2];
    end
    B= kron(F{l},B);
    C = kron(P,C);
end

switch z
       
    case 6
        fun = @(block_struct) Biharmonic_fun_6(block_struct,B,C); % mixed
    case 5
        fun = @(block_struct) Biharmonic_fun_5(block_struct,C); % only polynomial
    case 4
        fun = @(block_struct) Biharmonic_fun_4(block_struct,B,F); % only exp

    otherwise
        error('Indefined function "z", check help')
end

%parallel - true or false
Y = blockproc(W,[16 16],fun,'BorderSize',[5 5]);%,'UseParallel',true);


end


