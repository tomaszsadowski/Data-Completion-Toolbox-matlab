function [ A,B ] = lrmf( Y, alg, rank, varargin )
%LRMF Summary of this function goes here
%   Detailed explanation goes here

%% Set algorithm parameters from input or by using defaults
params = inputParser;
params.addParamValue('tol',1e-4,@isscalar);
params.addParamValue('maxiters',50,@(x) isscalar(x) & x > 0);
%params.addParamValue('init', 'random', @(x) (iscell(x) || ismember(x,{'random','nvecs'})));
params.addParamValue('printitn',1,@isscalar);
params.parse(varargin{1}{:});

%%
if strcmp(alg, 'svd')
    [U,S,V] = svd(Y, 'econ');
    V=V*S;
    A = U(:,1:rank);
    B = V(:,1:rank)';
else
    error('Add new method')
end
    
end

