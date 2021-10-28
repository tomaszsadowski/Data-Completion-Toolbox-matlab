function [Y,SIR,time_elapsed,Structure_data]=tests_toolbox(data, Methods,bench,varargin)
%function [Y,SIR,time_elapsed,procent]=tests_toolbox(data, Methods,bench);
%Inputs:
% data - select image or other kind of dataset
%   avaiable datasets: 'barbara','lena', 'monarch'
% Methods - choose method: 
%   1 - SmNMF-MC, 2 - HALS, 3- BSA-IC, 4- ICSA 
%   5 - Tucker ALS, 6- PTT, 7 - KA-TT, 8 - fALS, 
%   9 - Tensor Intrepolation, 10-12- GDC, 
% bench - choose benchmark type - select type of distortion:
%   * numeric ratio - holes (how many lost), example: 0.5,
%   * grid distortion: 'grid' ,
%   * cell ratio - noise (how many lost): {0.5}
%
%Outputs:
% Y - cell data with restored data,
% SIR - Signal to Interferance Ratio, efficiency measurment 
% time_elapsed - times of calculations for all methods and MC runs
% Structure_data - structure with all data, used in result_print.m
%
%After bench you can use additional parameters: 'r' - rank, 'tol' -
%tolerance, 'MC' - Monte Carlo runs, 'maxiter' - max number of iterations, etc.
%Example: run test for image lena, testing method HALS and ICSA, for 15% pixels avaiable, with rank = 30, 
%tolerance = 10^-5 for 10 Monte Carlo runs and max iteration = 100:
%[Y,SIR,time_elapsed,Structure_data]=tests_toolbox('lena', [2,4], 0.5, 'r', 30, 'tol', 10e-5, 'MC',10, 'maxiter',100);
%
%To use Matlab Parallel Pools, add option: 'para',1 .
%
%To save results in mat file use option: 'saving',1 , to name it after that
%add option 'name', 'name you want'
%
%To see results, evaluate:
% [SIR_table,~,time_table]=results_print(Y,SIR,time_elapsed,Structure_data)
%
%Implemented by Tomasz Sadowski 2016-2021

% Set algorithm parameters from input or by using defaults
params = inputParser;
params.addParameter('r',50,@isscalar);
params.addParameter('tol',1e-12,@isscalar);
params.addParameter('maxiter',100,@(x) isscalar(x) & x > 0);
params.addParameter('MC',1,@isscalar);
params.addParameter('show_inx',1,@isscalar);
params.addParameter('name','none',@ischar);
params.addParameter('saving',0,@isscalar);
params.addParameter('para',0,@isscalar);
params.addParameter('tau',3,@isscalar);
params.addParameter('H',2,@isscalar);
%params.addParameter('res_p',0,@isscalar);

%params.addParameter('init', 'random', @(x) (iscell(x) || ismember(x,{'random','nvecs','eigs'})));
params.parse(varargin{:});

% Copy from params object
r=params.Results.r;
tol=params.Results.tol;
maxiter=params.Results.maxiter;
MC=params.Results.MC;
show_inx=params.Results.show_inx;
para=params.Results.para;
saving=params.Results.saving;
name=params.Results.name;

H=params.Results.H;
tau=params.Results.tau;
%res_p=params.Results.res_p;


% Parallel pools
if para==1
    poolobj = gcp('nocreate');
    if isempty(poolobj)
        parpool('local');
    end
end

%select methods

Structure_data.Methods=Methods;
labels_methods_NMF = {'SmNMF-MC','HALS','BSA-IC','ICSA','tucker-als','PTT', 'KA-TT','fALS','TI-TC' ...
     'GDC(H)', 'GDC(V)', 'GDC(HV)'};
Structure_data.labels=labels_methods_NMF;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%prepare data path, adding functions
DataFolder = pwd;
addpath(genpath(fullfile(DataFolder,'methods')));
addpath(genpath(fullfile(DataFolder,'datasets')));


%prepare benchamrk
[X,ohm,X0,procent]= choose_benchmark(data, bench ,DataFolder);

%Specify Parameters

Structure_data.r=r;
Structure_data.show_inx=show_inx;
Structure_data.maxiter=maxiter;
Structure_data.tol= tol;
Structure_data.MC=MC;
Structure_data.X=X;
Structure_data.X0=X0;
Structure_data.ohm=ohm;

Structure_data.tau=tau;
Structure_data.H=H;

clear X X0 ohm;
% perform test with specified parameters

if para==1
    %parrarel version
    [Y,SIR,time_elapsed,Structure_data] = IC_methods_p(Structure_data,Methods);
else
    [Y,SIR,time_elapsed,Structure_data] = IC_methods(Structure_data,Methods);
end


if saving==1
    

    if strcmp(name,'none')
        
        outputFile = ['results/',data,'_MCruns_',num2str(Structure_data.MC),procent,'.mat'];
        save(outputFile)
        
    else
        outputFile = ['results/',name,'.mat'];
        save(outputFile)
        
    end
    
end
disp('Fin without errors')
end