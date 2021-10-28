%Data completion toolbox 

clear
close all

%select methods: 
% 1 - SmNMF-MC, 2 - HALS, 3- BSA-IC, 4- ICSA 
% 5 - Tucker ALS, 6- PTT, 7 - KA-TT, 8 - fALS, 
% 9 - Tensor Intrepolation, 10-12- GDC,

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%prepare data path, adding functions 
DataFolder = pwd;
addpath(genpath(fullfile(DataFolder,'methods')));
addpath(genpath(fullfile(DataFolder,'datasets')));

%choose benchmark
%ben - select image or other kind of dataset, avaiable datasets: 'lena', 'barbara', 'monarch'
%type - select type of distortion ( numeric ratio - holes (how many lost), grid, cell ratio - noise (how many lost).... )

%main function, for more info check 'help tests_toolbox'

[Y,SIR,time_elapsed,Structure_data]=tests_toolbox('lena', 1:12, 0.5);

%show results
[SIR_table,~,time_table]=results_print(Y,SIR,time_elapsed,Structure_data);

eval('SIR_table')
eval('time_table')
% outputFile = ['results.mat'];
% save(outputFile)
disp('Fin without errors')
