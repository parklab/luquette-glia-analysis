%% Clearing workspaces
close all;
clearvars;
clc;
addpath('source');
addpath(['source' filesep 'text_import']);

%% Starting the default cluster (PLEASE USE AT LEAST 100 CORES)
if ( ~isempty(gcp('nocreate')) )
    delete(gcp);
end
c = parcluster; 
job = parpool(c);

%% Perform mutational signatures analysis
tic
    inputFolder  = 'input';
    inputFile    = '21_breast_WGS_substitutions'; % please do not use an extension 
    outputFolder = ['output' filesep '21_WGS_BRCA_try2'];
    analysisSignatures(inputFolder, outputFolder, inputFile, ... % input folder; output folder; input file
                       2, 2, 10, ... % start number of signatures; end number of signatures; total iterations per core (please make at least 1,000 iterations)
                       job); % parallel job variable
toc
