function run_one(nsigs, ctx, iters, inputPath, outputPath, sourcePath)
    %% Clearing workspaces
    close all;
    addpath(sourcePath);
    addpath([sourcePath filesep 'text_import']);
    
    % Prevents automatic pool creation on spmdn deconvolute.m.
    % Parallel pooling is for some reason completely broken on O2.
    ps = parallel.Settings;
    ps.Pool.AutoCreate = false;
    % Fool decipherMutationalProcesses into thinking a pool was created.
    job.NumWorkers = 1;

    %% Perform mutational signatures analysis
    tic
        inputFolder  = inputPath; %strcat('nsigs', int2str(nsigs), '/SBS96/');
        inputFile    = 'input.csv';
        outputFolder = outputPath; %strcat('nsigs', int2str(nsigs), '/SBS96/');
        analysisSignaturesTextIO(inputFolder, outputFolder, sourcePath, inputFile, ctx, ...
                        nsigs, nsigs, iters, ... % start; end; iterations (per core?)
                        job); % parallel job variable
    toc
