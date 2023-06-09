%% Example for running SigProfilerSingleSample 
%% The example examines 10 biliary adenocarcinoma whole-genomes from ICGC PCAWG and
%% assigns mutational signatures using upcoming release of PCAWG mutational signatures
%% Clearing all data
close all;
clearvars;
clc;
addpath('source');
addpath(['source' filesep 'text_import']);

%% Starting the default cluster
if ( ~isempty(gcp('nocreate')) )
    delete(gcp);
end
c = parcluster;
job = parpool(c);

%% Analysis of all signatures in individual samples
tic
    analysis_individual_samples(['input' filesep 'common' filesep 'Consensus_subs_mutational_signatures.mat'], ... % set of signatures
                                ['input' filesep 'common' filesep 'signatures_in_samples_and_cancer_types.mat'], ... % set of signatures in samples
                                ['input' filesep 'matlab' filesep  'Biliary-AdenoCA_example_cancer_samples.mat'], ... % set of individual samples for examination
                                ['output' filesep 'Biliary-AdenoCA'], ... % output folder
                                'Biliary-AdenoCA', ... % output file
                                1, ... % boolean variable indicating whether to use rules or not (1==use rules; 0==do not use rules)
                                [1 5], ... % IDs of signatures to be included in all samples regardless of rules or sparsity  
                                {[2 17], [7 8 9 10], [13 14], [21 22]}); % connected signatures (e.g., signatures SBS-2 and SBS-13)
toc