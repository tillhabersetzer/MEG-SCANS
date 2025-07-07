%--------------------------------------------------------------------------
% Till Habersetzer, 05.07.2025
% Communication Acoustics, CvO University Oldenburg
% till.habersetzer@uol.de 
%
% To run from the command line:
% matlab -nodisplay -nosplash -r "analysis_pipeline_decoding; exit;"
%--------------------------------------------------------------------------

close all
clearvars
clc 

%% Import main settings 
%--------------------------------------------------------------------------
current_dir = pwd;
cd(fullfile('..'))
settings_speech
cd(current_dir)

%% Script settings 
%--------------------------------------------------------------------------
% Start recording output and timer
diary('analysis_pipeline_decoding_log.txt');  
tic; 
fprintf('--- Starting Decoding Pipeline ---\n');
fprintf('Date: %s\n\n', datetime('now', 'Format', 'dd-MMM-yyyy HH:mm:ss'));

% Select subjects
subjects                    = 1:24;
n_subj                      = length(subjects);
settings.crosscorr.subjects = subjects; % Copy into settings

% Add paths
%----------
addpath(fullfile(settings.path2project,'analysis','helper_functions'))

% Addpath for FieldTrip - is only added for Server
if contains(settings.rootpath,'/mnt/localSSDPOOL')
    addpath(settings.path2fieldtrip)
    ft_defaults
end

% run script for auditory toolbox
run(fullfile(settings.path2amtoolbox,'amtstart'));

% Addpath for mTRF Toolbox
addpath(genpath(settings.path2mtrftoolbox));

%% Analysis
%--------------------------------------------------------------------------
max_workers    = parcluster('local').NumWorkers; % Get the maximum number of workers available on your machine
workers_to_use = min(n_subj, max_workers);       % Determine the number of workers to use (the smaller of the two)

% Preprocessing of olsa envelopes
% preprocessing_olsa(settings)
% Preprocessing of audiobook envelopes
% preprocessing_audiobooks(settings)

parpool(workers_to_use);
parfor sub_idx = 1:n_subj
    % Select subject
    subject = sprintf('sub-%02d',subjects(sub_idx));

    % Propecessing of audiobook MEG recordings 
    % preprocessing_audiobooks_decoding(subject,settings)

    % Propecessing of olsa MEG recordings 
    % preprocessing_olsa_decoding(subject,settings)

    % Training decoder
    training_decoding(subject,settings)

end % subjects

% Computation grand-average cross-correlation function
% computation_gavg_crosscorr(subjects,settings)

% Stop timer and log the elapsed time
elapsed_time = toc;
fprintf('\n--- Pipeline Finished ---\n');
fprintf('Total processing time: %.2f minutes.\n', elapsed_time / 60);

diary off; % Stop logging output

% Use 'plot_crosscorr.m' to visualize the results separetely

