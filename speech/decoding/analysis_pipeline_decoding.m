%--------------------------------------------------------------------------
% Till Habersetzer, 05.07.2025
% Communication Acoustics, CvO University Oldenburg
% till.habersetzer@uol.de 
%
% Description:
%   Master script for running the complete MEG-audio neural decoding
%   analysis pipeline (without plotting). It orchestrates the preprocessing
%   of audio and MEG data and computes subject-level decoding models in
%   parallel.
%
%   Workflow:
%   - Initializes settings, paths, and starts a log file.
%   - Preprocesses the shared audio stimuli (OLSA and audiobooks).
%   - Processes each subject in parallel (parfor):
%       1. Preprocesses the subject's MEG data for both tasks.
%       2. Trains and evaluates the subject-level decoding model.
%   - Logs the entire process and total execution time.
%
% To run from the command line:
% matlab2024b -nodisplay -nosplash -r "analysis_pipeline_decoding; exit;"
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
fname_log = 'analysis_pipeline_decoding_log.txt';

if exist(fname_log, 'file') == 2
    delete(fname_log);
    fprintf('Existing log file "%s" deleted.\n', fname_log);
end

% Start recording output and timer
diary(fname_log);  

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
%--------------------------------
preprocessing_olsa(settings)
% Preprocessing of audiobook envelopes
%-------------------------------------
preprocessing_audiobooks(settings)

parpool(workers_to_use);
parfor sub_idx = 1:n_subj
    % Select subject
    %---------------
    subject = sprintf('sub-%02d',subjects(sub_idx));

    % Propecessing of audiobook MEG recordings 
    %-----------------------------------------
    preprocessing_audiobooks_decoding(subject,settings)

    % Propecessing of olsa MEG recordings 
    %------------------------------------
    preprocessing_olsa_decoding(subject,settings)

    % Training decoder
    %-----------------
    training_decoding(subject,settings)

end % subjects

% Stop timer and log the elapsed time
elapsed_time = toc;
fprintf('\n--- Pipeline Finished ---\n');
fprintf('Total processing time: %.2f minutes.\n', elapsed_time / 60);

diary off; % Stop logging output

% Use 'plot_decoding.m' to visualize the results separetely

