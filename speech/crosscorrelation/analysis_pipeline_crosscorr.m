%--------------------------------------------------------------------------
% Till Habersetzer, 05.07.2025
% Communication Acoustics, CvO University Oldenburg
% till.habersetzer@uol.de 
%
% Description:
%   Master script for running the 'complete' MEG-audio cross-correlation
%   analysis pipeline (without plotting).
%   It orchestrates the preprocessing of audio and MEG data, computes 
%   subject-level results in parallel, and finishes with a group-level 
%   grand average.
%
%   Workflow:
%   - Initializes settings, paths, and starts a log file.
%   - Preprocesses the shared audio stimuli.
%   - Processes each subject in parallel (parfor):
%       1. Preprocesses the subject's MEG data.
%       2. Computes the subject-level cross-correlation.
%   - Computes the group-level grand average across all subjects.
%   - Logs the entire process and total execution time.
%
% To run from the command line:
% matlab -nodisplay -nosplash -r "analysis_pipeline_crosscorr; exit;"
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
diary('analysis_pipeline_crosscorr_log.txt');  
tic; 
fprintf('--- Starting Cross-Correlation Pipeline ---\n');
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

%% Analysis
%--------------------------------------------------------------------------
max_workers    = parcluster('local').NumWorkers; % Get the maximum number of workers available on your machine
workers_to_use = min(n_subj, max_workers);       % Determine the number of workers to use (the smaller of the two)

% Preprocessing of audio envelopes
preprocessing_audio(settings)

parpool(workers_to_use);
parfor sub_idx = 1:n_subj
    % Select subject
    subject = sprintf('sub-%02d',subjects(sub_idx));
 
    % Propecessing of audiobook MEG recordings 
    preprocessing_crosscorr(subject,settings)
 
    % Computation cross-correlation function
    computation_crosscorr(subject,settings)

end % subjects

% Computation grand-average cross-correlation function
computation_gavg_crosscorr(subjects,settings)

% Stop timer and log the elapsed time
elapsed_time = toc;
fprintf('\n--- Pipeline Finished ---\n');
fprintf('Total processing time: %.2f minutes.\n', elapsed_time / 60);

diary off; % Stop logging output

% Use 'plot_crosscorr.m' to visualize the results separetely

