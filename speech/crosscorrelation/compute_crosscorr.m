%--------------------------------------------------------------------------
% Till Habersetzer, 26.06.2025
% Communication Acoustics, CvO University Oldenburg
% till.habersetzer@uol.de 
%
%   Description:
%   This script calculates the temporal cross-correlation between preprocessed
%   MEG data and the corresponding audio envelope. It performs this
%   computation for individual subjects and for a group-level grand average.
%
%   The script executes the following workflow:
%   1.  Iterates through each subject, loading their preprocessed data.
%   2.  For each trial, computes the cross-correlation between every MEG
%       channel and the aligned audio envelope.
%   3.  As a control, it repeats the analysis using randomly shuffled
%       pairings of MEG and audio trials to create a null distribution.
%   4.  Averages the cross-correlations across trials for each subject and
%       saves the individual results.
%   5.  After processing all subjects, it computes and saves a grand
%       average across the entire group for both the actual and shuffled
%       conditions.
%
% To run from the command line (linux server):
% matlab -nodisplay -nosplash -r "compute_crosscorr; exit;"
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
% For test purposes of eeg reference, only select single session and
% condition

% Select seubjects
subjects = 1:24; % all files available for sub-02-sub-23
% subjects = 4;
n_subj   = length(subjects);

% Apply additional zscoring of all trials
apply_zscore = settings.crosscorr.zscore;

% Downsampling frequency
fs_down = settings.crosscorr.fs_down;

% desired lags for crosscorrelation
timelags = settings.crosscorr.timelags;

% Add fieldtrip
addpath(settings.path2fieldtrip)
ft_defaults

%% Import data for each subject and compute crosscorrelation
%--------------------------------------------------------------------------
subjectnames = cell(1,n_subj);
n_trials_all = zeros(1,n_subj);

avgs_crosscorr          = cell(1,n_subj);
avgs_crosscorr_shuffled = cell(1,n_subj);

for sub_idx = 1:n_subj

    %% Import data
    %----------------------------------------------------------------------
    subject               = sprintf('sub-%02d',subjects(sub_idx));
    subjectnames{sub_idx} = subject;

    data         = importdata(fullfile(settings.path2project,'derivatives',subject,'speech',sprintf('%s_preprocessed_crosscorr.mat',subject)));   
    epochs_audio = data.epochs_audio;
    epochs_neuro = data.epochs_neuro;
    n_trials     = length(epochs_neuro.trial);

    if ~isequal(epochs_neuro.fsample,fs_down)
        error('%s: Unexpected sampling frequency %i Hz!',subject,epochs_neuro.fsample)
    end

    fprintf('%s loaded.\n',subject)
    clear data

    %% Apply zscoring
    %----------------------------------------------------------------------
    if apply_zscore
        for trl_idx = 1:n_trials
            epochs_audio{trl_idx}       = zscore(epochs_audio{trl_idx},0,2);
            epochs_neuro.trial{trl_idx} = zscore(epochs_neuro.trial{trl_idx},0,2);
        end
    end

    %% Compute Crosscorrelation over channels and trials
    %----------------------------------------------------------------------
    n_chan = length(epochs_neuro.label);

    % get correct index for lag vector
    timelags_samples = round(timelags/1000*fs_down); % in samples
    [~,lags_samples] = xcorr(epochs_neuro.trial{1}(1,:),epochs_audio{1},timelags_samples(2),'coeff'); % limits the lag range from -maxlag to maxlag

    lags_idx = dsearchn(lags_samples',timelags_samples');
    lags_sec = lags_samples(lags_idx(1):lags_idx(2))*(1/fs_down); % timelags in seconds
    n_lags   = length(lags_sec);

    % Create random permutation of trials
    %------------------------------------
    rng('shuffle')
    % rng(30062025);
    idx_shuffled = randperm(n_trials); % for random mapping of audiodata
    % check that elements aren't identical, so the difference should never be 0
    while ~all(idx_shuffled-(1:n_trials)) % check for nonzero elements
        rng('shuffle')
        idx_shuffled = randperm(n_trials);
        fprintf('Shuffled again.\n')
    end
    
    % Initialize crosscorr data
    epochs_crosscorr               = epochs_neuro;
    epochs_crosscorr.time          = repmat({lags_sec},1,n_trials);
    epochs_crosscorr_shuffled      = epochs_neuro;
    epochs_crosscorr_shuffled.time = repmat({lags_sec},1,n_trials);
    
    % Loop over trials and channels
    for trl_idx = 1:n_trials
        cfc          = zeros(n_chan,n_lags);
        cfc_shuffled = zeros(size(cfc));
        for ch_idx = 1:n_chan
            % xcorr(x(n+m),y(n),maxlag), x arrives latter for m>0 (so neuro and postiv m)
            [r,~]          = xcorr(epochs_neuro.trial{trl_idx}(ch_idx,:),epochs_audio{trl_idx},timelags_samples(2),'coeff');
            [r_shuffled,~] = xcorr(epochs_neuro.trial{trl_idx}(ch_idx,:),epochs_audio{idx_shuffled(trl_idx)},timelags_samples(2),'coeff');
     
            cfc(ch_idx,:)          = r(lags_idx(1):lags_idx(2));
            cfc_shuffled(ch_idx,:) = r_shuffled(lags_idx(1):lags_idx(2));
            % fprintf('trial: %i | channel: %i processed.\n',trl_idx,ch_idx);
        end
        epochs_crosscorr.trial{trl_idx}          = cfc;
        epochs_crosscorr_shuffled.trial{trl_idx} = cfc_shuffled;   
        clear cfc cfc_shuffled
    end

    % Compute average 
    %----------------
    cfg                              = [];
    cfg.channel                      = 'all';
    avgs_crosscorr{sub_idx}          = ft_timelockanalysis(cfg,epochs_crosscorr);
    avgs_crosscorr_shuffled{sub_idx} = ft_timelockanalysis(cfg,epochs_crosscorr_shuffled);
    clear epochs_crosscorr epochs_crosscorr_shuffled

    %% Save results
    %--------------
    dir2save = fullfile(settings.path2derivatives,subject,'speech');
    if ~exist(dir2save,'dir')
        mkdir(dir2save)
    end
    fname = sprintf('%s_crosscorr.mat',subject);
    
    results                        = struct;
    results.settings               = settings;
    results.n_trials               = n_trials;
    results.avg_crosscorr          = avgs_crosscorr{sub_idx};
    results.avg_crosscorr_shuffled = avgs_crosscorr_shuffled{sub_idx};
    
    save(fullfile(dir2save,fname),'results','-v7.3'); 
    fprintf("\n%s from %s saved.\n",fname,subject)

end % Loop over subjects

%% Add grand average over all subjects
%--------------------------------------------------------------------------
subject = 'grandaverage';

cfg                     = [];
cfg.latency             = 'all';
gavg_crosscorr          = ft_timelockgrandaverage(cfg,avgs_crosscorr{:});
gavg_crosscorr_shuffled = ft_timelockgrandaverage(cfg,avgs_crosscorr_shuffled{:});

dir2save = fullfile(settings.path2derivatives,subject,'speech');
if ~exist(dir2save,'dir')
    mkdir(dir2save)
end
fname    = sprintf('%s_crosscorr.mat',subject);

results                         = struct;
results.n_subjects              = n_subj;
results.gavg_crosscorr          = gavg_crosscorr;
results.gavg_crosscorr_shuffled = gavg_crosscorr_shuffled;

save(fullfile(dir2save,fname),'results','-v7.3'); 
fprintf("\n%s from %s saved.\n",fname,subject)