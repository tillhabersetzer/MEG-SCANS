function computation_crosscorr(subject,settings)
%--------------------------------------------------------------------------
% Till Habersetzer, 05.07.2025
% Communication Acoustics, CvO University Oldenburg
% till.habersetzer@uol.de 
%
% Description:
%   Performs a subject-level temporal cross-correlation analysis between
%   MEG channel data and a corresponding audio envelope.
%
%   Key Steps:
%   - Loads preprocessed and synchronized MEG and audio data for a subject.
%   - For each trial, computes the cross-correlation between every MEG
%     channel and the aligned audio envelope.
%   - Repeats the analysis using randomly shuffled MEG-audio trial pairings
%     to create a null distribution for statistical comparison.
%   - Averages the results across all trials for both the actual and
%     shuffled conditions.
%   - Saves the final averaged cross-correlation data to a .mat file.
%
% Inputs:
%   subject (char/string): The subject identifier (e.g., 'sub-01').
%   settings (struct):     Settings structure with all necessary paths 
%                          and parameters.
%--------------------------------------------------------------------------

%% Script settings 
%--------------------------------------------------------------------------

% Apply additional zscoring of all trials
apply_zscore = settings.crosscorr.zscore;

% Downsampling frequency
fs_down = settings.crosscorr.fs_down;

% desired lags for crosscorrelation
timelags = settings.crosscorr.timelags;

%% Import data
%----------------------------------------------------------------------
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
%-----------------------------------------------------------------------
if apply_zscore
    for trl_idx = 1:n_trials
        epochs_audio{trl_idx}       = zscore(epochs_audio{trl_idx},0,2);
        epochs_neuro.trial{trl_idx} = zscore(epochs_neuro.trial{trl_idx},0,2);
    end
end

%% Compute Crosscorrelation over channels and trials
%-----------------------------------------------------------------------
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
cfg                    = [];
cfg.channel            = 'all';
avg_crosscorr          = ft_timelockanalysis(cfg,epochs_crosscorr);
avg_crosscorr_shuffled = ft_timelockanalysis(cfg,epochs_crosscorr_shuffled);
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
results.avg_crosscorr          = avg_crosscorr;
results.avg_crosscorr_shuffled = avg_crosscorr_shuffled;

save(fullfile(dir2save,fname),'results','-v7.3'); 
fprintf("\n%s from %s saved.\n",fname,subject)

end % end of function