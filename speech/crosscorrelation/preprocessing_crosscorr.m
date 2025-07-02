%--------------------------------------------------------------------------
% Till Habersetzer, 26.06.2025
% Communication Acoustics, CvO University Oldenburg
% till.habersetzer@uol.de 
%
%   Description:
%   This script preprocesses magnetoencephalography (MEG) and corresponding
%   audiobook data to prepare it for a cross-correlation analysis.
%
%   The script performs the following main steps for each subject:
%   1.  Loads continuous MEG data and stimulus audio files.
%   2.  Filters MEG data (band-pass) and computes the audio envelope.
%   3.  Defines and creates epochs of a fixed duration for both modalities
%       based on experimental triggers.
%   4.  Downsamples the epoched MEG and audio data to a common frequency.
%   5.  Synchronizes the trial counts between MEG and audio, discarding
%       invalid or mismatched epochs.
%   6.  Appends data from all experimental runs.
%   7.  Saves the final preprocessed data structures ('epochs_neuro',
%       'epochs_audio') into a single .mat file per subject.
%
% To run from the command line (linux server):
% matlab -nodisplay -nosplash -r "preprocessing_crosscorr; exit;"
%--------------------------------------------------------------------------

close all
clearvars
clc 

% To Do
% include other stories for sub-01

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
% subjects = 3;
n_subj   = length(subjects);

% Filter settings
bpfreq = settings.crosscorr.meg.bpfreq; % meg
lpfreq = settings.crosscorr.audio.lpfreq.audio; % audio

% Epoch length
trialdur = settings.crosscorr.trialdur;

% Downsampling frequency
fs_down = settings.crosscorr.fs_down;

% Audio frequency
fs_audio = settings.fs_audio;

% check whether maxfiltered data should be analyzed
use_maxfilter = settings.use_maxfilter;

% Latency Correction
apply_latency_correction = settings.apply_latency_correction;
audio_latency            = settings.audio_latency;

bids_dir = settings.path2bids;
stim_dir = fullfile(bids_dir,'stimuli','audiobooks'); % for audio data

% Add fieldtrip
addpath(settings.path2fieldtrip)
ft_defaults

% Addpath for additional functions
addpath(fullfile(settings.path2project,'analysis','helper_functions'))

%% Processing
%--------------------------------------------------------------------------

% Init audiobook filter
%----------------------
[b,a] = butter(3,lpfreq/(fs_audio/2),'low');
% h     = fvtool(b,a);
% h.Fs  = fsamp;

% Loop over subjects, both audiobooks, including both runs
for sub_idx = 1:n_subj

    subject = sprintf('sub-%02d',subjects(sub_idx));
    
    % Initialize data
    epochs_neuro   = cell(1,4);
    epochs_audio   = cell(1,4);
    mapping_epochs = zeros(4,2); % mapping between audiobooks and epoch trials (audiobooks x [begin end])
    mapping_labels = cell(1,4); % contains audiobook labels

    epoch_counter  = 0;

    f_idx = 1; % file counter
    for task_idx = 1:2     
        task = sprintf('audiobook%i',task_idx);
        for run_idx = 1:2
            run = sprintf('0%i',run_idx);

            % Define datapath
            %----------------
            if use_maxfilter
                rawdata_path = fullfile(settings.path2derivatives,subject,'maxfilter');
                % Use movement-corrected data first
                filename     = sprintf('%s_task-%s_run-%s_proc-tsss-mc_meg.fif',subject,task,run); 
                data_path    = fullfile(rawdata_path,filename);
        
                % use bids data for events
                event_path = fullfile(settings.path2bids,subject,'meg',sprintf('%s_task-%s_run-%s_meg.fif',subject,task,run));
            else
                rawdata_path = fullfile(settings.path2bids,subject,'meg');
                filename     = sprintf('%s_task-%s_run-%s_meg.fif',subject,task,run);
                data_path    = fullfile(rawdata_path,filename);
                % badchannels  = settings.badchannels; 
                event_path   = data_path;
            end

            % Process data if it exists
            if isfile(data_path)

                %% Define trials
                %------------% based on bids events
                cfg                     = [];
                cfg.trialfun            = 'my_trialfun_audiobook'; 
                cfg.dataset             = event_path; % path should always point to original bids-path
                cfg.trialdef.method     = 'fieldtrip';
                cfg.trialdef.trialdur   = 10; % 10 (s)
                cfg.trialdef.eventvalue = 1;      
                % bids
                % cfg.trialdef.eventtype  = {'audiobook 1s'};
                % fieldtrip
                cfg.trialdef.eventtype  = 'STI101';   
                cfg                     = ft_definetrial(cfg);

                trl = cfg.trl;

                %% Add latency correction 
                %------------------------
                if apply_latency_correction
                    trl(:,1:2) = trl(:,1:2) + round(audio_latency*1000); % in samples
                end

                %% Filter continuous data to avoid edge artifacts
                %------------------------------------------------
                cfg            = [];
                cfg.dataset    = data_path;
                cfg.bpfilter   = 'yes';
                cfg.bpfreq     = bpfreq;
                cfg.channel    = 'meg';
                data           = ft_preprocessing(cfg); 

                % plot spectrum
                %--------------
                % cfg         = [];
                % cfg.output  = 'pow';
                % cfg.channel = 'meg';
                % cfg.method  = 'mtmfft';
                % cfg.taper   = 'hanning';
                % cfg.foi     = 0:100;    
                % spectr      = ft_freqanalysis(cfg,data);   
                % 
                % figure
                % plot(spectr.freq,db(spectr.powspctrm(:,:)))
                % xlabel('f / Hz')    
                % xlim([0,80])

                %% Epoch data
                %------------
                cfg     = [];
                cfg.trl = trl;  
                epoch   = ft_redefinetrial(cfg,data);
                % clear data

                % Set new offset for time: time 0 is the first trigger onset and
                % from there on the time moves forward in each trial. It 
                % doesn't start always from 0.
                cfg        = [];
                cfg.offset = trl(:,1)-trl(1,1);
                epoch      = ft_redefinetrial(cfg, epoch);   

                %% Downsampling
                %--------------
                cfg                 = [];
                cfg.resamplefs      = fs_down;
                cfg.detrend         = 'no';
                epochs_neuro{f_idx} = ft_resampledata(cfg,epoch);
                % clear epoch

                %% Preprocess audio data
                %-----------------------

                % Count epochs
                n_epochs                = length(epoch.trial);
                mapping_epochs(f_idx,:) = [epoch_counter+1,epoch_counter + n_epochs];
                epoch_counter           = epoch_counter + n_epochs;
                mapping_labels{f_idx}   = sprintf('task-%s_run-%s',task,run);

                %% Import raw audio and compute envelope
                %---------------------------------------
                fname_audio = sprintf('task-%s_run-%s_stim.wav',task,run);

                if strcmp(subject,'sub-01') % pilot had audiobooks
                    [raw_audiodata,fs] = audioread(fullfile(stim_dir,'pilot',fname_audio)); 
                else
                    [raw_audiodata,fs] = audioread(fullfile(stim_dir,fname_audio)); 
                end
                
                % Compute onset envelope
                envelope = cal_envelope(raw_audiodata,fs,b,a,'onset_envelope');

                % Cut into epochs and downsampling
                cfg             = [];
                cfg.fs_audio    = fs_audio;
                cfg.trialdur    = trialdur;
                cfg.fs_down     = fs_down;
                envelope_epochs = epoch_audiobook_crosscorr(cfg, envelope);
                fprintf('%s preprocessed.\n',fname_audio)

                epochs_audio{f_idx} = envelope_epochs;
                clear envelope_epochs raw_audiodata envelope

                % Increment counter
                f_idx = f_idx + 1;
                

            end % file exists
        end % loop over runs
    end % loop over tasks

    %% Append data over trials and correct for differences in number of trials and trial durations
    %---------------------------------------------------------------------------------------------
    % Delete empty cells if data is missing - otherwise appenddata won't
    % work
    nonEmptyIdx    = ~cellfun('isempty',epochs_neuro);
    % Keep only non-empty cells
    epochs_neuro   = epochs_neuro(nonEmptyIdx);
    epochs_audio   = epochs_audio(nonEmptyIdx);
    mapping_epochs = mapping_epochs(nonEmptyIdx,:);
    mapping_labels = mapping_labels(nonEmptyIdx);

    % Remove additional audio trials & trials that are shorter than
    % trial_dur
    %--------------------------------------------------------------
    n_files = length(epochs_neuro);

    epoch_counter = 0;
    min_dur       = inf; % minimal durations for trials
    for f_idx = 1:n_files
        % Audio is always longer than neuro due to experimental design of
        % trigger -> remove extra audio trials
        %------------------------------------------------------------------
        trials_neuro = epochs_neuro{f_idx}.trial;
        trials_audio = epochs_audio{f_idx}.trials;
        n_neuro      = length(trials_neuro);

        % Remove additional audio trials
        trials_audio = trials_audio(1:n_neuro);

        % Check duration for epochs that are shorter than 90% of trial_dur
        %------------------------------------------------------------------
        keep_idx     = true(1,n_neuro);
        dur_expected = fs_down*trialdur; % desired trial duration in samples
        for trl_idx = 1:n_neuro
            dur_neuro = size(trials_neuro{trl_idx},2);
            dur_audio = length(trials_audio{trl_idx}); % end column: audio

            % Flag as bad if trial duration < 90% of expected value
            % Check for last trial
            
            if dur_neuro < 0.99 * dur_expected   
                keep_idx(trl_idx) = false; 
            else % Check remaining trials similar durations
                % expcept 1 samples difference per second duration
                % (arbitrarly chosen)
                if abs(dur_neuro-dur_audio) > trialdur * 10
                    error("%s: Trial durations don't match",subject)
                end

                % Update minimum duration
                min_dur = min([min_dur, dur_neuro, dur_audio]);
            end

        end % loop over trials

        % Remove unwanted epochs
        %-----------------------
        epochs_audio{f_idx}.trials = trials_audio(keep_idx);

        cfg                 = [];
        cfg.trials          = keep_idx;
        epochs_neuro{f_idx} = ft_selectdata(cfg, epochs_neuro{f_idx});

        n_neuro = length(epochs_neuro{f_idx}.trial);
        n_audio = length(epochs_audio{f_idx}.trials);
        if ~isequal(n_neuro,n_audio)
            error('%s: Unexpected number of audio trials (%i) and neuro trials (%i)',subject,n_audio,n_neuro)
        end

        % Update mapping
        mapping_epochs(f_idx,:) = [epoch_counter+1,epoch_counter + n_neuro];
        epoch_counter           = epoch_counter + n_neuro; 

    end

    % Append data over all recordings
    %--------------------------------
    % Append MEG recordings
    cfg                = [];
    cfg.keepsampleinfo = 'no'; 
    epochs_neuro       = ft_appenddata(cfg, epochs_neuro{:});
    
    % Append Audio recordings
    epochs_audio_all = cell(1,length(epochs_audio));
    for f_idx = 1:n_files
        epochs_audio_all{f_idx} = epochs_audio{f_idx}.trials;
    end
    epochs_audio = [epochs_audio_all{:}];
    clear epochs_audio_all

    % Match length of trials
    %-----------------------
    n_trials = length(epochs_audio);

    for trl_idx = 1:n_trials
        epochs_neuro.trial{trl_idx} = epochs_neuro.trial{trl_idx}(:,1:min_dur);
        epochs_neuro.time{trl_idx}  = epochs_neuro.time{trl_idx}(1:min_dur);
        epochs_audio{trl_idx}       = epochs_audio{trl_idx}(1:min_dur);
    end

    %% Save results
    %--------------
    dir2save = fullfile(settings.path2derivatives,subject,'speech');
    if ~exist(dir2save,'dir')
        mkdir(dir2save)
    end
    fname = sprintf('%s_preprocessed_crosscorr.mat',subject);
    
    results                = struct;
    results.settings       = settings;
    results.epochs_neuro   = epochs_neuro;
    results.epochs_audio   = epochs_audio;
    results.mapping_epochs = mapping_epochs;
    results.mapping_label  = mapping_labels;
    
    save(fullfile(dir2save,fname),'results','-v7.3'); 
    % save(fullfile(dir2save,fname),"-fromstruct",results); 
    fprintf("\n%s from %s saved.\n",fname,subject)

end % Loop over subjects

%% Functions
%--------------------------------------------------------------------------

function [audiodata_epoched] = epoch_audiobook_crosscorr(cfg, audiodata)
%--------------------------------------------------------------------------
% Till Habersetzer, 27.06.2025
% Communication Acoustics, CvO University Oldenburg
% till.habersetzer@uol.de 
%
% EPOCH_AUDIOBOOK_CROSSCORR Segments and resamples raw audio data.
%
%   [audiodata_epoched] = epoch_audiobook_crosscorr(cfg, audiodata)
%   segments a continuous audio signal into consecutive epochs (trials) of a
%   fixed duration. Each resulting epoch is then resampled to a new target
%   sampling frequency. The function returns a data structure similar to the
%   FieldTrip raw data format.
%
%   Input:
%       cfg (struct): Configuration structure with the following fields:
%       .fs_audio (numeric) - Original sampling rate of the audio data (in Hz).
%       .fs_down  (numeric) - Target sampling rate for downsampling (in Hz).
%       .trialdur (numeric) - Desired duration of each epoch (in seconds).
%       audiodata (vector): Raw audio data, expected as a 1xT or Tx1 time-series.
%
%   Output:
%       audiodata_epoched (struct): A structure containing the epoched data:
%       .trials (cell) - A 1xN cell array where each cell contains one epoched 
%                        and resampled audio trial.
%       .trl (matrix)  - An Nx2 matrix where each row contains the start and
%                        end sample indices of an epoch [start, end]
%                        relative to the original audiodata.
%       .cfg (struct)  - The original configuration structure used.
%--------------------------------------------------------------------------

% Ensure audiodata is a row vector (1xT) for consistent processing.
if iscolumn(audiodata)
    audiodata = audiodata';
end

% Define Epoch Boundaries 
%------------------------
fs_audio = cfg.fs_audio;
trialdur = cfg.trialdur;
L        = length(audiodata);

% Calculate the length of each epoch in samples of the original audio
epoch_len_samples = round(trialdur * fs_audio);

% Determine the start sample for each trial
trl_starts    = (1:epoch_len_samples:L)';
% Determine the end sample for each trial
trl_ends      = trl_starts + epoch_len_samples - 1;
% Ensure the very last sample of the last trial does not exceed the audio length
trl_ends(end) = min(trl_ends(end), L);

% Create the trial definition matrix [start_sample, end_sample]
trl        = [trl_starts, trl_ends];
num_trials = size(trl, 1);

% Create and Resample Epochs 
%---------------------------
trials = cell(1, num_trials);
for trl_idx = 1:num_trials
    % Extract the audio segment for the current trial
    original_trial = audiodata(trl(trl_idx, 1) : trl(trl_idx, 2));

    % Resample the trial to the target frequency and store it
    trials{trl_idx} = resample(original_trial, cfg.fs_down, fs_audio);
end

% Package Data into Output Structure 
%-----------------------------------
audiodata_epoched        = struct();
audiodata_epoched.trials = trials;
audiodata_epoched.trl    = trl;
audiodata_epoched.cfg    = cfg;

end
