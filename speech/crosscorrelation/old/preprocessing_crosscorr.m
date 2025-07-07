%--------------------------------------------------------------------------
% Till Habersetzer, 05.07.2025
% Communication Acoustics, CvO University Oldenburg
% till.habersetzer@uol.de 
%
% Description:
%   Prepares multi-run MEG and corresponding audio envelope data for a
%   subject-level cross-correlation analysis.
%
%   Key Processing Steps:
%   -   Loads continuous MEG data (.fif) and pre-processed audio envelopes.
%   -   Filters the MEG data and epochs it based on experimental triggers
%       using the FieldTrip toolbox.
%   -   Downsamples the epoched MEG data to a target frequency.
%   -   Synchronizes the MEG and audio trials by rejecting pairs that are
%       too short or mismatched.
%   -   Appends the cleaned data from all runs and truncates every trial
%       to a uniform length for sample-by-sample consistency.
%   -   Saves the final 'epochs_neuro' and 'epochs_audio' data structures
%       into a single .mat file per subject.
%
% To run from the command line (linux server):
% matlab -nodisplay -nosplash -r "preprocessing_crosscorr; exit;"
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
% subjects = 1;
n_subj   = length(subjects);

% Filter settings
filtertype = settings.crosscorr.filtertype;
bpfreq     = settings.crosscorr.meg.bpfreq; % meg
bpfiltord  = settings.crosscorr.meg.bpfiltord; % meg

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

% Loop over subjects, both audiobooks, including both runs
for sub_idx = 1:n_subj

    subject = sprintf('sub-%02d',subjects(sub_idx));
    
    % Initialize data
    epochs_neuro   = cell(1,4);
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
                %---------------
                % based on fieldtrip events - searching in STI101 rather
                % than importing bids trigger
                cfg                     = [];
                cfg.trialfun            = 'my_trialfun_audiobook'; 
                cfg.dataset             = event_path; % path should always point to original bids-path
                cfg.trialdef.method     = 'fieldtrip';
                cfg.trialdef.trialdur   = trialdur; % 10 (s)
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
                cfg              = [];
                cfg.channel      = 'meg';
                cfg.dataset      = data_path;
                cfg.bpfilter     = 'yes';
                cfg.bpfilttype   = filtertype;
                cfg.bpfiltord    = 4; %bpfiltord; 
                cfg.bpfreq       = bpfreq;
                cfg.plotfiltresp = 'no';
                data             = ft_preprocessing(cfg); 

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

                % Count epochs
                n_epochs                = length(epoch.trial);
                mapping_epochs(f_idx,:) = [epoch_counter+1,epoch_counter + n_epochs];
                epoch_counter           = epoch_counter + n_epochs;
                mapping_labels{f_idx}   = sprintf('task-%s_run-%s',task,run);

                %% Downsampling
                %--------------
                cfg                 = [];
                cfg.resamplefs      = fs_down;
                cfg.detrend         = 'no';
                epochs_neuro{f_idx} = ft_resampledata(cfg,epoch);
                % clear epoch

            end % file exists

            % Increment counter
            f_idx = f_idx + 1;

        end % loop over runs
    end % loop over tasks

    %% Import audio data
    %-------------------------------------------------------------------
    if strcmp(subject,'sub-01')
        fname = 'sub-01_preprocessed_audio_crosscorr.mat';
    else
        fname = 'sub-others_preprocessed_audio_crosscorr.mat';
    end
    audio            = importdata(fullfile(settings.path2derivatives,'stimuli',fname));
    epochs_audio     = audio.epochs_audio;
    audiobook_labels = audio.audiobook_labels;
    clear audio

    %% Append data over trials and correct for differences in number of trials and trial durations
    %---------------------------------------------------------------------------------------------
    % Keep only non-empty cells
    nonEmptyIdx = ~cellfun('isempty',mapping_labels);

    % Compare mapping between audio and neura data and continue only if identical
    check_mapping = isequal(mapping_labels(nonEmptyIdx),audiobook_labels(nonEmptyIdx));
    if ~check_mapping
        error('%s: Mismatch in mapping betwenn neuro and audio files!');
    end

    % Delete empty cells if data is missing - otherwise appenddata won't
    % work, Keep only non-empty cells
    epochs_neuro   = epochs_neuro(nonEmptyIdx);
    epochs_audio   = epochs_audio(nonEmptyIdx);
    mapping_epochs = mapping_epochs(nonEmptyIdx,:);
    mapping_labels = mapping_labels(nonEmptyIdx);

    % Remove additional audio trials & trials that are shorter than
    % trial_dur
    %----------------------------------------------------------------------
    n_files = length(epochs_neuro);

    epoch_counter = 0;
    min_dur       = inf; % minimal durations for trials
    for f_idx = 1:n_files
        % Audio is always longer than neuro due to experimental design of
        % trigger -> remove extra audio trials
        %------------------------------------------------------------------
        trials_neuro = epochs_neuro{f_idx}.trial;
        trials_audio = epochs_audio{f_idx}.trial;
        n_neuro      = length(trials_neuro);

        % Remove additional audio trials
        trials_audio = trials_audio(1:n_neuro);

        % Check duration for epochs that are shorter than 99% of trial_dur
        %------------------------------------------------------------------
        is_trial_valid = true(1,n_neuro);
        dur_expected   = fs_down*trialdur; % desired trial duration in samples
        for trl_idx = 1:n_neuro
            dur_neuro = size(trials_neuro{trl_idx},2);
            dur_audio = length(trials_audio{trl_idx}); % end column: audio

            % Condition 1: Trial is too short
            if dur_neuro < 0.99 * dur_expected   
                is_trial_valid(trl_idx) = false;
                continue; % Skip to the next trial
            end

             % Condition 2: Sanity check for audio/neuro async (arbitrary threshold)
            if abs(dur_neuro-dur_audio) > trialdur * 2
                error('%s: Mismatch in trial duration between audio and neuro is too large for trial %d.', subject, trl_idx);
            end

            % Update minimum duration
            min_dur = min([min_dur, dur_neuro, dur_audio]);

        end % loop over trials

        % Remove unwanted epochs
        %-----------------------
        epochs_audio{f_idx}.trial = trials_audio(is_trial_valid);

        cfg                 = [];
        cfg.trials          = is_trial_valid;
        epochs_neuro{f_idx} = ft_selectdata(cfg, epochs_neuro{f_idx});

        n_neuro = length(epochs_neuro{f_idx}.trial);
        n_audio = length(epochs_audio{f_idx}.trial);
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
        epochs_audio_all{f_idx} = epochs_audio{f_idx}.trial;
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
    settings.min_dur       = min_dur;
    results.epochs_neuro   = epochs_neuro;
    results.epochs_audio   = epochs_audio;
    results.mapping_epochs = mapping_epochs;
    results.mapping_label  = mapping_labels;
    
    save(fullfile(dir2save,fname),'results','-v7.3');  
    fprintf("\n%s from %s saved.\n",fname,subject)

end % Loop over subjects