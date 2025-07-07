function preprocessing_audio(settings)
%--------------------------------------------------------------------------
% Till Habersetzer, 05.07.2025
% Communication Acoustics, CvO University Oldenburg
% till.habersetzer@uol.de 
%
% Description:
%   Prepares audio stimuli for a neural cross-correlation analysis by
%   converting raw audiobooks into epoched onset envelopes.
%
%   Key Steps:
%   -   Computes the onset envelope from the raw audio (.wav) files.
%   -   Segments the continuous envelope into fixed-length epochs and
%       downsamples each one to a common target frequency.
%   -   Saves the final data to two separate .mat files: one for the pilot
%       subject (sub-01), and another for all other subjects who shared
%       the same audio material.
%
% Input:
%   settings (struct): Settings structure containing all necessary paths 
%                      and processing parameters.
%--------------------------------------------------------------------------

%% Script settings 
%--------------------------------------------------------------------------

% Filter settings
filtertype = settings.crosscorr.filtertype;
lpfreq     = settings.crosscorr.audio.lpfreq; % audio
lpfiltord  = settings.crosscorr.audio.lpfiltord; % audio

% Epoch length
trialdur = settings.crosscorr.trialdur;

% Downsampling frequency
fs_down = settings.crosscorr.fs_down;

% Audio frequency
fs_audio = settings.fs_audio;

bids_dir = settings.path2bids;
stim_dir = fullfile(bids_dir,'stimuli','audiobooks'); % for audio data

%% Preprocess audio - same for most of the subjects
%-----------------------------------------------------------------------

% Audio from pilot 1 differs from rest
for sub_idx = 1:2

    if sub_idx == 1
        subject          = 'sub-01';
        stim_dir_subject = fullfile(stim_dir,'pilot');
    else
        subject          = 'sub-others';
        stim_dir_subject = stim_dir;
    end

    f_idx            = 1; % file counter
    audiobook_labels = cell(1,4);
    epochs_audio     = cell(1,4);

    for task_idx = 1:2     
        task = sprintf('audiobook%i',task_idx);
        for run_idx = 1:2
            run = sprintf('0%i',run_idx);
    
            % Import raw audio and compute envelope
            fname_audio             = sprintf('task-%s_run-%s',task,run);
            audiobook_labels{f_idx} = fname_audio;

            [raw_audiodata,fs] = audioread(fullfile(stim_dir_subject,[fname_audio,'_stim.wav'])); 
    
            if ~isequal(fs,fs_audio)
                error('Unexpected sampling frequency in audiofile (%i)!',fs)
            end

            % Compute onset envelope
            cfg              = [];
            cfg.type         = 'onset_envelope';
            cfg.fs           = fs_audio;
            cfg.lpfreq       = lpfreq;
            cfg.lpfiltord    = lpfiltord;
            cfg.filtertype   = filtertype;
            cfg.plotfiltresp = 'no';
            envelope         = cal_envelope(cfg, raw_audiodata);

            % Cut into epochs and downsampling
            cfg             = [];
            cfg.fs_audio    = fs_audio;
            cfg.trialdur    = trialdur;
            cfg.fs_down     = fs_down;
            envelope_epochs = epoch_audiobook_crosscorr(cfg, envelope);
            fprintf('%s: %s preprocessed.\n',subject,fname_audio)

            epochs_audio{f_idx} = envelope_epochs;
            clear envelope_epochs raw_audiodata envelope

            % Increment counter
            f_idx = f_idx + 1;

        end % loop over runs
    
    end % loop over tasks

    % Save results
    %--------------
    dir2save = fullfile(settings.path2derivatives,'stimuli');
    if ~exist(dir2save,'dir')
        mkdir(dir2save)
    end
    fname = sprintf('%s_preprocessed_audio_crosscorr.mat',subject);

    audio                  = struct();
    audio.epochs_audio     = epochs_audio;
    audio.audiobook_labels = audiobook_labels;

    save(fullfile(dir2save,fname),'audio','-v7.3'); 
    fprintf("%s from %s saved.\n",fname,subject)

end % loop over subjects

%% Functions
%--------------------------------------------------------------------------

function [audiodata_epoched] = epoch_audiobook_crosscorr(cfg,audiodata)
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
%       .trial (cell)  - A 1xN cell array where each cell contains one epoched 
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
L = length(audiodata);

% Calculate the length of each epoch in samples of the original audio
epoch_len_samples = round(cfg.trialdur * cfg.fs_audio);

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
    trials{trl_idx} = resample(original_trial, cfg.fs_down, cfg.fs_audio);
end

% Package Data into Output Structure 
%-----------------------------------
audiodata_epoched       = struct();
audiodata_epoched.trial = trials;
audiodata_epoched.trl   = trl;
audiodata_epoched.cfg   = cfg;

end

end % function end