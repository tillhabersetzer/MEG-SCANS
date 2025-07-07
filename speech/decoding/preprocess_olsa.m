%--------------------------------------------------------------------------
% Till Habersetzer, 04.07.2025
% Communication Acoustics, CvO University Oldenburg
% till.habersetzer@uol.de 
%
% To run from the command line (linux server):
% matlab -nodisplay -nosplash -r "preprocess_olsa; exit;"
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

% Select subjects
subjects = 1:24; % all files available for sub-02-sub-23
% subjects = 4;
n_subj   = length(subjects);

% Filter settings
filtertype = settings.decoding.filtertype; % windowed sinc type I linear phase FIR filter
bpfreq     = settings.decoding.bpfreq; % meg + audio

% Epoch length
trialdur = settings.decoding.trialdur;

% Downsampling frequency
fs_down = settings.decoding.fs_down;

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
% addpath(settings.path2fieldtrip)
% ft_defaults

% Addpath for additional functions
addpath(fullfile(settings.path2project,'analysis','helper_functions'))

% run script for auditory toolbox
run(fullfile(settings.path2amtoolbox,'amtstart'));



% Preprocessing Settings
%-----------------------
fs_audio  = settings.olsa.fs_audio;
fs_neuro  = settings.olsa.fs_neuro;
bp_freq   = settings.olsa.bp_freq;
fs_down   = settings.olsa.fs_down;
prestim   = settings.olsa.prestim;
poststim  = settings.olsa.poststim;

audio_dir = fullfile(settings.path2project,'experiment','audio');

%% Loop over files
%--------------------------------------------------------------------------

neuro_data = cell(1,n_runs);
audio_data = cell(1,n_runs);
snr_data   = cell(1,n_runs);
triallabel = zeros(n_runs,2);
filenames  = cell(1,n_runs);

for run_idx = 1:n_runs % loop over files

    run = sprintf('run-0%i',runs(run_idx));

    % Define datapath
    %----------------
    if use_maxfilter
        rawdata_path = fullfile(settings.path2project,'derivatives',subject,'maxfilter');
        filename     = sprintf('%s_%s_task-%s_%s_proc-tsss_meg.fif',subject,session,'olsa',run); 
        datapath     = fullfile(rawdata_path,filename);
    else
        rawdata_path = fullfile(settings.path2project,'rawdata',subject,session,'meg');
        filename     = sprintf('%s_%s_task-%s_%s_meg.fif',subject,session,'olsa',run); 
        datapath     = fullfile(rawdata_path,filename);
        % badchannels  = settings.badchannels; 
    end

    filenames{run_idx} = filename;

    settings_path = fullfile(settings.path2project,'rawdata',subject,session,'meg',sprintf('%s_%s_task-%s_%s_results.json',subject,session,'olsa',run));
    
    % Load continuous data
    %---------------------
    cfg              = [];
    cfg.dataset      = datapath;
    cfg.channel      = sensors; 
    cfg.continuous   = 'yes';
    cfg.coilaccuracy = 0; % non-empty value ensures that sensors are expressed in SI units
    data             = ft_preprocessing(cfg);

    if ~isequal(data.fsample,fs_neuro)
        error('Unexpected sampling frequency of neurophysiological data!')
    end

    % Filter data
    %------------
    
    % lowpass
    cfg            = [];
    % cfg.detrend    = 'yes';
    % cfg.demean     = 'yes';
    cfg.channel    = sensors;               
    cfg.continuous = 'yes';
    cfg.lpfilter   = 'yes';
    cfg.lpfreq     = bp_freq(2); 
    data           = ft_preprocessing(cfg, data);

    % highpass
    cfg            = [];
    % cfg.detrend    = 'yes';
    % cfg.demean     = 'yes'; 
    cfg.channel    = sensors;         
    cfg.continuous = 'yes';
    cfg.hpfilter   = 'yes';
    cfg.hpfreq     = bp_freq(1); 
    cfg.hpfiltord  = 4; % default: 6
    data           = ft_preprocessing(cfg, data);
      
    % Check spectrum
    %---------------
    if check_spectrum
        
        cfg         = [];
        cfg.subject = subject;
        cfg.channel = 'meg';
        plot_spectrum(cfg,data)

        if isfield(data,'elec')
            cfg        = [];
            cfg.elec   = data.elec;
            layout_eeg = ft_prepare_layout(cfg);
        end

        cfg            = [];
        cfg.subject    = subject;
        cfg.channel    = 'eeg';
        cfg.layout_eeg = layout_eeg;
        plot_spectrum(cfg,data)

    end

    % Epoch data
    %-----------
    % Load experiment settings and extract olsa sentence durations
    str          = fileread(settings_path); % dedicated for reading files as text
    exp_settings = jsondecode(str); % using the jsondecode function to parse JSON from string 
    durationlist = exp_settings.settings.olsa.run_cfg.(sprintf('run_%i',runs(run_idx))).durationlist; % in samples
    durationlist = durationlist/exp_settings.settings.soundmexpro.sampling_frequency; 
    triggerlist  = exp_settings.settings.olsa.run_cfg.(sprintf('run_%i',runs(run_idx))).triggerlist.('triggerlist_int');

    % Load and check intelligibilities and SNR
    % and correct snr label if necessary
    snrs_new                  = exp_settings.settings.olsa.trial_params.snr_computation.snrs;
    snrs_new(snrs_new == 100) = 0; % for sub-00

    if runs(run_idx) == runs(1)
        intelligibilities = exp_settings.settings.olsa.trial_params.snr_computation.intelligibilities;
        snrs              = snrs_new;
    else 
        if ~isequal(intelligibilities,exp_settings.settings.olsa.trial_params.snr_computation.intelligibilities)
            error("Unexpected intelligibilities!")
        end
        if ~isequal(snrs,snrs_new)
            error("Unexpected snrs!")
        end
    end
                
    % Epoching
    %---------
    cfg                             = [];
    cfg.dataset                     = datapath;
    cfg.trialfun                    = 'my_trialfun_olsa'; 
    cfg.trialdef.eventtype          = 'STI101';
    cfg.trialdef.prestim            = prestim;                
    cfg.trialdef.poststim           = poststim;  
    cfg.trialdef.sentence_durations = durationlist;  
    cfg.trialdef.triggerlist        = triggerlist;
    cfg                             = ft_definetrial(cfg);
    trl                             = cfg.trl;

    % Add latency correction 
    %-----------------------
    if apply_latency_correction
        trl(:,1:2) = trl(:,1:2) + round(audio_latency*1000); % in samples
    end

    cfg     = [];
    cfg.trl = trl;
    data    = ft_redefinetrial(cfg, data);

    % Downsample data
    %----------------
    cfg            = [];
    cfg.resamplefs = fs_down;
    cfg.detrend    = 'no';
    data           = ft_resampledata(cfg,data);

    % Collect data and create label for trials
    %-----------------------------------------
    % add information about triallabels
    if runs(run_idx) == runs(1)
        triallabel(run_idx,:) = [1,length(data.trial)];
    else
        triallabel(run_idx,:) = [triallabel(run_idx-1,2)+1,triallabel(run_idx-1,2)+length(data.trial)];
    end
    
    neuro_data{run_idx} = data;
    clear data 

    % Proprocess audio
    %-----------------
    playlist          = exp_settings.settings.olsa.run_cfg.(sprintf('run_%i',runs(run_idx))).playlist; % in samples
    snr_data{run_idx} = exp_settings.settings.olsa.run_cfg.(sprintf('run_%i',runs(run_idx))).snrlist; % in samples

    % Olsa - Compute Envelops
    %------------------------
    n_sentences = length(playlist);
    sentences   = cell(n_sentences,1);

    for stn_idx = 1:n_sentences % loop over sentences

        [audio, fs] = audioread(fullfile(audio_dir,'olsa','sentences',playlist{stn_idx}));
        if ~isequal(fs,fs_audio)
            error('Unexpected sampling frequency in audiofile!')
        end

        if size(audio, 2) == 2
            % Check if left and right channels are equal
            if ~isequal(audio(:,1), audio(:,2))
                error('Left and right channel are not the same!');
            end
            audio = audio(:,1);
        end

        % Calculate envelope
        [outsig,~] = auditoryfilterbank(audio,fs_audio,'flow',50,'fhigh',5000);
        audio      = mean(abs(outsig).^(0.6),2);
        % Add prestim/poststim interval
        audio      = [zeros(1,round(prestim*fs_audio)),audio',zeros(1,round(poststim*fs_audio))];

        % Lowpass against aliasing
        [b,a] = butter(6,300/(fs_audio/2),'low'); % 6
        % h     = fvtool(b,a);
        % h.Fs  = fsamp;
        audio = filtfilt(b,a,audio);
        % downsampling
        audio = resample(audio,fs_neuro,fs_audio);
        % filtering with same settings as in
        [b,a] = butter(6,max(bp_freq)/(fs_neuro/2),'low'); % 6
        audio = filtfilt(b,a,audio);
        [b,a] = butter(4,min(bp_freq)/(fs_neuro/2),'high'); % 4
        audio = filtfilt(b,a,audio);
        % Downsampling
        audio = resample(audio,fs_down,fs_neuro);

        sentences{stn_idx} = audio;
        clear audio
        fprintf('%s/ %s / Sentence %i processed.\n',session,run,stn_idx)
        
    end % sentences

    audio_data{run_idx} = sentences;
    clear sentences

end % runs

% Append data over trials and correct for differences in durations
%-----------------------------------------------------------------

% Append data
%------------
cfg                = [];
cfg.keepsampleinfo = 'no';
neuro_data         = ft_appenddata(cfg,neuro_data{:});

audio_data = vertcat(audio_data{:});
snr_data   = vertcat(snr_data{:});

% Match durations
%----------------
n_neuro = length(neuro_data.trial);
n_audio = length(audio_data);

if ~isequal(n_neuro,n_audio)
    error('Unexpected number of audio trials (%i) and neuro trials (%s)!',n_audio,n_neuro)
end

check_length = zeros(n_audio,2);
for trl_idx = 1:n_neuro
    check_length(trl_idx,1) = size(neuro_data.trial{trl_idx},2);
    check_length(trl_idx,2) = length(audio_data{trl_idx});
end

% Check whether distance is smaller than 1 sample
if any(abs(check_length(:,1)-check_length(:,2))>1)
    error("Trial durations don't match")
end

% Get index of minimum length for each trial
[min_durs,~] = min(check_length,[],2);

% Match length of trials
for trl_idx = 1:n_neuro
    min_dur                   = min_durs(trl_idx);
    neuro_data.trial{trl_idx} = neuro_data.trial{trl_idx}(:,1:min_dur);
    neuro_data.time{trl_idx}  = neuro_data.time{trl_idx}(1:min_dur);
    audio_data{trl_idx}       = audio_data{trl_idx}(1:min_dur);
end

% Save data
%-----------
data_preprocessed                            = struct();
data_preprocessed.neuro_data                 = neuro_data;
data_preprocessed.audio_data                 = audio_data;
data_preprocessed.snr_data                   = snr_data;
data_preprocessed.cfg_olsa.triallabel        = triallabel;
data_preprocessed.cfg_olsa.filenames         = filenames;
data_preprocessed.cfg_olsa.intelligibilities = intelligibilities;
data_preprocessed.cfg_olsa.snrs              = snrs;

if save_data
    fname    = sprintf('%s_%s_olsa_preproc.mat',subject,session);
    dir2save = fullfile(settings.path2derivatives,subject,'speech');
    if ~exist(dir2save,'dir')
        mkdir(dir2save)
    end
    save(fullfile(dir2save,fname),'data_preprocessed','-v7.3')  
end

clear data_preprocessed neuro_data audio_data snr_data triallabel intelligibilities filenames

end % end function