%--------------------------------------------------------------------------
% Till Habersetzer, 19.06.2025
% Communication Acoustics, CvO University Oldenburg
% till.habersetzer@uol.de 
%
% This script performs a complete Auditory Evoked Field (AEF) analysis on 
% MEG data using the FieldTrip toolbox. For a list of subjects, it 
% automatically loads raw data from a BIDS-like structure, defines trials 
% around chirp stimuli, and applies a full preprocessing pipeline: artifact 
% rejection (jump, clip, z-value), filtering, and baseline correction. 
% The script processes magnetometers and gradiometers separately. It concludes 
% by saving the time-locked average for each subject and a final grand 
% average across all participants.
%
% To run from the command line (linux server):
% matlab -nodisplay -nosplash -r "compute_aefs; exit;"
%--------------------------------------------------------------------------

close all
clearvars
clc 

%% Import main settings 
%--------------------------------------------------------------------------
current_dir = pwd;
cd(fullfile('..'))
settings_chirps
cd(current_dir)

%% Script settings 
%--------------------------------------------------------------------------
% For test purposes of eeg reference, only select single session and
% condition

% Select seubjects
subjects = 2:24; % all files available for sub-02-sub-23
% subjects = 4;
n_subj   = length(subjects);

% Event Value
eventvalue = settings.trig_id.chirp;

% Filter settings
bpfreq  = settings.bpfreq;

% Epoch length
prestim  = settings.trialdef.prestim;
poststim = settings.trialdef.poststim;

% Threshold
z_value_threshold = settings.z_value_threshold;

% Use maxfiltered files
% check whether maxfiltered data should be analyzed
use_maxfilter = settings.use_maxfilter;

% Latency Correction
apply_latency_correction = settings.apply_latency_correction;
audio_latency            = settings.audio_latency;

bids_dir = settings.path2bids;

% Add fieldtrip
addpath(settings.path2fieldtrip)
ft_defaults

%% Processing
%--------------------------------------------------------------------------
% Save data for grand average
avgs_mag  = cell(1,n_subj);
avgs_grad = cell(1,n_subj);

for sub_idx = 1:n_subj

    subject = sprintf('sub-%02d',subjects(sub_idx));

    % Collect data
    epochs        = cell(1,4);
    n_trials      = zeros(2,4); % 4 audiobooks (before, after rejection), 3rd row: 
    n_trials_mag  = zeros(1,2); % combined, after rejection
    n_trials_grad = zeros(1,2); % combined, after rejection
 
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
                % if ~isfile(datapath)
                %     filename = sprintf('%s_%s_task-transient_proc-tsss_meg.fif',subject,session); 
                %     datapath = fullfile(rawdata_path,filename);
                % end
        
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
                cfg                     = [];
                cfg.dataset             = event_path;
                cfg.trialfun            = 'ft_trialfun_general'; 
                cfg.trialdef.eventtype  = {'up-chirp'};
                cfg.trialdef.eventvalue = eventvalue; 
                cfg.trialdef.prestim    = prestim;         
                cfg.trialdef.poststim   = poststim;         
                cfg                     = ft_definetrial(cfg);
            
                trl = cfg.trl;
                % Before rejection
                n_trials(1,f_idx) = size(trl,1); % count trials

                %% Detect bad trials
                %-------------------
                cfg                        = [];
                cfg.trl                    = trl;
                cfg.dataset                = data_path;
                cfg.artfctdef.jump.channel = 'meg'; 
                [~, artifact_jump]         = ft_artifact_jump(cfg);
                
                cfg                        = [];
                cfg.trl                    = trl;
                cfg.dataset                = data_path;
                cfg.artfctdef.clip.channel = 'meg'; 
                [~, artifact_clip]         = ft_artifact_clip(cfg);
                
                cfg                         = [];
                cfg.trl                     = trl;
                cfg.dataset                 = data_path;
                cfg.artfctdef.jump.artifact = artifact_jump;
                cfg.artfctdef.clip.artifact = artifact_clip;
                cfg                         = ft_rejectartifact(cfg);
                trl                         = cfg.trl;
                % After rejection
                n_trials(2,f_idx) = size(trl,1); % count trials

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

                % SSP only applicable to non-maxfiltered data
                if ~use_maxfilter
                    % cfg            = [];
                    % cfg.ssp        = 'all';
                    % cfg.trials     = 'all';
                    % cfg.updatesens = 'yes';
                    % data_ssp       = ft_denoise_ssp(cfg, data);
                    % rename to data instead of data_ssp to be useful
                end

                %% Epoch data
                %------------
                cfg           = [];
                cfg.trl       = trl;  
                epochs{f_idx} = ft_redefinetrial(cfg,data);

                % Increment counter
                f_idx = f_idx + 1;

            end % ffile exists
        end % loop over runs
    end % loop over tasks

    %% Append data
    %----------------------------------------------------------------------
    cfg                = [];
    cfg.keepsampleinfo = 'no'; 
    epochs_all         = ft_appenddata(cfg, epochs{:});

    % clear epochs

    %% Split data in magetometers and gradiometers
    %---------------------------------------------
    cfg              = [];
    cfg.channel      = 'megplanar';
    epochs_grad      = ft_selectdata(cfg, epochs_all);
    n_trials_grad(1) = length(epochs_grad.trial); % count trials

    cfg.channel      = 'megmag';
    epochs_mag       = ft_selectdata(cfg, epochs_all);
    n_trials_mag(1)  = length(epochs_mag.trial); % count trials

    % clear epochs_all

    %% (Semi-)automatic artifact rejection
    %-------------------------------------
    % separately for gradiometers
    % cfg             = [];
    % cfg.metric      = 'maxzvalue';
    % cfg.channel     = 'megplanar';
    % cfg.keepchannel = 'yes';  % This keeps those channels that are not displayed in the data
    % epochs_grad      = ft_rejectvisual(cfg,epochs_grad);
    
    % separately for magnetometers
    % cfg.channel     = 'megmag';
    % epochs_mag      = ft_rejectvisual(cfg,epochs_mag);

    cfg                              = [];
    cfg.artfctdef.zvalue.channel     = 'megplanar';
    cfg.artfctdef.zvalue.interactive = 'no'; % make plot interactive
    cfg.artfctdef.zvalue.cutoff      = z_value_threshold;
    [~,artifacts_grad]               = ft_artifact_zvalue(cfg,epochs_grad);

    cfg.artfctdef.zvalue.channel     = 'megmag';
    [~,artifacts_mag]                = ft_artifact_zvalue(cfg,epochs_mag);

    cfg                              = [];
    cfg.artfctdef.zvalue.artifact    = artifacts_grad;
    epochs_grad                      = ft_rejectartifact(cfg,epochs_grad);
    n_trials_grad(2)                 = length(epochs_grad.trial); % count trials

    cfg.artfctdef.zvalue.artifact    = artifacts_mag;
    epochs_mag                       = ft_rejectartifact(cfg,epochs_mag);
    n_trials_mag(2)                  = length(epochs_mag.trial); % count trials

    %% Apply baseline correction
    %---------------------------
    cfg                = [];
    cfg.baselinewindow = [-prestim 0];
    cfg.demean         = 'yes';
    epochs_grad        = ft_preprocessing(cfg,epochs_grad);
    epochs_mag         = ft_preprocessing(cfg,epochs_mag);
    
    % %% downsample data
    % %-----------------
    % if downsample_status
    %     cfg               = [];
    %     cfg.resamplefs    = fs_down;
    %     cfg.detrend       = 'no';
    %     epochs_grad = ft_resampledata(cfg,epochs_grad);
    %     epochs_mag  = ft_resampledata(cfg,epochs_mag);
    % end

    %% Add timelocked average
    %-----------------------
    cfg                = [];
    avgs_grad{sub_idx} = ft_timelockanalysis(cfg, epochs_grad);
    avgs_mag{sub_idx}  = ft_timelockanalysis(cfg, epochs_mag);

    %% Save results
    %--------------
    dir2save = fullfile(settings.path2derivatives,subject,'chirp');
    if ~exist(dir2save,'dir')
        mkdir(dir2save)
    end
    fname    = sprintf('%s_avgs.mat',subject);
    
    results               = struct;
    results.settings      = settings;
    results.avg_grad      = avgs_grad{sub_idx};
    results.avg_mag       = avgs_mag{sub_idx};
    results.n_trials      = n_trials;
    results.n_trials_grad = n_trials_grad;
    results.n_trials_mag  = n_trials_mag;
    
    save(fullfile(dir2save,fname),'results','-v7.3'); 
    fprintf("\n%s from %s saved.\n",fname,subject)

end % Loop over subjects

%% Add grand average over all subjects
%--------------------------------------------------------------------------
subject = 'grandaverage';

cfg         = [];
cfg.latency = 'all';
gavg_grad   = ft_timelockgrandaverage(cfg,avgs_grad{:});
gavg_mag    = ft_timelockgrandaverage(cfg,avgs_mag{:});

dir2save = fullfile(settings.path2derivatives,subject,'chirp');
if ~exist(dir2save,'dir')
    mkdir(dir2save)
end
fname    = sprintf('%s_avgs.mat',subject);

results           = struct;
results.settings  = settings;
results.gavg_grad = gavg_grad;
results.gavg_mag  = gavg_mag;

save(fullfile(dir2save,fname),'results','-v7.3'); 
fprintf("\n%s from %s saved.\n",fname,subject)