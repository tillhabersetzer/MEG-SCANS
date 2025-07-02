% settings

% ensure that we don't mix up settings
clear settings

% Define paths
%-------------
settings.path2project     = fullfile('/mnt','localSSDPOOL','fiko7761','masterthesis');
% settings.path2project     = fullfile('M:\masterthesis');
% settings.path2project     = fullfile('/Volumes','ssdpool','masterthesis');
settings.path2derivatives = fullfile(settings.path2project,'derivatives');
settings.path2bids        = fullfile(settings.path2project,'bids_conversion','bids_data');

settings.path2fieldtrip   = fullfile('/mnt','localSSDPOOL','fiko7761','toolboxes','fieldtrip-20250523');
% settings.path2fieldtrip   = fullfile('M:','toolboxes','fieldtrip-20250523');
% settings.path2fieldtrip   = fullfile('/Volumes','ssdpool','toolboxes','fieldtrip-20250523');

% Other settings
%---------------
settings.use_maxfilter            = true;
settings.apply_latency_correction = true;
settings.audio_latency            = 3/1000; % 3 ms (see measurements 28.06.24)

% Speech
%-------
settings.fs_audio                     = 44100;
settings.crosscorr.trialdur           = 10; % 10s
settings.crosscorr.meg.bpfreq         = [0.5,45];
settings.crosscorr.audio.lpfreq.audio = 25;
settings.crosscorr.fs_down            = 250;
settings.crosscorr.zscore             = false;
settings.crosscorr.timelags           = [-100,900]; % in ms



