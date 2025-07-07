% settings

% ensure that we don't mix up settings
clear settings

% Define paths
%-------------
% settings.path2project     = fullfile('/mnt','localSSDPOOL','fiko7761','masterthesis');
settings.path2project     = fullfile('M:\masterthesis');
settings.path2derivatives = fullfile(settings.path2project,'derivatives');
settings.path2bids        = fullfile(settings.path2project,'bids_conversion','bids_data');

% settings.path2fieldtrip   = fullfile('/mnt','localSSDPOOL','fiko7761','toolboxes','fieldtrip-20250523');
settings.path2fieldtrip   = fullfile('M:','toolboxes','fieldtrip-20250523');

% Other settings
%---------------
settings.use_maxfilter            = true;
settings.apply_latency_correction = true;
settings.audio_latency            = 3/1000; % 3 ms 
settings.apply_epoch_rejection    = false;

% Chirp
%------
settings.trig_id.chirp     = 2;
settings.trialdef.prestim  = 0.1;  
settings.trialdef.poststim = 0.5;  
settings.bpfreq            = [1,40];
settings.z_value_threshold = 6;

% Forward model / dipolfit
%-------------------------
% settings.path2template_grid = importdata(fullfile(settings.path2fieldtrip,'template','sourcemodel','standard_sourcemodel3d10mm.mat'));
settings.path2template_grid = fullfile(settings.path2fieldtrip,'template','sourcemodel','standard_sourcemodel3d5mm.mat');
settings.task_ref_head      = "audiobook1_run-01";
settings.timewindow_dipfit  = [0.090 0.130];

