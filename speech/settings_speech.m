% settings

% ensure that we don't mix up settings
clear settings

% Define paths
%-------------
settings.rootpath = fullfile('/mnt','localSSDPOOL','fiko7761'); % server
% settings.rootpath = fullfile('M:'); % local computer

settings.path2project     = fullfile(settings.rootpath,'masterthesis');
settings.path2derivatives = fullfile(settings.path2project,'derivatives');
settings.path2bids        = fullfile(settings.path2project,'bids_conversion','bids_data');
settings.path2fieldtrip   = fullfile(settings.rootpath,'toolboxes','fieldtrip-20250523');
settings.path2amtoolbox   = fullfile(settings.rootpath,'toolboxes','amtoolbox-1.6.0');
settings.path2mtrftoolbox = fullfile(settings.rootpath,'toolboxes','mTRF-Toolbox-2.7','mtrf');

% Other settings
%---------------
settings.use_maxfilter            = true;
settings.apply_latency_correction = true;
settings.audio_latency            = 3/1000; % 3 ms (see measurements 28.06.24)

% Speech
%-------
% Cross-correlation
settings.fs_audio                  = 44100;
settings.crosscorr.trialdur        = 10; % 10s
settings.crosscorr.filtertype      = 'but';
settings.crosscorr.meg.bpfreq      = [0.5,45];
settings.crosscorr.meg.bpfiltord   = 4; % default FieldTrip (6 reported in Pertersen paper leads to unstable filter)
settings.crosscorr.audio.lpfreq    = 25;
settings.crosscorr.audio.lpfiltord = 3;
settings.crosscorr.fs_down         = 250;
settings.crosscorr.zscore          = false;
settings.crosscorr.timelags        = [-100,900]; % in ms
% Decoder
settings.decoding.trialdur         = 120;
settings.decoding.filtertype       = 'firws';
settings.decoding.bpfreq           = [0.5,4];
settings.decoding.fs_neuro         = 1000;
settings.decoding.fs_down          = 64;
settings.decoding.zscore           = true;
settings.decoding.olsa.prestim     = 0;
settings.decoding.olsa.poststim    = 0.5;
settings.decoding.olsa.padding     = 10; %s zeropad olsa wav files before filtering (in sec)

% mTRF-model hyperparameters
settings.decoding.decoder.drct        = -1; % direction: backward model
settings.decoding.decoder.tmin        = 0; % minimum time lag (ms)
settings.decoding.decoder.tmax        = 75; % maximum time lag (ms)
settings.decoding.decoder.lambda_exp  = -14:2:10; % exponent regularization values
settings.decoding.decoder.zeropad     = 1; % zero-pad the outer rows of the design matrix or delete them
settings.decoding.decoder.corr_metric = 'Spearman'; % 'Pearson' specifiy accuracy metric
settings.decoding.decoder.fast        = 1; % use the fast cross-validation method (requires more memory)
settings.decoding.decoder.dim         = 2; % work along the rows, observations in columns
settings.decoding.decoder.type        = 'multi'; % use all lags simultaneously to fit a multi-lag model


