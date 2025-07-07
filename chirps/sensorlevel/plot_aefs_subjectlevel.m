%--------------------------------------------------------------------------
% Till Habersetzer, 20.06.2025
% Communication Acoustics, CvO University Oldenburg
% till.habersetzer@uol.de 
%
% This script serves as a visualization tool for pre-computed Auditory 
% Evoked Field (AEF) data, using the FieldTrip toolbox. It loads the 
% time-locked average data for a single, specified subject and generates a 
% comprehensive set of figures to inspect the results. Visualizations 
% include multi-channel sensor layouts (ft_multiplotER), single-channel 
% waveforms (ft_singleplotER), and 2D topographic maps (ft_topoplotER). 
% Each plot type is generated for magnetometers, planar gradiometers, and 
% their combined planar gradiometer representation.
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
% Choose subject for plotting
subject = 'sub-21'; % 

%% Load data
%--------------------------------------------------------------------------
data     = importdata(fullfile(settings.path2project,'derivatives',subject,'chirp',[subject,'_avgs.mat']));   
avg_grad = data.avg_grad;
avg_mag  = data.avg_mag;

%% Add combined planar gradiometers
%--------------------------------------------------------------------------
cfg        = [];
cfg.method = 'sum';
avg_cmb    = ft_combineplanar(cfg,avg_grad);

%% Plot of event related fields for all sensors arranged topographically 
%--------------------------------------------------------------------------

% Magnetometer
%-------------
figure
cfg            = [];
cfg.showlabels = 'yes';
cfg.fontsize   = 6;
cfg.layout     = 'neuromag306mag_helmet.lay';
ft_multiplotER(cfg, avg_mag);
title(sprintf('%s: magnetometer',subject))

% Gradiometer
%------------
figure
cfg            = [];
cfg.showlabels = 'yes';
cfg.fontsize   = 6;
cfg.layout     = 'neuromag306planar_helmet.lay';
ft_multiplotER(cfg, avg_grad);
title(sprintf('%s: gradiometer',subject))

% Combined Planar Gradiometers
%-----------------------------
figure
cfg            = [];
cfg.showlabels = 'yes';
cfg.fontsize   = 6;
cfg.layout     = 'neuromag306cmb_helmet.lay';
ft_multiplotER(cfg, avg_cmb);
title(sprintf('%s: combined gradiometer',subject))

%% Plot one specific channel
%--------------------------------------------------------------------------
% Same condition as above
chan2plot = 'MEG1221'; % e.g. magnetometer

cfg         = [];
cfg.channel = chan2plot;
ft_singleplotER(cfg, avg_mag);
title(sprintf('%s: %s',subject,chan2plot))

%% Plot the topographic distribution of the data
%--------------------------------------------------------------------------
timewin = [0.1, 0.13]; % time window to look at in more detail

figure
cfg          = [];
cfg.xlim     = timewin;
cfg.style    = 'straight';
cfg.comment  = 'no';
cfg.marker   = 'off';
cfg.colorbar = 'southoutside';
cfg.style    = 'both'; % both colormap and contour lines

cfg.layout   = 'neuromag306mag_helmet.mat';
cfg.figure   = subplot(1,3,1);
ft_topoplotER(cfg, avg_mag);
title('Magnetometer');

cfg.layout   = 'neuromag306planar_helmet.mat';
cfg.figure   = subplot(1,3,2);
ft_topoplotER(cfg, avg_grad);
title('Gradiometer');

cfg.layout   = 'neuromag306cmb_helmet.mat';
cfg.figure   = subplot(1,3,3);
ft_topoplotER(cfg, avg_cmb);
title('Combined Gradiometer');

sgtitle(sprintf('%s: [%d-%d] ms',subject, timewin(1)*1000,timewin(2)*1000))
