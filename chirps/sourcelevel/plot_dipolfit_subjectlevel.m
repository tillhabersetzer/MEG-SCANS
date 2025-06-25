%--------------------------------------------------------------------------
% Till Habersetzer, 24.06.2025
% Communication Acoustics, CvO University Oldenburg
% till.habersetzer@uol.de 
%
% This script is for the visualization and quality control of computed
% dipole fitting results for a single subject. It loads the output from the
% main dipole fitting analysis and generates a series of diagnostic plots to
% inspect the model's validity and interpret the findings. The
% visualization can be done for magnetometers and gradiometers.
%
% KEY VISUALIZATIONS:
%
% 1.  **Dipole Location (Native Space):** Plots the dipole locations from both
%     the initial symmetric fit and the refined non-linear fit on the
%     subject's anatomical MRI.
%
% 2.  **Anatomical Labeling (MNI Space):** Transforms the dipole coordinates
%     to MNI space and uses an atlas (AAL) to identify the corresponding
%     anatomical brain regions.
%
% 3.  **Source Time Courses:** Plots the reconstructed dipole activity for each
%     hemisphere, allowing for comparison between:
%       - A freely-oriented vector model (x, y, and z moments).
%       - A fixed-orientation scalar model.
%
% USAGE:
% The user needs to specify a single subject ID in the "Script settings"
% section. The output consists of several MATLAB figures for inspection.
%
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
subject   = 'sub-04'; % 'sub-02'
chan2plot = 'megplanar'; % 'megmag' 'megplanar'
time2plot = [-100 500]; % for dipole timecourse

% Addpath for additional functions
addpath(fullfile(settings.path2project,'analysis','helper_functions'))

%% Load data
%--------------------------------------------------------------------------
data = importdata(fullfile(settings.path2project,'derivatives',subject,'chirp',sprintf('%s_dipolfits.mat',subject))); 

% Select correct data based on channel selection
switch chan2plot
    case 'megmag'
        source_sym      = data.source_sym_mag;
        source_nosym    = data.source_nosym_mag;
        source_vec      = data.source_vec_mag;
        source_sca_mean = data.source_sca_mean_mag;
    case 'megplanar'
        source_sym      = data.source_sym_grad;
        source_nosym    = data.source_nosym_grad;
        source_vec      = data.source_vec_mag;
        source_sca_mean = data.source_sca_mean_mag;
end
clear data

sourcemodel   = importdata(fullfile(settings.path2derivatives,subject,'chirp',sprintf('%s_sourcemodel-volumetric.mat',subject)));
sourcemodel   = ft_convert_units(sourcemodel, 'm'); 
template_grid = importdata(settings.path2template_grid);
template_grid = ft_convert_units(template_grid,'mm'); 

%% Rescale units - optional for plotting (Am -> nAm)
%--------------------------------------------------------------------------
source_vec.dip.mom = 10^9*source_vec.dip.mom;
source_sca_mean    = 10^9*source_sca_mean;

%% (1) Inspect dipole locations
%--------------------------------------------------------------------------

% load mri
%---------
mri = importdata(fullfile(settings.path2derivatives,subject,'chirp',sprintf('%s_T1w-segmented.mat',subject)));
% Keep it in SI-units
mri = ft_convert_units(mri,'m'); 

figure('Name',subject)
hold on

% symmetric dipole fit 
ft_plot_dipole(source_sym.dip.pos(1,:), mean(source_sym.dip.mom(1:3,:),2), 'color', 'r', 'unit','m')
ft_plot_dipole(source_sym.dip.pos(2,:), mean(source_sym.dip.mom(4:6,:),2), 'color', 'r', 'unit','m')

% refinement of symmetric dipole fit 
ft_plot_dipole(source_nosym.dip.pos(1,:), mean(source_nosym.dip.mom(1:3,:),2), 'color', 'g', 'unit','m')
ft_plot_dipole(source_nosym.dip.pos(2,:), mean(source_nosym.dip.mom(4:6,:),2), 'color', 'g', 'unit','m')

pos = mean(source_sym.dip.pos,1);
ft_plot_slice(mri.anatomy, 'transform', mri.transform, 'location', pos, 'orientation', [1 0 0], 'resolution', 0.001)
ft_plot_slice(mri.anatomy, 'transform', mri.transform, 'location', pos, 'orientation', [0 1 0], 'resolution', 0.001)
ft_plot_slice(mri.anatomy, 'transform', mri.transform, 'location', pos, 'orientation', [0 0 1], 'resolution', 0.001)

ft_plot_crosshair(pos, 'color', [1 1 1]/2);

axis tight
axis off

view(12, -10)
title(subject)

%% (2) Inspect dipole locations in mni space
%--------------------------------------------------------------------------
dippos_loc_sym   = cell(1,2); % earphones x left/right
dippos_loc_nosym = cell(1,2);

% Positons
%---------
% positions in mni space
dippos_sym     = source_sym.dip.pos;
idx            = dsearchn(sourcemodel.pos,dippos_sym);
dippos_sym_mni = template_grid.pos(idx,:);

% positions in mni space
dippos_nosym     = source_nosym.dip.pos;
idx              = dsearchn(sourcemodel.pos,dippos_nosym);
dippos_nosym_mni = template_grid.pos(idx,:);

% this path needs to point to your fieldtrip path
atlas = ft_read_atlas(fullfile(settings.path2fieldtrip,'template','atlas','aal','ROI_MNI_V4.nii')); % mm

cfg            = [];
cfg.atlas      = atlas;
cfg.inputcoord = 'mni';
cfg.output     = 'label';

positions = vertcat(dippos_sym_mni,dippos_nosym_mni);
for p=1:2 % left / right
    positions           = dippos_sym_mni;
    cfg.roi             = positions(p,:); 
    labels              = ft_volumelookup(cfg, atlas);
    [~, indx]           = max(labels.count);
    dippos_loc_sym{p}   = labels.name(indx);

    positions           = dippos_nosym_mni;
    cfg.roi             = positions(p,:); 
    labels              = ft_volumelookup(cfg, atlas);
    [~, indx]           = max(labels.count);
    dippos_loc_nosym{p} = labels.name(indx);
end  

%% (3) Inspect timecourses 
%--------------------------------------------------------------------------

% (3.1) All moments, xyz-directions and both hemispheres
%-------------------------------------------------------

% Mapping between dipole locations and hemispheres
%-------------------------------------------------
pos     = dippos_nosym;
% [1,2] if first dipole belongs to left hemisphere, [2,1] if first dipole belongs to right hemisphere
% 1: left, 2: right
mapping = check_diploc(pos); 

idxs      = [];
idxs(1,:) = 1:3; % left
idxs(2,:) = 4:6; % right

name    = {'left hemisphere','right hemisphere'};
timevec = source_vec.time*1000; 
minval  = min(source_vec.dip.mom,[],'all');
maxval  = max(source_vec.dip.mom,[],'all');

figure('Name','Freely oriented dipole (fixed positions, loose orientation')
for i = 1:2
    idx     = idxs(mapping(i),:);
    axisvec = horzcat(time2plot,[minval,maxval]);
    
    subplot(2,1,i); 
    plot(timevec, source_vec.dip.mom(idx,:), '-')
    if i ==2
        xlabel('t / ms')
    end
    ylabel('dipole moment / nAm'); 
    legend({'x', 'y', 'z'});
    axis(axisvec) 
    grid on
    title(name{i})
    
end
sgtitle(sprintf('%s: %s',subject, chan2plot))

% (3.2) Visualize fixed dipole (fixed orientation) 
%-------------------------------------------------
% mean dipolmoment orientation has been used as orientation constraint

% Mapping between dipole locations and hemispheres
%-------------------------------------------------
pos     = dippos_nosym;
mapping = check_diploc(pos); 

name    = {'left hemisphere','right hemisphere'};
timevec = source_vec.time*1000; 
minval  = min(source_sca_mean,[],'all');
maxval  = max(source_sca_mean,[],'all');

figure('Name','Fixed oriented dipole (fixed positions, fixed orientation')
for i = 1:2

    idx     = mapping(i);
    axisvec = horzcat(time2plot,[minval,maxval]);
    
    subplot(2,1,i); 
    plot(timevec, source_sca_mean(idx,:), '-')
    if i ==2
        xlabel('t / ms')
    end
    ylabel('dipole moment / nAm')
    axis(axisvec)
    grid on
    title(name{i})
end
sgtitle(sprintf('%s: %s',subject, chan2plot))