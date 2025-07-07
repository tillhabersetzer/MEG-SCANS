%--------------------------------------------------------------------------
% Till Habersetzer, 23.06.2025
% Communication Acoustics, CvO University Oldenburg
% till.habersetzer@uol.de 
%
% This script performs a quality control (QC) check by visualizing the
% outputs of the MRI processing and MEG/MRI coregistration pipeline. It is
% based on the FieldTrip toolbox coregistration quality control example.
% https://www.fieldtriptoolbox.org/example/coregistration_quality_control/
%
% For a single subject, the script loads the pre-computed headmodel,
% segmented MRI, sourcemodel, and the corresponding raw MEG data (for sensor
% and headshape information). It then generates a series of figures to
% visually verify the alignment between these different components, such as:
%
%   - Sensor locations, headshape, headmodel, and sourcemodel.
%   - Brain segmentation and headmodel boundaries on the anatomical MRI.
%   - Digitized headshape points on the MRI-derived scalp surface.
%
% The script does not process any data; it is intended solely for visual
% inspection and verification.
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

sub_idx = 13;
subject = sprintf('sub-%02d', sub_idx);

% Maxfilter
use_maxfilter = settings.use_maxfilter;
% Use task that has been used as a reference head position
task_ref_head = settings.task_ref_head;
parts         = split(task_ref_head,"_");
task          = parts(1);
run           = parts(2);

%% Import data
%--------------------------------------------------------------------------
% Choose MEG dataset for coregistration, only for later visualization here
if use_maxfilter
    megfile = fullfile(settings.path2derivatives,subject,'maxfilter',sprintf('%s_task-%s_%s_proc-tsss-mc_meg.fif',subject,task,run));
else
    megfile = fullfile(settings.path2bids,subject,'meg',sprintf('%s_task-%s_%s.fif',subject,task,run));
end

% All in mm
headmodel   = importdata(fullfile(settings.path2derivatives,subject,'chirp',sprintf('%s_headmodel-singleshell.mat',subject)));
mri         = importdata(fullfile(settings.path2derivatives,subject,'chirp',sprintf('%s_T1w-segmented.mat',subject)));
sourcemodel = importdata(fullfile(settings.path2derivatives,subject,'chirp',sprintf('%s_sourcemodel-volumetric.mat',subject)));
headshape   = ft_read_headshape(megfile,'unit','mm'); % has to be the file for reference head position
grad        = ft_convert_units(ft_read_sens(megfile,'senstype','meg'),'mm'); 

%% figure 1, All together

figure
ft_plot_sens(grad)
ft_plot_headshape(headshape)
ft_plot_headmodel(headmodel)
alpha 0.5;  % camlight;
alpha 0.4;  % make the surface transparent
ft_plot_mesh(sourcemodel.pos(sourcemodel.inside,:)); 
ft_plot_axes([], 'unit', 'mm');

%% figure 2, MRI anatomy and brain segmentation

cfg              = [];
cfg.anaparameter = 'anatomy';
cfg.funparameter = 'brain';
% cfg.funparameter = 'scalp';
cfg.location     = [0 0 60];
ft_sourceplot(cfg, mri)

%% figure 3 and 4, MRI anatomy and headmodel

location = [0 0 60];
figure
ft_plot_ortho(mri.anatomy, 'transform', mri.transform, 'location', location, 'intersectmesh', headmodel.bnd,'unit','mm')

figure
ft_plot_montage(mri.anatomy, 'transform', mri.transform, 'intersectmesh', headmodel.bnd,'unit','mm')

%% figure 5, MRI scalp surface and polhemus headshape
% problems with surface extraction due to "worse" mri segmentation likely

cfg             = [];
cfg.tissue      = {'scalp'};
cfg.method      = 'isosurface';
cfg.numvertices = 10000;
scalp           = ft_prepare_mesh(cfg, mri);

figure
ft_plot_mesh(scalp, 'facecolor', 'skin')
lighting phong
camlight left
camlight right
material dull
alpha 0.5
ft_plot_headshape(headshape, 'vertexcolor', 'k');

%% figure 6, MRI and anatomical landmarks
% seems not to work :/
figure
for i=1:3
  subplot(2,2,i)
  figure
  title(headshape.fid.label{i});
  location = headshape.fid.pos(i,:);
  ft_plot_ortho(mri.anatomy, 'transform', mri.transform, 'style', 'intersect', 'location', location, 'plotmarker', location, 'markersize', 5, 'markercolor', 'y')
end

%% figure 7, MRI scalp surface and anatomical landmarks

figure
ft_plot_mesh(scalp, 'facecolor', 'skin')
lighting phong
camlight left
camlight right
material dull
alpha 0.3
ft_plot_mesh(headshape.fid, 'vertexcolor', 'k', 'vertexsize', 10);