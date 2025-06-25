%--------------------------------------------------------------------------
% Till Habersetzer, 23.06.2025
% Communication Acoustics, CvO University Oldenburg
% till.habersetzer@uol.de 
%  
% This script automates the creation and visual validation of individual
% headmodels and sourcemodels for MEG source analysis. It processes
% T1-weighted anatomical MRIs for a batch of subjects.
%
% WORKFLOW:
% For each subject, the script performs the following sequential steps:
%   1.  **Coregistration**: Aligns the T1w MRI with the MEG sensor space.
%   2.  **Reslicing**: Standardizes the MRI volumes.
%   3.  **Segmentation**: Segments the volume into brain, skull, and scalp.
%   4.  **Headmodel**: Creates a single-shell conductive headmodel.
%   5.  **Sourcemodel**: Generates a subject-specific volumetric grid.
%   6.  **Visualization**: Creates a 4-panel figure to visually confirm the
%       alignment of the sensors, headshape, headmodel, and sourcemodel.
%
% REPORTING:
% Using the MATLAB Report Generator, the script compiles the visualization
% from every subject into a single HTML report file. This allows for
% efficient quality control of the coregistration for the entire batch.
%
% OUTPUTS:
%   - For each subject, three .mat files are saved: segmented MRI,
%     headmodel, and sourcemodel.
%   - A single 'report-coregistration.html' file containing all alignment plots.
%
% To run from the command line (linux server):
% matlab -nodisplay -nosplash -r "compute_headmodel_sourcemodel; exit;"
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

subjects = 2:24; % all files available for sub-02-sub-23
% subjects = 4;
n_subj   = length(subjects);

% Maxfilter
use_maxfilter = settings.use_maxfilter;
% Use task that has been used as a reference head position
task_ref_head = settings.task_ref_head;
parts         = split(task_ref_head,"_");
task          = parts(1);
run           = parts(2);

% Add fieldtrip
addpath(settings.path2fieldtrip)
ft_defaults

% Load template grid for source model
%--------------------------------------------------------------------------
template_grid = importdata(settings.path2template_grid);
template_grid = ft_convert_units(template_grid,'mm'); 
template_grid.coordsys = 'mni';

% Use atlas to create a binary mask
%----------------------------------
atlas = ft_read_atlas(fullfile(settings.path2fieldtrip,'template','atlas','aal','ROI_MNI_V4.nii'));

% Select all temporal labels including Heschl
% idx    = contains(atlas.tissuelabel,{'Temporal','Heschl'});
% labels = atlas.tissuelabel(idx);

cfg       = [];
cfg.atlas = atlas;
cfg.roi   = atlas.tissuelabel;  % here you can also specify a single label, i.e. single ROI
% cfg.roi   = labels;  % here you can also specify a single label, i.e. single ROI
mask      = ft_volumelookup(cfg, template_grid);

inside_grid          = false(template_grid.dim);
inside_grid(mask==1) = true;
template_grid.inside = reshape(inside_grid, [], 1);

% sum(template_grid.inside) % number of dipole positions
%--------------------------------------------------------------------------

% figure
% ft_plot_mesh(template_grid.pos(template_grid.inside,:));

%% Create a report
%-----------------
import mlreportgen.report.*
import mlreportgen.dom.*

report_dir            = fullfile(settings.path2project,'analysis','chirps','sourcelevel');
report_filename       = fullfile(report_dir,'report-coregistration');
report_coregistration = Report(report_filename,'html-file');

% Delete existing report file if it exists
if exist(report_filename, 'file')
    delete(report_filename);
end

fig_counter = 0; % counter for figures
files2print = cell(1,n_subj);

%% Loop over subjects
%--------------------------------------------------------------------------
% MRI and MEG have already been coregistered during BIDS conversion.
% Therefore, this step can be automated.

for sub_idx = 1:n_subj

    subject = sprintf('sub-%02d',subjects(sub_idx));

    % Add new subject to report
    para_subj          = Paragraph(sprintf('Coregistration: %s', subject));
    para_subj.Bold     = true;
    para_subj.Style    = {HAlign('center')}; % Center the title
    para_subj.FontSize = '50pt';
    append(report_coregistration,para_subj);
    spacer            = Paragraph(' ');
    spacer.Style      = {Height('100pt')};  

    % Choose MEG dataset for coregistration, only for later visualization here
    if use_maxfilter
        megfile = fullfile(settings.path2derivatives,subject,'maxfilter',sprintf('%s_task-%s_%s_proc-tsss-mc_meg.fif',subject,task,run));
    else
        megfile = fullfile(settings.path2bids,subject,'meg',sprintf('%s_task-%s_%s.fif',subject,task,run));
    end

    % Read in Mri
    %------------
    mrifile  = fullfile(settings.path2bids,subject,'anat',[subject,'_T1w.nii.gz']);
    mri_orig = ft_read_mri(mrifile,'readbids','no'); 
   
    % 1. Coregistration to Neuromag Coordsystem based on Anatomical Landmarks
    %------------------------------------------------------------------------
    % Load fiducial voxel coordinates from .json file
    str            = fileread(strrep(mrifile,'.nii.gz','.json')); 
    json_sidecar   = jsondecode(str); 
    anat_landmarks = json_sidecar.AnatomicalLandmarkCoordinates;

    % Coregistration
    cfg              = [];
    cfg.method       = 'fiducial';
    cfg.coordsys     = 'neuromag';
    cfg.fiducial.nas = anat_landmarks.NAS;
    cfg.fiducial.lpa = anat_landmarks.LPA;
    cfg.fiducial.rpa = anat_landmarks.RPA;
    mri_coreg        = ft_volumerealign(cfg, mri_orig);

    % You could check it as well here
    % headshape                 = ft_read_headshape(megfile,'unit','mm');
    % cfg                       = [];
    % cfg.method                = 'headshape';
    % cfg.headshape.headshape   = headshape;
    % cfg.coordsys              = 'neuromag';
    % cfg.headshape.interactive = 'yes';
    % cfg.headshape.icp         = 'no'; 
    % mri_coreg                 = ft_volumerealign(cfg, mri_coreg);

    % 2. Reslice
    %-----------
    if strcmp(subject,'sub-09') || strcmp(subject,'sub-20')
        % Segmentation seems to fail for resliced mri with [256 256 256]
        cfg            = [];
        % cfg.resolution = 1; % in mm 
        cfg.dim        = [320 320 320];
        cfg.coordsys   = 'neuromag';
        mri_resliced   = ft_volumereslice(cfg,mri_coreg); 
    else
        cfg            = [];
        % cfg.resolution = 1; % in mm 
        cfg.dim        = [256 256 256];
        cfg.yrange     = [-130,125]; % could be adjusted, see above remark, sum up to 255 -> 256 voxel
        cfg.coordsys   = 'neuromag';
        mri_resliced   = ft_volumereslice(cfg,mri_coreg); 
    end

    % Check 
    %------
    % cfg = [];
    % ft_sourceplot(cfg, mri_coreg);
    % cfg = [];
    % ft_sourceplot(cfg, mri_resliced);

    % 3. Segmentation (most time consuming)
    %--------------------------------------
    cfg           = [];
    % cfg.output    = {'brain'};
    cfg.output    = {'brain','skull','scalp'}; % you only need brain for meg
    mri_segmented = ft_volumesegment(cfg, mri_resliced);
    
    % Copy the Anatomy into the segmented Mri
    mri_segmented.anatomy = mri_resliced.anatomy;
   
    % Check
    %------
    cfg              = [];
    cfg.funparameter = 'brain';
    ft_sourceplot(cfg, mri_segmented);

    % 4. Headmodel
    %-------------
    cfg        = [];
    cfg.method = 'singleshell';
    vol        = ft_prepare_headmodel(cfg, mri_segmented);

    %% Create volumetric sourcemodel
    %----------------------------------------------------------------------
    % Coregistration of subject specific grid to the atlas based template grid
    cfg           = [];
    cfg.warpmni   = 'yes';
    cfg.template  = template_grid;
    cfg.nonlinear = 'yes';
    cfg.mri       = mri_segmented;
    cfg.unit      = 'mm';
    sourcemodel   = ft_prepare_sourcemodel(cfg);

    %% Save data
    %----------------------------------------------------------------------
    save(fullfile(settings.path2derivatives,subject,'chirp',sprintf('%s_T1w-segmented.mat',subject)),'mri_segmented'); % in mm
    save(fullfile(settings.path2derivatives,subject,'chirp',sprintf('%s_headmodel-singleshell.mat',subject)),'vol'); % in mm
    save(fullfile(settings.path2derivatives,subject,'chirp',sprintf('%s_sourcemodel-volumetric.mat',subject)),'sourcemodel'); % in mm

    %% Check alignment and add plot to report
    %----------------------------------------------------------------------
    % Same Units for Plot
    grad_mm      = ft_convert_units(ft_read_sens(megfile,'senstype','meg'),'mm'); 
    vol_mm       = ft_convert_units(vol,'mm');
    headshape_mm = ft_read_headshape(megfile,'unit','mm');

    % Add outer shape of head
    cfg             = [];
    cfg.tissue      = {'scalp'};
    cfg.method      = 'isosurface';
    cfg.numvertices = 10000;
    scalp           = ft_prepare_mesh(cfg, mri_segmented);

    angles = [-90, 0; 90, 0; 180, 0; 180, 90]; % left, right, frontal, from above

    fig_coregistration = figure('Visible', 'off', 'Units', 'pixels', 'Position', [0, 0, 1920, 1200]); % Create an invisible figure
    drawnow; % Ensures MATLAB processes the window state change
    t = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    for i = 1:4
        ax = nexttile; % Move to the next tile in the layout
        hold(ax, 'on'); % Use hold on for the specific axes  
        ft_plot_mesh(scalp, 'facecolor', 'skin')
        ft_plot_headmodel(vol_mm,  'facecolor', 'cortex', 'edgecolor', 'none');
        ft_plot_mesh(sourcemodel.pos(sourcemodel.inside,:)); % plot only locations inside the volume
        ft_plot_sens(grad_mm, 'style', '*b');
        ft_plot_headshape(headshape_mm);
        alpha 0.5;  % make the surface transparent
        view(angles(i,:)); 
        camlight; % Add a light source based on the camera position
        lighting phong % This command chooses the type of light bulb
        material dull % This command sets the surface finish of the objects
        hold(ax, 'off');
    end
    sgtitle(subject)

    % Export graphics
    %----------------
    fig_counter              = fig_counter + 1;
    fname                    = sprintf('%s_coregistration.png',subject);
    files2print{fig_counter} = fullfile(report_dir,fname);
    exportgraphics(fig_coregistration, files2print{fig_counter}, 'Resolution', 300); % High-quality PNG
    fprintf('\nGraphic %s exported.',fname)

    % Add figure to report 
    %---------------------
    img       = Image(files2print{fig_counter});
    % img.Style = {ScaleToFit(true), Width('100%')};
    append(report_coregistration, img);
    append(report_coregistration, spacer);

end

%% Close report
%--------------------------------------------------------------------------
close(report_coregistration);
% rptview(report_coregistration);

% Delete the image file after adding it to the report
for idx = 1:fig_counter
    delete(files2print{idx});
end