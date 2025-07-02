%--------------------------------------------------------------------------
% Till Habersetzer, 24.06.2025
% Communication Acoustics, CvO University Oldenburg
% till.habersetzer@uol.de 
%
% This script performs dipole model fitting on averaged MEG data for a
% batch of subjects to estimate the location and time course of underlying
% neural sources, presumably in the bilateral auditory cortices.
%
% The script implements a multi-step fitting strategy for both magnetometer
% and planar gradiometer data separately:
%
% 1.  **Symmetric Grid Search:** An initial localization is performed by fitting
%     two dipoles simultaneously with a symmetry constraint across the x-axis.
%     This provides a robust starting point for the auditory sources.
%
% 2.  **Non-linear Refinement:** Using the positions from the grid search, a
%     non-linear fit is performed without the symmetry constraint to find the
%     optimal dipole locations more precisely.
%
% 3.  **Time Course Reconstruction:** With the final dipole positions fixed, the
%     script reconstructs the source time courses using two approaches:
%     a) **Vector Model:** Fits a dipole with a freely rotating orientation at
%        each time point using `ft_dipolefitting`.
%     b) **Scalar Model:** Projects the data onto a single, fixed orientation
%        (derived from the mean orientation over the fit time window) using a
%        custom `constrained_dipolfitting` function.
%
% OUTPUTS:
% For each subject, a single .mat file is saved, containing a structure with
% the results from all fitting stages for both sensor types.
%  
% Checkout for additional information: https://www.fieldtriptoolbox.org/tutorial/source/dipolefitting/
%
% To run from the command line (linux server):
% matlab -nodisplay -nosplash -r "compute_dipolfit; exit;"
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

subjects = 22:24; % all files available for sub-02-sub-23
% subjects = 4;
n_subj   = length(subjects);

% Time window for dipolfit
timewin = settings.timewindow_dipfit;

% Add fieldtrip
addpath(settings.path2fieldtrip)
ft_defaults

% Addpath for additional functions
addpath(fullfile(settings.path2project,'analysis','helper_functions'))

%% Dipolfit - Loop over subjects
%--------------------------------------------------------------------------

for sub_idx = 1:n_subj

    subject = sprintf('sub-%02d',subjects(sub_idx));

    % Import data
    %------------
    % Proceed in SI-units. Leads apparantly to correct units for dipole moments in Am
    sourcemodel = importdata(fullfile(settings.path2derivatives,subject,'chirp',sprintf('%s_sourcemodel-volumetric.mat',subject)));
    sourcemodel = ft_convert_units(sourcemodel, 'm'); 

    headmodel = importdata(fullfile(settings.path2derivatives,subject,'chirp',sprintf('%s_headmodel-singleshell.mat',subject)));
    headmodel = ft_convert_units(headmodel, 'm');  

    data      = importdata(fullfile(settings.path2derivatives,subject,'chirp',sprintf('%s_avgs.mat',subject)));
    avg_mag   = data.avg_mag;
    avg_grad  = data.avg_grad;
    clear data

    % Symmetrical gridsearch
    %-----------------------
    cfg             = [];
    cfg.latency     = timewin;
    cfg.numdipoles  = 2; % fit 2 dipoles
    cfg.symmetry    = 'x'; % expect activity in both auditory cortices, fit symmetric dipole
    cfg.sourcemodel = sourcemodel;
    cfg.headmodel   = headmodel;
    cfg.gridsearch  = 'yes';
    cfg.senstype    = 'meg';
    % Gradiometer
    cfg.channel     = 'megplanar';
    source_sym_grad = ft_dipolefitting(cfg, avg_grad);
    % Magnetometer
    cfg.channel     = 'megmag';
    source_sym_mag  = ft_dipolefitting(cfg, avg_mag);

    % Now that we have a better starting point for the dipole fit, we can release the symmetry constraint
    %----------------------------------------------------------------------------------------------------
    cfg                  = [];
    cfg.latency          = timewin;
    cfg.numdipoles       = 2;
    cfg.symmetry         = []; % no symmetry constraint
    cfg.gridsearch       = 'no';
    cfg.nonlinear        = 'yes'; % non-linear search for optimal dipole parameters
    cfg.sourcemodel.unit = sourcemodel.unit; % defines units of dipole positions
    cfg.headmodel        = headmodel;
    cfg.senstype         = 'meg';
    % Gradiometer
    cfg.dip.pos          = source_sym_grad.dip.pos;
    cfg.channel          = 'megplanar';
    source_nosym_grad    = ft_dipolefitting(cfg, avg_grad);
    % Magnetometer
    cfg.dip.pos          = source_sym_mag.dip.pos;
    cfg.channel          = 'megmag';
    source_nosym_mag     = ft_dipolefitting(cfg, avg_mag);

    % Reconstruct time courses of all conditions with different approaches
    %---------------------------------------------------------------------
    % 1. Approach: Fit freely oriented dipole (fixed positions, loose orientation)
    cfg                  = [];
    cfg.latency          = 'all';
    cfg.numdipoles       = 2;
    cfg.symmetry         = [];
    cfg.nonlinear        = 'no'; % use a fixed position
    cfg.gridsearch       = 'no';
    cfg.sourcemodel.unit = sourcemodel.unit; 
    cfg.headmodel        = headmodel;
    cfg.senstype         = 'meg';
    % Gradiometer
    cfg.dip.pos          = source_nosym_grad.dip.pos; 
    cfg.channel          = 'megplanar';
    source_vec_grad      = ft_dipolefitting(cfg, avg_grad); % estimation of amplitude and orientation
    % Magnetometer
    cfg.dip.pos          = source_nosym_mag.dip.pos; 
    cfg.channel          = 'megmag';
    source_vec_mag       = ft_dipolefitting(cfg,avg_mag); % estimation of amplitude and orientation

    % 2. Approach: Fit fixed oriented dipole (fixed positions, fixed orientation)
    %----------------------------------------------------------------------------
    % dipole orientations are estimated based on previously computed
    % dipolmoments (xyz-components)
    % mean orientation is used to return single dipolar timecourse
    cfg                  = [];
    cfg.headmodel        = headmodel;  
    cfg.unit             = sourcemodel.unit;
    cfg.ori              = 'mean';
    cfg.timewin          = timewin;
    % Gradiometer
    cfg.channel          = 'megplanar'; 
    cfg.dippos           = source_nosym_grad.dip.pos;
    cfg.dipmom           = source_nosym_grad.dip.mom;
    cfg.diptime          = source_nosym_grad.time;
    source_sca_mean_grad = constrained_dipolfitting(cfg, avg_grad); % estimation of amplitude and orientation
    % Magnetometer
    cfg.channel          = 'megmag'; 
    cfg.dippos           = source_nosym_mag.dip.pos;
    cfg.dipmom           = source_nosym_mag.dip.mom;
    cfg.diptime          = source_nosym_mag.time;
    source_sca_mean_mag  = constrained_dipolfitting(cfg, avg_mag); % estimation of amplitude and orientation

    %% Save results
    %--------------
    dir2save = fullfile(settings.path2derivatives,subject,'chirp');
    if ~exist(dir2save,'dir')
        mkdir(dir2save)
    end
    fname    = sprintf('%s_dipolfits.mat',subject);
    
    results                      = struct;
    results.source_sym_grad      = source_sym_grad;
    results.source_sym_mag       = source_sym_mag;
    results.source_nosym_grad    = source_nosym_grad;
    results.source_nosym_mag     = source_nosym_mag;
    results.source_vec_grad      = source_vec_grad;
    results.source_vec_mag       = source_vec_mag;
    results.source_sca_mean_grad = source_sca_mean_grad;
    results.source_sca_mean_mag  = source_sca_mean_mag;

    save(fullfile(dir2save,fname),'results','-v7.3'); 
    fprintf("\n%s from %s saved.\n",fname,subject)

end % Loop over subjects