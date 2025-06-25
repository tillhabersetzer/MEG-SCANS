%--------------------------------------------------------------------------
% Till Habersetzer, 28.04.2025
% Communication Acoustics, CvO University Oldenburg
% till.habersetzer@uol.de
% 
% This script performs dipolfitting with a loop over all subjects
%--------------------------------------------------------------------------

close all
clear
clc

%% Import main settings 
%--------------------------------------------------------------------------
current_dir = pwd;
cd(fullfile('..','..'))
settings_analysis
cd(current_dir)

%% Script settings
%--------------------------------------------------------------------------
subjects = [1,2,3];
n_subj   = length(subjects);

% Switch between erf-types
%-------------------------
timewin = [0.02 0.05]; % 20-50ms covering N19m-P30m-Peak2Peak interval

chantype = 'megplanar'; % 'megmag','megplanar','meg'

conditions_audio     = {'Click','Up'};
conditions_earphones = {'Tip300','Sensimetrics'};

% Maxfilter
use_maxfilter = settings.use_maxfilter;

% Filter settings
bpfreq  = settings.bpfreq;
dftfreq = settings.dftfreq;

%% Preprocess data
%--------------------------------------------------------------------------

for sub_idx = 1:n_subj % loop over subjects

    subject = sprintf('sub-%02d', subjects(sub_idx));

    % Load sourcemodel and headmodel
    %-------------------------------
    % Proceed in SI-units. Leads apparantly to correct units for dipole moments in Am

    sourcemodel = importdata(fullfile(settings.path2project,'derivatives',subject,'forward_modelling',[subject,'_sourcemodel-volumetric.mat']));
    sourcemodel = ft_convert_units(sourcemodel, 'm'); 

    headmodel = importdata(fullfile(settings.path2project,'derivatives',subject,'forward_modelling',[subject,'_headmodel-singleshell.mat']));
    headmodel = ft_convert_units(headmodel, 'm');  

    % Load data
    %----------
    data       = importdata(fullfile(settings.path2project,'derivatives',subject,'sensorlevel',[subject,'_erfs.mat']));
    avgs       = data.avgs;
    conditions = data.conditions;
    clear data

    % Loop over two earphones
    %------------------------
    source_vec          = cell(2,2); % 2 ear-phones x 2 conditions (up,click)
    source_sca_mean     = cell(2,2); 
    source_sca_svd      = cell(2,2);  
    source_pooled_sym   = cell(1,2); % 2 ear-phones 
    source_pooled_nosym = cell(1,2);

    % Apply Sphering if magnetometers and gradiometers are used
    % together
    %----------------------------------------------------------
    sphering = strcmp(chantype,'meg');
    if sphering     
        if use_maxfilter
            noisefile = fullfile(settings.path2derivatives,subject,'maxfilter',sprintf('%s_task-%s_proc-sss_meg.fif',subject,'emptyroom'));
        else
            noisefile = fullfile(settings.path2project,'rawdata',subject,'meg',sprintf('%s_task-%s.fif',subject,'emptyroom'));
        end

        cfg              = [];
        cfg.dataset      = noisefile;
        cfg.channel      = 'meg'; 
        cfg.continuous   = 'yes';
        cfg.bpfilter     = 'yes';
        cfg.bpfreq       = bpfreq;
        cfg.dftfilter    = 'yes';      % enable notch filtering to eliminate power line noise
        cfg.dftfreq      = dftfreq;    % set up the frequencies for notch filtering
        cfg.dftreplace   = 'zero';     % 'zero' is default'
        cfg.coilaccuracy = 0;          % ensure that sensors are expressed in SI units
        noise            = ft_preprocessing(cfg); 

        cfg                  = [];
        cfg.removemean       = 'yes'; % default for covariance computation
        cfg.covariance       = 'yes';
        cfg.covariancewindow = 'all';
        avg_noise            = ft_timelockanalysis(cfg,noise);

        % imagesc(avg_noise.cov)
    end

    for ep_idx = 1:2 % Loop over earphones

        earphone = conditions_earphones{ep_idx};
    
        %% Fit initial dipole model to pooled conditions
        %------------------------------------------------------------------     

        % Get pooled condition
        con_idx = (contains(conditions,earphone) & contains(conditions,'pooled'));
    
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
        cfg.channel     = chantype;
    
        if sphering; cfg.dipfit.noisecov = avg_noise.cov; end 
        source_pooled_sym{ep_idx} = ft_dipolefitting(cfg, avgs{con_idx});

        % Now that we have a better starting point for the dipole fit, we can release the symmetry constraint
        %----------------------------------------------------------------------------------------------------
        cfg                  = [];
        cfg.latency          = timewin;
        cfg.numdipoles       = 2;
        cfg.symmetry         = []; % no symmetry constraint
        cfg.gridsearch       = 'no';
        cfg.nonlinear        = 'yes'; % non-linear search for optimal dipole parameters
        cfg.dip.pos          = source_pooled_sym{ep_idx}.dip.pos;
        cfg.sourcemodel.unit = sourcemodel.unit; % defines units of dipole positions
        cfg.headmodel        = headmodel;
        cfg.senstype         = 'meg';
        cfg.channel          = chantype;
    
        if sphering; cfg.dipfit.noisecov = avg_noise.cov; end 
        source_pooled_nosym{ep_idx} = ft_dipolefitting(cfg, avgs{con_idx});

        %% Reconstruct time courses of all conditions with different approaches
        %----------------------------------------------------------------------
        
        % 1. Approach: Fit freely oriented dipole (fixed positions, loose orientation)
        %-----------------------------------------------------------------------------
        cfg                  = [];
        cfg.latency          = 'all';
        cfg.numdipoles       = 2;
        cfg.symmetry         = [];
        cfg.nonlinear        = 'no';  % use a fixed position
        cfg.gridsearch       = 'no';
        cfg.dip.pos          = source_pooled_nosym{ep_idx}.dip.pos; % use estimated dipole positions from pooled data as reference
        cfg.sourcemodel.unit = sourcemodel.unit; % defines units of dipole positions
        cfg.headmodel        = headmodel;
        cfg.channel          = chantype; 
        cfg.senstype         = 'meg';
        
        if sphering; cfg.dipfit.noisecov = avg_noise.cov; end 
        con_idx              = (contains(conditions,earphone) & contains(conditions,conditions_audio{1}));
        source_vec{ep_idx,1} = ft_dipolefitting(cfg, avgs{con_idx}); % estimation of amplitude and orientation
        con_idx              = (contains(conditions,earphone) & contains(conditions,conditions_audio{2}));
        source_vec{ep_idx,2} = ft_dipolefitting(cfg, avgs{con_idx});

        % 2. Approach: Fit fixed oriented dipole (fixed positions, fixed orientation)
        %----------------------------------------------------------------------------
        % dipole orientations are estimated based on previously computed
        % dipolmoments (xyz-components)
        % mean orientation is used to return single dipolar timecourse
        cfg           = [];
        cfg.headmodel = headmodel;
        cfg.channel   = chantype;     
        cfg.dippos    = source_pooled_nosym{ep_idx}.dip.pos;
        cfg.dipmom    = source_pooled_nosym{ep_idx}.dip.mom;
        cfg.diptime   = source_pooled_nosym{ep_idx}.time;
        cfg.unit      = sourcemodel.unit;
        cfg.ori       = 'mean';
        cfg.timewin   = timewin;
        
        if sphering; cfg.noisecov = avg_noise.cov; end 
        con_idx                   = (contains(conditions,earphone) & contains(conditions,conditions_audio{1}));
        source_sca_mean{ep_idx,1} = constrained_dipolfitting(cfg, avgs{con_idx}); % estimation of amplitude and orientation
        con_idx                   = (contains(conditions,earphone) & contains(conditions,conditions_audio{2}));
        source_sca_mean{ep_idx,2} = constrained_dipolfitting(cfg, avgs{con_idx});

        % 3. Approach: Fit fixed oriented dipole (fixed positions, fixed orientation)
        %----------------------------------------------------------------------------
        % dipole orientations are estimated based on previously computed
        % dipolmoments (xyz-components) 
        % maximum variance orientation via SVD is used to return single dipolar timecourse
        cfg           = [];
        cfg.headmodel = headmodel;
        cfg.channel   = chantype;     
        cfg.dippos    = source_pooled_nosym{ep_idx}.dip.pos;
        cfg.dipmom    = source_pooled_nosym{ep_idx}.dip.mom;
        cfg.diptime   = source_pooled_nosym{ep_idx}.time;
        cfg.unit      = sourcemodel.unit;
        cfg.ori       = 'svd';
        cfg.timewin   = timewin;
        
        if sphering; cfg.noisecov = avg_noise.cov; end 
        con_idx                  = (contains(conditions,earphone) & contains(conditions,conditions_audio{1}));
        source_sca_svd{ep_idx,1} = constrained_dipolfitting(cfg, avgs{con_idx}); % estimation of amplitude and orientation
        con_idx                  = (contains(conditions,earphone) & contains(conditions,conditions_audio{2}));
        source_sca_svd{ep_idx,2} = constrained_dipolfitting(cfg, avgs{con_idx});

    end % Loop over earphones

    %% Create structure for saving data
    %----------------------------------------------------------------------

    % Save reconstructed waveforms and dipole positions
    %--------------------------------------------------
    dir2save = fullfile(settings.path2project,'derivatives',subject,'sourcelevel');
    if ~exist(dir2save, 'dir')
        mkdir(dir2save)
    end

    data                      = [];
    data.chantype             = chantype;
    data.conditions_audio     = conditions_audio;
    data.conditions_earphones = conditions_earphones;
    % dipole moments
    data.source_pooled_sym   = source_pooled_sym;
    data.source_pooled_nosym = source_pooled_nosym;
    % vector dipol moments of conditions
    data.source_vec          = source_vec;
    % scalar dipol moments of conditions
    data.source_sca_mean     = source_sca_mean;
    % scalar dipol moments of conditions
    data.source_sca_svd      = source_sca_svd;
    
    save(fullfile(dir2save,[subject,'_dipolefits.mat']),'-struct','data');

    clear source_pooled_sym source_pooled_nosym source_vec source_sca_mean source_sca_svd data

end % loop over subjects
