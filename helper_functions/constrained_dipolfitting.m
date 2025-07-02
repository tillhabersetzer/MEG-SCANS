function [scalar_dipmom] = constrained_dipolfitting(cfg,data)
%--------------------------------------------------------------------------
% Till Habersetzer (10.02.2022)
% Communication Acoustics, CvO University Oldenburg
% till.habersetzer@uol.de 
%
% This function performs dipolfitting for two dipoles. It computes scalar 
% dipolmoments for two dipoles with given locations and orientations. 
% The current fieldtrip version doesn't support dipolfitting with fixed
% orientations and consequently returns freely oriented dipoles with xyz-
% directions for the dipolmoment. 
% The orientation is computed as a mean over a precalculated dipolmoment
% timecourse (3,T) or can be predefined in cfg.ori.
% The solution of the overconstrained linear equation can be found in:
% Brain signals: Physics and Mathematics of MEG and EEG, p.97, eq. (5.23)
%
% Output:
%   scalar_dipmom: scalar dipolmoment for 2 dipols (2,T)
%
% Input: 
%   cfg.dipmom     dipol moments (6,T), 2 dipoles with xyz-components. Used
%                  for orientation estimation.
%   cfg.dippos:    dipol positions (2,3), each row contains a dipole 
%                  position
%   cfg.unit:      unit of dipol positions
%   cfg.headmodel: precomputed headmodel
%   cfg.channel:   channeltype for leadfield computation 'meg' 'megmag'
%                  'megplanar' 'eeg' etc.
%   cfg.noisecov:  optional, if provided: sphere the date using the noise 
%                  covariance matrix supplied
%   cfg.ori:       'mean' or 'svd' as options for dipole orientation 
%                  computation.
%                  'mean': dipole orientations are computed as a mean over 
%                          the freely oriented dipole timecourses (6,T) 
%                          over a specified time window in cfg.timewin. 
%                  'svd': dipole orientations are computed as the maximum 
%                         variance orientations via SVD over a specified
%                         time window in cfg.timewin.     
%   cfg.diptime    time axis of dipoles
%   cfg.timewin:   timewindow if dipole moment used for orientation 
%                  estimation. If not specified, the full time course
%                  of the dipole moment is used.
%   data:          averaged timecourse of data, usually computed via
%                  ft_timelockanalysis
%
% T:               number of timesamples
%--------------------------------------------------------------------------

dipmom    = cfg.dipmom;
dippos    = cfg.dippos; % dipole locations
diptime   = cfg.diptime;
headmodel = cfg.headmodel;
chantype  = cfg.channel;
unit      = cfg.unit;
ori       = cfg.ori;

% Extract Sphering information
%--------------------------------------------------------------------------
if isfield(cfg,'noisecov')
    noisecov = cfg.noisecov;
else
    noisecov = [];
end

% Extract dipole orientations 
%--------------------------------------------------------------------------
if isfield(cfg,'timewin')
    timewin = cfg.timewin;
    idx     = dsearchn(diptime',timewin');
    dipmom  = cfg.dipmom(:,idx(1):idx(2)); % 6xT dipole moments used to extract orientation
else
    dipmom  = cfg.dipmom;
end

switch ori
    case 'mean'
        ori1   = mean(dipmom(1:3,:),2); % mean dipolmoment
        ori1   = ori1/norm(ori1);       % make unit vector
        ori2   = mean(dipmom(4:6,:),2); % mean dipolmoment
        ori2   = ori2/norm(ori2);       % make unit vector
    case 'svd'
        moments   = dipmom(1:3,:);
        [u, ~, ~] = svd(moments, 'econ'); 
        ori1      = u(:,1); % already normalized
        moments   = dipmom(4:6,:);
        [u, ~, ~] = svd(moments, 'econ'); 
        ori2      = u(:,1); % already normalized
    otherwise 
        error('Selected orientation method is not supported!')
end

% Extract channel type of interest
%--------------------------------------------------------------------------
if contains(chantype,'meg')
    cfg         = [];
    cfg.channel = chantype;
    data        = ft_selectdata(cfg,data);
    sens        = data.grad;
elseif contains (chantype,'eeg')
    cfg         = [];
    cfg.channel = 'eeg*';
    data        = ft_selectdata(cfg,data);
    sens        = data.elec;
end

avg = data.avg;

% Sphere data using the noise covariance matrix
%--------------------------------------------------------------------------
if ~isempty(noisecov) && strcmp(chantype,'meg')
  disp('Sphering data with noise covariance.')
  [u, s]        = svd(noisecov);
  tol           = max(size(noisecov)) * eps(norm(s, inf));
  s             = diag(s);
  r1            = sum(s > tol) + 1;
  s(1:(r1 - 1)) = 1 ./ sqrt(s(1:(r1 - 1)));
  s(r1:end)     = 0;
  sphere        = diag(s) * u';
  % apply the sphering to the data
  avg           = sphere * avg;
  % apply the sphering as a pre-multiplication to the sensor definition
  montage          = [];
  montage.labelold = data.label;
  montage.labelnew = data.label;
  montage.tra      = sphere;
  sens             = ft_apply_montage(sens, montage, 'balancename', 'sphering');
end

% Computation of leadfield for predefined dipole positions
%--------------------------------------------------------------------------
cfg                       = [];
cfg.grad                  = sens;  
cfg.channel               = chantype;                       
cfg.sourcemodel.pos       = dippos;
cfg.sourcemodel.unit      = unit;
cfg.headmodel             = headmodel;                
% cfg.singleshell.batchsize = 5000;     
cfg.sourcemodel.inside    = [1,1]; % assume both dipoles are inside even if they are not
lf                        = ft_prepare_leadfield(cfg);

% Apply orientation constraint to leadfield
%--------------------------------------------------------------------------
lf1 = lf.leadfield{1}*ori1;
lf2 = lf.leadfield{2}*ori2;
L   = [lf1,lf2];

% Solve overconstrained linear equation for scalar dipolmoments
%--------------------------------------------------------------------------
% avg = L*mom ,  avg (Nch,T), L (Nch,2), mom (2,T)
% mom = inv(L'*L)*L' * avg
% operator = inv(L'*L)*L';
operator      = (L'*L)\L';
scalar_dipmom = operator*avg; % scalar dipolmoments

end