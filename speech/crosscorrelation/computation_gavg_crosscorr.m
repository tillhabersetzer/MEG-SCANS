function computation_gavg_crosscorr(subjects,settings)
%--------------------------------------------------------------------------
% Till Habersetzer, 05.07.2025
% Communication Acoustics, CvO University Oldenburg
% till.habersetzer@uol.de 
% 
% Description:
%   Computes a group-level grand average from individual-subject temporal
%   cross-correlation results.
%
%   Key Steps:
%   -   Loads the averaged cross-correlation data for each specified
%       subject from their pre-computed result files.
%   -   Calculates the grand average across all subjects for both the actual
%       and the shuffled control conditions using 'ft_timelockgrandaverage'.
%   -   Saves the final grand-averaged results to a single .mat file.
%
% Inputs:
%   subjects (vector): A list of subject numbers to include in the average.
%   settings (struct): Settings structure with all necessary paths and 
%                      parameters.
%--------------------------------------------------------------------------

%% Script settings 
%--------------------------------------------------------------------------
n_subj = length(subjects);

%% Import data for each subject
%--------------------------------------------------------------------------
avgs_crosscorr          = cell(1,n_subj);
avgs_crosscorr_shuffled = cell(1,n_subj);

for sub_idx = 1:n_subj
    subject = sprintf('sub-%02d',subjects(sub_idx));

    data                             = importdata(fullfile(settings.path2derivatives,subject,'speech',sprintf('%s_crosscorr.mat',subject)));   
    avgs_crosscorr{sub_idx}          = data.avg_crosscorr;
    avgs_crosscorr_shuffled{sub_idx} = data.avg_crosscorr_shuffled;
    clear data
    fprintf('%s loaded.\n',subject)

end

%% Compute grand-average cross-correlation 
%--------------------------------------------------------------------------
subject = 'grandaverage';

cfg                     = [];
cfg.latency             = 'all';
gavg_crosscorr          = ft_timelockgrandaverage(cfg,avgs_crosscorr{:});
gavg_crosscorr_shuffled = ft_timelockgrandaverage(cfg,avgs_crosscorr_shuffled{:});

dir2save = fullfile(settings.path2derivatives,subject,'speech');
if ~exist(dir2save,'dir')
    mkdir(dir2save)
end
fname = sprintf('%s_crosscorr.mat',subject);

results                         = struct;
results.n_subjects              = n_subj;
results.gavg_crosscorr          = gavg_crosscorr;
results.gavg_crosscorr_shuffled = gavg_crosscorr_shuffled;

save(fullfile(dir2save,fname),'results','-v7.3'); 
fprintf("\n%s from %s saved.\n",fname,subject)

end % end of function