%--------------------------------------------------------------------------
% Till Habersetzer, 01.07.2025
% Communication Acoustics, CvO University Oldenburg
% till.habersetzer@uol.de 
% 
% Description:
%   This script performs a comprehensive analysis and visualization of MEG
%   cross-correlation data. It focuses on comparing a primary experimental
%   condition ('sorted') against a shuffled null condition to identify
%   statistically significant differences. The analysis uses FieldTrip's
%   cluster-based permutation statistics and generates a wide range of plots
%   to inspect the data at both the grand-average and single-subject level.
%
%   Key Features & Visualizations:
%   1.  Data Loading: Imports pre-computed cross-correlation results for
%       individual subjects and the grand average.
%   2.  Trial Count Summary: Plots a bar chart of the number of trials
%       retained for each subject.
%   3.  Waveform Visualization: Creates topographical plots of grand-average
%       waveforms for both magnetometers and gradiometers, overlaying the
%       sorted, shuffled, and difference waves.
%   4.  Statistical Analysis: Conducts a cluster-based permutation t-test
%       to find significant time-channel clusters where the sorted and
%       shuffled conditions differ.
%   5.  Statistical Visualization:
%       -   Generates topographical plots of t-values, highlighting
%           significant clusters within specific time windows.
%       -   Plots the full time-course of the statistical results.
%   6.  Detailed Channel Plots: For a selection of channels, it visualizes
%       the time-series showing the grand-average with SEM error bands and
%       the data from all individual subjects.
%   7.  Butterfly Plots: Overlays the waveforms from all sensors to provide
%       a global view, including a version where channels are color-coded
%       by their anatomical region (e.g., temporal, frontal).
%--------------------------------------------------------------------------

close all
clearvars
clc 

%% Import main settings 
%--------------------------------------------------------------------------
current_dir = pwd;
cd(fullfile('..'))
settings_speech
cd(current_dir)

% Addpath for additional functions
addpath(fullfile(settings.path2project,'analysis','helper_functions'))

%% Script settings
%--------------------------------------------------------------------------
% Choose subject for plotting
subjects = 1:24;
n_subj   = length(subjects);

%% Import data
%--------------------------------------------------------------------------
avgs_crosscorr          = cell(1,n_subj);
avgs_crosscorr_shuffled = cell(1,n_subj);
n_trials                = zeros(1,n_subj);
subjectnames            = cell(1,n_subj);

for sub_idx = 1:n_subj
    subject               = sprintf('sub-%02d',subjects(sub_idx));
    subjectnames{sub_idx} = subject;

    data                             = importdata(fullfile(settings.path2derivatives,subject,'speech',sprintf('%s_crosscorr.mat',subject)));   
    avgs_crosscorr{sub_idx}          = data.avg_crosscorr;
    avgs_crosscorr_shuffled{sub_idx} = data.avg_crosscorr_shuffled;
    n_trials(sub_idx)                = data.n_trials;
    clear data
    fprintf('%s loaded.\n',subject)

end

% Import grand average over all subjects
%---------------------------------------
subject                 = 'grandaverage';
data                    = importdata(fullfile(settings.path2derivatives,subject,'speech',sprintf('%s_crosscorr.mat',subject)));   
gavg_crosscorr          = data.gavg_crosscorr;
gavg_crosscorr_shuffled = data.gavg_crosscorr_shuffled;
n_subj_gavg             = data.n_subjects;
clear data

if ~isequal(n_subj_gavg,n_subj)
    error('Unexpected number of subjects for group analysis! Recompute grandaverage!')
end

% Compute raweffect
%------------------
cfg                      = [];
cfg.operation            = 'subtract';
cfg.parameter            = 'avg';
gavg_crosscorr_raweffect = ft_math(cfg,gavg_crosscorr,gavg_crosscorr_shuffled);

%% Plot kept trials
%--------------------------------------------------------------------------

figure; 
b1 = bar(1:n_subj, n_trials, 'FaceColor', [0.8 0.8 0.8]); % Light grey
ax                    = gca; 
ax.XTick              = 1:n_subj; 
ax.XTickLabel         = subjectnames;
ax.XTickLabelRotation = 45; 
title(sprintf('Number of trials (Trialduration: %is)',settings.crosscorr.trialdur));
ylabel('Number of Trials');
xlabel('Subject');
grid on; 

%% Plot crosscorrelations for all sensors arranged topographically 
%--------------------------------------------------------------------------
figure
cfg            = [];
cfg.showlabels = 'yes';
cfg.fontsize   = 6;
cfg.layout     = 'neuromag306mag_helmet.lay';
cfg.linestyle  = {'--','-.','-'};
cfg.linecolor  = 'rmb';
ft_multiplotER(cfg, gavg_crosscorr,gavg_crosscorr_shuffled,gavg_crosscorr_raweffect);
title('gavg: magnetometer')

% Gradiometer
%------------
figure
cfg            = [];
cfg.showlabels = 'yes';
cfg.fontsize   = 6;
cfg.layout     = 'neuromag306planar_helmet.lay';
cfg.linestyle  = {'--','-.','-'};
cfg.linecolor  = 'rmb';
ft_multiplotER(cfg, gavg_crosscorr,gavg_crosscorr_shuffled,gavg_crosscorr_raweffect);
title('gavg: gradiometer')

%% Calculate cluster-based permutation statistic and show clusters
%--------------------------------------------------------------------------
% Within-subjects experiment (within-UO design)
% null hypothesis: the probability distribution of the condition-specific averages 
%                  is independent of the experimental conditions.

% choose sensors
%---------------
% you may need to change alpha value and cfg.minnbchan 
sensortype = 'megmag';
cfg        = [];

switch sensortype
    case 'megmag'
        layout       = 'neuromag306mag_helmet.mat';
        cfg.channel  = sensortype;
        cfg.method   = 'template';
        cfg.template = 'neuromag306mag_neighb.mat';
        neighbours   = ft_prepare_neighbours(cfg);
    case 'megplanar'
        layout       = 'neuromag306planar_helmet.mat';
        cfg.channel  = sensortype;
        cfg.method   = 'template';
        cfg.template = 'neuromag306planar_neighb.mat';
        neighbours   = ft_prepare_neighbours(cfg);
end
cfg.neighbours = neighbours;
cfg.grad       = avgs_crosscorr{1}.grad;
ft_neighbourplot(cfg) 
 
cfg                         = [];
cfg.channel                 = sensortype;
cfg.method                  = 'montecarlo';
cfg.statistic               = 'depsamplesT'; % the dependent samples T-statistic (within-UO)
cfg.correctm                = 'cluster'; %  MCP by calculating a so-called cluster-based test statistic 
cfg.clusteralpha            = 0.05; % alpha level of the sample-specific test statistic that will be used for thresholding
cfg.clusterstatistic        = 'maxsum';
cfg.minnbchan               = 2; % megmag = 2, megplanar 3, meg = 4 / % minimum number of neighborhood channels that is required for a selected sample to be included
cfg.tail                    = 0; % -1, 1 or 0 (default = 0); one-sided or two-sided test
cfg.clustertail             = 0;
cfg.alpha                   = 0.025; % alpha level of the permutation test
cfg.numrandomization        = 5000;
cfg.neighbours              = neighbours;
design                      = zeros(2,2*n_subj);
design(1,:)                 = [1:n_subj 1:n_subj];
design(2,:)                 = [ones(1,n_subj) ones(1,n_subj)*2];

cfg.design = design;
cfg.uvar   = 1; % row of design matrix that contains unit of observation
cfg.ivar   = 2; % row of design matrix that contains independent variable 

stat = ft_timelockstatistics(cfg,avgs_crosscorr{:},avgs_crosscorr_shuffled{:});

% Show a specific cluster 
%--------------------------------------------------------------------------
cluster_number = 1;
pos            = any(ismember(stat.posclusterslabelmat,cluster_number),2);
neg            = any(ismember(stat.negclusterslabelmat,cluster_number),2);
% stat.selected_cluster = pos;
stat.selected_cluster = neg;

figure('Color', 'w');
cfg                  = [];
cfg.xlim             = [0.06,0.1];
cfg.parameter        = 'stat';
cfg.comment          = 'no';
cfg.highlight        = 'on';
cfg.highlightchannel = stat.label(stat.selected_cluster);
cfg.highlightcolor   = 'k';
cfg.highlightsize    = 60;
cfg.highlightsymbol  = '.';
cfg.layout           = layout;
cfg.interactive      = 'no';
ft_topoplotER(cfg,stat);
colormap(bluewhitered)

%% Visualize Statistic over time
%--------------------------------------------------------------------------

% Set alpha-value for mask
%-------------------------
alpha_mask = 0.025;
mask       = give_stat_mask(stat,alpha_mask);
stat.mask  = mask;

% Plot mask
% imagesc(mask);
% colormap(gray);

% Time windows for topoplots
%---------------------------
timewindows_statistic = [50, 90; 160, 230]/1000; % sec

time_idx = dsearchn(stat.time',reshape(timewindows_statistic,1,numel(timewindows_statistic))');
time_idx = reshape(time_idx,size(timewindows_statistic,1),size(timewindows_statistic,2));

minval = min(gavg_crosscorr_raweffect.avg(:));
maxval = max(gavg_crosscorr_raweffect.avg(:));

% Plot over statistic
%--------------------
% First ensure the channels to have the same order in the average and in the statistical output.
% This might not be the case, because ft_math might shuffle the order
[i1,i2] = match_str(gavg_crosscorr_raweffect.label, stat.label);

figure('Color', 'w')
% Pre-allocate an array to store the handles of the subplots
ax = gobjects(1, size(timewindows_statistic,1));

for p_idx = 1:size(timewindows_statistic,1) 

   chan2highlight       = zeros(numel(gavg_crosscorr_raweffect.label),1);
   chan2highlight(i1)   = all(stat.mask(i2, time_idx(p_idx,1):time_idx(p_idx,2)), 2);
   % less strict
   % chan2highlight(i1)   = any(stat.mask(i2,time_idx(p_idx,1):time_idx(p_idx,2)), 2); 

   % Store the handle of the current subplot
   ax(p_idx)            = subplot(1,2,p_idx);

   cfg                  = [];
   cfg.xlim             = timewindows_statistic(p_idx,:);  
   cfg.zlim             = [minval,maxval];
   cfg.highlight        = 'on';
   cfg.highlightchannel = find(chan2highlight);
   cfg.highlightcolor   = 'k';
   cfg.highlightsize    = 50;
   cfg.highlightsymbol  = '.';
   cfg.comment          = 'xlim';
   cfg.commentpos       = 'title';
   cfg.layout           = layout;
   cfg.interactive      = 'no';
   cfg.figure           = gca; % Use the current axes handle
   ft_topoplotER(cfg, gavg_crosscorr_raweffect);
   colormap(bluewhitered)
end

% Get the position of the second subplot using its stored handle
pos2        = get(ax(2), 'Position'); % [left, bottom, width, height]
c           = colorbar;
c.Position  = [pos2(1) + pos2(3) + 0.02, pos2(2), 0.03, pos2(4)]; % [left, bottom, width, height]
c.LineWidth = 1;
c.FontSize  = 15;
title(c, '$\hat{R}$', 'fontweight', 'bold', 'Interpreter', 'latex');

sgtitle(sprintf('%s: Correlation-values condition contrast', sensortype), 'fontweight', 'bold', 'FontSize', 20);

% Crosscorrelation timeseries with statistic
%-------------------------------------------
cfg               = [];
cfg.layout        = 'neuromag306all_helmet.mat';
cfg.parameter     = 'stat';
cfg.maskparameter = 'mask';
cfg.graphcolor    = 'r';
figure('Color', 'w')
ft_multiplotER(cfg,stat);

sgtitle(sprintf('%s: Cluster-based permutation test ',sensortype),'fontweight','bold','FontSize',20)

%% Visualize cross correlation timeseries of both conditions
%--------------------------------------------------------------------------

% Compute standard error of mean 
gavg_crosscorr_sem          = sqrt(gavg_crosscorr.var) ./ sqrt(gavg_crosscorr.dof);
gavg_crosscorr_shuffled_sem = sqrt(gavg_crosscorr_shuffled.var) ./ sqrt(gavg_crosscorr_shuffled.dof);

% Plot selected channels for grand average (order not important) with 
% standard error of mean (sem) and individual time courses (optional)
%--------------------------------------------------------------------------

% channels to plot
chan2plot              = {'MEG0341','MEG1221','MEG0231','MEG1341','MEG1611','MEG2421'};
n_chan                 = length(chan2plot);
error_alpha            = 0.2; % Shaded error area transparency
timevec                = gavg_crosscorr.time*1000;
subject_alpha          = 0.6;
subject_color          = [0.8, 0.8, 0.8]; 
subject_color_shuffled = [0.5, 0.5, 0.5]; 
axis_font_size         = 16; 
legend_font_size       = 18;
title_font_size        = 14; 

xlims  = [-100,900];
idx    = find(contains(gavg_crosscorr.label,chan2plot));
minval = min(gavg_crosscorr.avg(idx,:)-gavg_crosscorr_sem(idx,:),[],'all');
maxval = max(gavg_crosscorr.avg(idx,:)+gavg_crosscorr_sem(idx,:),[],'all');
% ylims  = [-0.02,0,0.02];

figure('Color', 'w')
t = tiledlayout(3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
axis_handles = gobjects(1, n_chan); % Pre-allocate an array for axis handles
for ch_idx = 1:n_chan

    chan_name = chan2plot{ch_idx};
    idx       = find(contains(gavg_crosscorr.label,chan_name));

    % Go to the next tile (subplot) 
    axis_handles(ch_idx) = nexttile;
    hold on;

    % Plot individual subjects
    for sub_idx = 1:n_subj
        plot(timevec, avgs_crosscorr_shuffled{sub_idx}.avg(idx,:),'LineWidth', 1,'Color', [subject_color_shuffled,subject_alpha], 'LineStyle', '-','HandleVisibility', 'off');
        plot(timevec, avgs_crosscorr{sub_idx}.avg(idx,:),'LineWidth', 1,'Color', [subject_color,subject_alpha], 'LineStyle', '-','HandleVisibility', 'off');
        % update minval, maxval
        % minval = min([minval,min([avgs_crosscorr_shuffled{sub_idx}.avg(idx,:),avgs_crosscorr{sub_idx}.avg(idx,:)])]);
        % maxval = max([maxval,max([avgs_crosscorr_shuffled{sub_idx}.avg(idx,:),avgs_crosscorr{sub_idx}.avg(idx,:)])]);
    end

    ft_plot_vector(timevec, [gavg_crosscorr_shuffled.avg(idx,:) + gavg_crosscorr_shuffled_sem(idx,:); gavg_crosscorr_shuffled.avg(idx,:) - gavg_crosscorr_shuffled_sem(idx,:)], 'highlightstyle', 'difference', 'facealpha', error_alpha, 'facecolor', 'r')
    plt1 = ft_plot_vector(timevec, gavg_crosscorr_shuffled.avg(idx,:), 'color', 'r', 'linewidth', 2, 'style', '--');

    ft_plot_vector(timevec, [gavg_crosscorr.avg(idx,:) + gavg_crosscorr_sem(idx,:); gavg_crosscorr.avg(idx,:) - gavg_crosscorr_sem(idx,:)], 'highlightstyle', 'difference', 'facealpha', error_alpha, 'facecolor', 'b')
    plt2 = ft_plot_vector(timevec, gavg_crosscorr.avg(idx,:), 'color', 'b', 'linewidth', 2, 'style', '-');

    if ch_idx > 4
        xlabel('t / ms')
    end
    ylabel('$\hat{R}$','Interpreter','Latex')
    title(chan_name, 'FontSize', title_font_size)

    xticks(xlims(1):200:xlims(2)); % Set x-axis ticks in steps of 100 ms
    % yticks(ylims);
    set(gca, 'FontSize', axis_font_size); % Set font size for axis values
    
    grid on;
    grid minor;
    box on; 

    % Add legend only for the first subplot
    if ch_idx == 1
        legend([plt1, plt2],{'shuffled','sorted'},'Interpreter', 'none','Location','Southeast','FontSize',legend_font_size)
    end

    if ch_idx <= 4
        xticklabels([]); % Hides the numbers, but keeps grid lines 
    end
end
% Link the X-axes of all subplots
linkaxes(axis_handles, 'xy');
% Setting it on one linked axis will propagate to all others
ylims  = [minval, maxval]; % Centralized Y-axis limits
set(axis_handles, 'XLim', xlims, 'YLim', ylims);

% Cross-correlation timeseries
%-----------------------------
figure('Color', 'w')
cfg             = [];
cfg.showlabels  = 'yes';
cfg.fontsize    = 6;
cfg.layout      = 'neuromag306mag_helmet.lay';
% cfg.layout     = 'neuromag306planar_helmet.lay';
cfg.linestyle   = {'-','-'};
cfg.linecolor   = 'rb'; % shuffled: red, sorted: blue
cfg.showlabels  = 'on';
cfg.showcomment = 'off';
cfg.showscale   = 'off';
cfg.linewidth   = 1;
ft_multiplotER(cfg,gavg_crosscorr_shuffled,gavg_crosscorr);

%% Visualize cross correlations
%--------------------------------------------------------------------------

sensortype = 'megmag';
% sensortype = 'megplanar';

switch sensortype
    case 'megmag'
        is_sensor  = endsWith(gavg_crosscorr.label, '1');
        is_sensor2 = endsWith(gavg_crosscorr_shuffled.label, '1');

    case 'megplanar'
        is_sensor  = endsWith(gavg_crosscorr.label, '2') | endsWith(gavg_crosscorr.label, '3');
        is_sensor2 = endsWith(gavg_crosscorr_shuffled.label, '2') | endsWith(gavg_crosscorr_shuffled.label, '3');
end
% Check if the two logical arrays are identical. If not, throw an error.
if ~isequal(is_sensor, is_sensor2)
    error('Mismatch in sensor labels: The order or type of channels in the two datasets does not match.');
end

n_channels = sum(is_sensor);
timevec    = gavg_crosscorr.time*1000; 
xlims      = [-100,900];

axis_font_size   = 16; 
legend_font_size = 16;
title_font_size  = 18; 

% All channels same color
%--------------------------------------------------------------------------
figure;
h1 = plot(timevec, gavg_crosscorr.avg(is_sensor,:), 'b', 'LineWidth', 1);
hold on;
h2 = plot(timevec, gavg_crosscorr_shuffled.avg(is_sensor,:), 'r', 'LineWidth', 1);
hold off; 
xlabel('t / ms')
ylabel('$\hat{R}$','Interpreter','Latex')
title('Comparison of Original and Shuffled Cross-Correlations');
xticks(xlims(1):200:xlims(2)); % Set x-axis ticks in steps of 100 ms
set(gca, 'FontSize', axis_font_size); % Set font size for axis values

legend([h1(1), h2(1)], {'Sorted', 'Shuffled'}, 'Location', 'Southeast');
grid on;
grid minor;
box on; 

% Color grouped by sensor location
%--------------------------------------------------------------------------
channel_layout_meg    = importdata('channel_layout_meg.mat');
channel_subsets       = cell(1,4);
channel_subsets{4}    = unique([channel_layout_meg.occipital_left(:); channel_layout_meg.occipital_right(:)]); % occipital_channels
channel_subsets{3}    = unique([channel_layout_meg.parietal_left(:); channel_layout_meg.parietal_right(:)]); % parietal_channels
channel_subsets{2}    = unique([channel_layout_meg.frontal_left(:); channel_layout_meg.frontal_right(:)]); % frontal channels
channel_subsets{1}    = unique([channel_layout_meg.temporal_left(:); channel_layout_meg.temporal_right(:)]); % temporal channels
channel_subsets_label = {'Sorted (temporal)','Sorted (frontal)','Sorted (parietal)','Sorted (occipital)'};

cmap = [
    0.0000    0.2980    0.5294;    % Dark Blue (#004C87)
    0.2000    0.7451    0.9333;    % Cyan (#33BEE7)
    0.4660    0.6740    0.1880;    % Green (#77AB30)
    0.9290    0.5940    0.1920     % Orange (#ED9731)
];

line_styles = {'-','--',':','-.'};
% line_styles = {'-','-','-','-'};

% Compute limits
%---------------
all_data = [gavg_crosscorr.avg(is_sensor,:); gavg_crosscorr_shuffled.avg(is_sensor,:)];

% 2. Find the absolute min and max across all the data
min_val = min(all_data(:));
max_val = max(all_data(:));
clear all_data

data_range = max_val - min_val;
padding    = 0.05 * data_range;
yLims      = [min_val - padding, max_val + padding];

figure;
subplot(2,1,1)
hold on;

% Pre-allocate a vector to store plot handles for the legend
h = gobjects(1, length(channel_subsets) + 2);

for loc_idx = 1:length(channel_subsets)

    subset_idx = contains(gavg_crosscorr.label, channel_subsets{loc_idx});
  
    chan_idx     = and(is_sensor,subset_idx); % channel type + subsets
    plot_handles = plot(timevec, gavg_crosscorr.avg(chan_idx,:), 'color', cmap(loc_idx,:), 'LineWidth', 1, 'LineStyle', line_styles{loc_idx});

    % Save only the FIRST handle from this group for the legend
    if ~isempty(plot_handles)
        h(loc_idx) = plot_handles(1);
    end

end
plot_handles = plot(timevec, gavg_crosscorr_shuffled.avg(is_sensor,:), 'color', 'r', 'LineWidth', 1);
h(length(channel_subsets)+1) = plot_handles(1);

% Define patch properties
patchColor = [1, 1, 0]; % A light, transparent yellow
patchAlpha = 0.2;

% Create the first patch for the 50-90ms window
p1 = patch([timewindows_statistic(1,1), timewindows_statistic(1,2), timewindows_statistic(1,2), timewindows_statistic(1,1)]*1000, [yLims(1), yLims(1), yLims(2), yLims(2)], ...
           patchColor, 'FaceAlpha', patchAlpha, 'EdgeColor', 'none');
h(end) = p1;
% Create the second patch for the 160-230ms window
p2 = patch([timewindows_statistic(2,1), timewindows_statistic(2,2), timewindows_statistic(2,2), timewindows_statistic(2,1)]*1000, [yLims(1), yLims(1), yLims(2), yLims(2)], ...
           patchColor, 'FaceAlpha', patchAlpha, 'EdgeColor', 'none');
           
% Send the patches to the background, behind the data lines
uistack([p1, p2], 'bottom');

hold off; 
xlabel('t / ms')
ylabel('$\hat{R}$','Interpreter','Latex')
title('Butterfly plot cross-correlation functions');
xticks(xlims(1):200:xlims(2)); % Set x-axis ticks in steps of 100 ms
set(gca, 'FontSize', axis_font_size); % Set font size for axis values
ylim(yLims)

legend(h, [channel_subsets_label, 'Shuffled', 'time window statistic'], 'Location', 'Southeast', 'NumColumns', 5);
grid on;
grid minor;
box on; 

% Create dummy data for mapping color coding of sensors
%------------------------------------------------------
dummy_data      = gavg_crosscorr_raweffect;
dummy_data.time = 0;
dummy_data.avg  = zeros(306,1);

for loc_idx = 1:length(channel_subsets)

    subset_idx                 = contains(dummy_data.label, channel_subsets{loc_idx});
    dummy_data.avg(subset_idx) = loc_idx;

end

% Highlight channels in topoplot
%-------------------------------
fig = figure('Color', 'w');
hold on
cfg                  = [];
cfg.comment          = 'no';
cfg.xlim             = [0,0]; 
cfg.markersize       = 20;
cfg.markersymbol     = '.';
cfg.highlight        = 'on';
cfg.highlightchannel  = chan2plot;
cfg.highlightcolor   = 'k';
cfg.highlightsize    = 50;
cfg.highlightsymbol  = '.';
cfg.layout           = 'neuromag306mag_helmet.lay';
cfg.interactive      = 'no';
cfg.colorbar         = 'yes';
cfg.interpolation    = 'nearest';
cfg.figure           = fig;
ft_topoplotER(cfg,dummy_data); % plot difference
colormap(cmap)


%% Functions
%--------------------------------------------------------------------------

function mask = give_stat_mask(stat,alpha)
%--------------------------------------------------------------------------
% Till Habersetzer, 01.07.2025
% Communication Acoustics, CvO University Oldenburg
% till.habersetzer@uol.de 
%
%   MASK = GIVE_STAT_MASK(STAT, ALPHA) takes the output of a FieldTrip
%   cluster-based permutation analysis (STAT) and a desired alpha-value.
%   It returns a binary mask (MASK) where significant clusters (p < ALPHA)
%   are marked with 'true' and non-significant areas with 'false'.
%   This mask can be used for visualizing statistically significant regions.
%
%   Inputs:
%     stat  - Structure array containing the results of ft_timelockstatistics
%             or ft_sourcestatistics, including cluster information.
%     alpha - Desired significance level (e.g., 0.05) for thresholding the clusters.
%
%   Output:
%     mask  - A logical array (binary mask) indicating the significant
%             clusters. True for significant, false for non-significant.
% 
% Example Usage for plotting the mask:
% Assuming 'stat' is the result of a FieldTrip cluster analysis
% and 'alpha' is your significance level (e.g., 0.05)
% mask = give_stat_mask(stat, alpha);
% figure;
% imagesc(mask);
% colormap(gray); % or any suitable colormap for binary masks
%--------------------------------------------------------------------------

% positive clusters
%------------------
pos_cluster_pvals = [stat.posclusters(:).prob];
pos_signif_clust  = find(pos_cluster_pvals < alpha);
pos               = ismember(stat.posclusterslabelmat, pos_signif_clust);
% negative clusters
%------------------
neg_cluster_pvals = [stat.negclusters(:).prob];
neg_signif_clust  = find(neg_cluster_pvals < alpha);
neg               = ismember(stat.negclusterslabelmat, neg_signif_clust);
% combine clusters
%-----------------
mask              = or(neg,pos);

end