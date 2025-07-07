%--------------------------------------------------------------------------
% Till Habersetzer, 01.07.2025
% Communication Acoustics, CvO University Oldenburg
% till.habersetzer@uol.de 
% 
% Description:
%   Analyzes and visualizes MEG cross-correlation results, focusing on the
%   statistical comparison between an actual condition and a shuffled null
%   condition.
%
%   Key Steps:
%   -   Loads pre-computed, subject-level and grand-average cross-correlation
%       data for both conditions.
%   -   Performs a cluster-based permutation t-test to identify significant
%       time-channel clusters where the conditions differ.
%   -   Generates a suite of visualizations, including:
%       - Topographical plots of grand-average waveforms.
%       - Topographical maps of statistical results (t-values, clusters).
%       - Time-series plots of selected channels showing individual and
%         grand-average data.
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
cfg.statistic               = 'ft_statfun_depsamplesT';
cfg.correctm                = 'cluster';
cfg.clusteralpha            = 0.05;
cfg.clusterstatistic        = 'maxsum';
cfg.minnbchan               = 3; % megmag = 2, megplanar 3, meg = 4 
cfg.tail                    = 0;
cfg.clustertail             = 0;
cfg.alpha                   = 0.025;
cfg.numrandomization        = 5000;
cfg.neighbours              = neighbours;
design                      = zeros(2,2*n_subj);
design(1,1:n_subj)          = 1;
design(1,n_subj+1:2*n_subj) = 2;
design(2,1:n_subj)          = 1:n_subj;
design(2,n_subj+1:2*n_subj) = 1:n_subj;

cfg.design = design;
cfg.ivar   = 1; % row of design matrix that contains independent variable 
cfg.uvar   = 2; % row of design matrix that contains unit of observation

stat = ft_timelockstatistics(cfg,avgs_crosscorr{:},avgs_crosscorr_shuffled{:});

% Show clusters
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

%% Visualize Statistic
%--------------------------------------------------------------------------

% Set alpha-value for mask
%-------------------------
% alpha_mask = 0.025;
alpha_mask = 0.001;
mask       = give_stat_mask(stat,alpha_mask);
stat.mask  = mask;

% Plot mask
% imagesc(mask);
% colormap(gray);

% Time windows for topoplots
%---------------------------
% windows = [-0.1,0;0,0.1;0.1,0.2;0.2,0.3;0.3,0.4;0.4,0.5;0.5,0.6;...
%            0.6,0.7;0.7,0.8;0.8,0.9];
windows = [-0.1,-0.05;-0.05,0;0,0.05;0.05,0.1;0.1,0.15;0.15,0.2;0.2,0.25;0.25,0.3;0.3,0.35;0.35,0.4;...
           0.4,0.5;0.5,0.6;0.6,0.7;0.7,0.8;0.8,0.9];     
idx     = dsearchn(stat.time',reshape(windows,1,numel(windows))');
idx     = reshape(idx,size(windows,1),size(windows,2));
maxval  = max(abs(stat.stat),[],'all');

% Plot over statistic
%--------------------
figure('Color', 'w')
for p_idx = 1:15
   subplot(3,5,p_idx);
   cfg      = [];
   cfg.xlim = windows(p_idx,:);  
   cfg.zlim = [-maxval,maxval];
   cfg.highlight        = 'on';
   chan2highlight       = all(stat.mask(:,idx(p_idx,1):idx(p_idx,2)), 2); % all samples must be significant during timeinterval
   % less severe
   % chan2highlight       = any(stat.mask(:,idx(p_idx,1):idx(p_idx,2)), 2); 
   cfg.highlightchannel = find(chan2highlight);
   cfg.highlightcolor   = 'k';
   cfg.highlightsize    = 50;
   cfg.highlightsymbol  = '.';
   cfg.comment          = 'xlim';
   cfg.commentpos       = 'title';
   cfg.layout           = layout;
   cfg.interactive      = 'no';
   cfg.parameter        = 'stat';
   cfg.figure           = gcf;
   ft_topoplotER(cfg, stat);
   set(gca,'ColorScale','linear') % if you want to plot on a logarithmic scale
   set(gca,'fontsize', 12)
   colormap(bluewhitered)
end
posi        = get(subplot(3,5,10),'Position'); %[left bottom width height]
posi(1)     = posi(1)+posi(3)+0.01; posi(3) = 0.25*posi(3);
c           = colorbar('Position',posi);
c.LineWidth = 1;
c.FontSize  = 15;
title(c,'t-Wert','fontweight','bold');
% set(c,'Ticks',[0.1,1,2,3,4,5,8],'TickLabels',{'0.1','1','2','3','4','5','8'})
sgtitle(sprintf('%s: t-values',sensortype),'fontweight','bold','FontSize',20)

% Crosscorrelation timeseries with statistic
%-------------------------------------------
cfg               = [];
cfg.layout        = 'neuromag306all_helmet.mat';
cfg.parameter     = 'stat';
cfg.maskparameter = 'mask';
cfg.graphcolor    = 'r';
figure('Color', 'w')
ft_multiplotER(cfg,stat);

sgtitle(sprintf('%s: t-values',sensortype),'fontweight','bold','FontSize',20)

%% Visualize cross correlation timeseries of both conditions
%--------------------------------------------------------------------------

% Compute standard error of mean 
gavg_crosscorr_sem          = sqrt(gavg_crosscorr.var) ./ sqrt(gavg_crosscorr.dof);
gavg_crosscorr_shuffled_sem = sqrt(gavg_crosscorr_shuffled.var) ./ sqrt(gavg_crosscorr_shuffled.dof);

% Plot selected channels for grand average (order not important) with 
% standard error of mean (sem) and individual time courses (optional)
%--------------------------------------------------------------------------

% channels to plot
chan2plot              = {'MEG0341','MEG0231','MEG1611','MEG1221','MEG1341','MEG2421'};
n_chan                 = length(chan2plot);
error_alpha            = 0.2; % Shaded error area transparency
timevec                = gavg_crosscorr.time*1000;
subject_alpha          = 0.6;
subject_color          = [0.8, 0.8, 0.8]; 
subject_color_shuffled = [0.5, 0.5, 0.5]; 
axis_font_size         = 18; 
legend_font_size       = 16;
title_font_size        = 18; 

xlims  = [-100,900];
idx    = find(contains(gavg_crosscorr.label,chan2plot));
minval = min(gavg_crosscorr.avg(idx,:)-gavg_crosscorr_sem(idx,:),[],'all');
maxval = max(gavg_crosscorr.avg(idx,:)+gavg_crosscorr_sem(idx,:),[],'all');

figure('Color', 'w')
axis_handles = gobjects(1, n_chan); % Pre-allocate an array for axis handles
for ch_idx = 1:n_chan

    chan_name = chan2plot{ch_idx};
    idx       = find(contains(gavg_crosscorr.label,chan_name));

    subplot(3,2,ch_idx);
    axis_handles(ch_idx) = gca; % Store the current axes handle
    hold on

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

    xlabel('t / ms')
    ylabel('$\hat{R}$','Interpreter','Latex')
    title(chan_name, 'FontSize', title_font_size)

    xticks(xlims(1):200:xlims(2)); % Set x-axis ticks in steps of 100 ms
    set(gca, 'FontSize', axis_font_size); % Set font size for axis values
    
    grid on;
    grid minor;
    box on; 

    % Add legend only for the first subplot
    % if ch_idx == 1
        legend([plt1, plt2],{'shuffled','sorted'},'Interpreter', 'none','Location','Southeast','FontSize',legend_font_size)
    % end
end
% Link the X-axes of all subplots
linkaxes(axis_handles, 'xy');
% Setting it on one linked axis will propagate to all others
ylims  = [minval, maxval]; % Centralized Y-axis limits
set(axis_handles, 'XLim', xlims, 'YLim', ylims);

% Highlight channels in topoplot
%-------------------------------
figure('Color', 'w')
cfg                  = [];
cfg.comment          = 'xlim';
cfg.xlim             = [0.06,0.1]; %[0.1,0.2]; ; %[0.095,0.125]; [0.08,0.08];
cfg.highlight        = 'on';
cfg.highlightchannel  = chan2plot;
cfg.highlightcolor   = 'k';
cfg.highlightsize    = 60;
cfg.highlightsymbol  = '.';
cfg.layout           = 'neuromag306all_helmet.lay';
cfg.interactive      = 'no';
cfg.colorbar         = 'yes';
ft_topoplotER(cfg,gavg_crosscorr_shuffled);
colormap(bluewhitered)

figure('Color', 'w')
cfg             = [];
cfg.showlabels  = 'yes';
cfg.fontsize    = 6;
cfg.layout      = 'neuromag306mag_helmet.lay';
% cfg.layout     = 'neuromag306planar_helmet.lay';
cfg.linestyle   = {'-','-'};
cfg.linecolor   = 'rb'; % shuffled: red, sorted: blue
cfg.showlabels  = 'off';
cfg.showcomment = 'off';
cfg.showscale   = 'off';
cfg.linewidth   = 1;
ft_multiplotER(cfg,gavg_crosscorr_shuffled,gavg_crosscorr);


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