%--------------------------------------------------------------------------
% Till Habersetzer, 20.06.2025
% Communication Acoustics, CvO University Oldenburg
% till.habersetzer@uol.de 
%
% This script generates the final summary figures for the AEF sensor level
% results, aggregating results across all participants. 
% It loads the pre-computed time-locked data for each subject and the grand 
% average to create a series of publication-quality visualizations. 
% The script produces a data quality summary (kept trials per subject), 
% FieldTrip plots of the grand average (ft_multiplotER, ft_topoplotER), and 
% a detailed analysis of selected channels. This final analysis includes an 
% overlay of all individual subject responses on the grand average 
% (butterfly plot) and a corresponding plot of the grand average with the 
% standard error, providing a comprehensive overview of the results and 
% their variability.
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

% Addpath for additional functions
addpath(fullfile(settings.path2project,'analysis','helper_functions'))

%% Script settings
%--------------------------------------------------------------------------
% Choose subject for plotting
subjects = 2:24;
n_subj   = length(subjects);

%% Import data
%--------------------------------------------------------------------------
avgs_mag      = cell(1,n_subj);
avgs_grad     = cell(1,n_subj);
avgs_cmb      = cell(1,n_subj);
n_trials_mag  = zeros(1,n_subj);
n_trials_grad = zeros(1,n_subj);
n_trials_tot  = zeros(1,n_subj);
subjectnames  = cell(1,n_subj);

for sub_idx = 1:n_subj
    subject = sprintf('sub-%02d',subjects(sub_idx));

    data                   = importdata(fullfile(settings.path2project,'derivatives',subject,'chirp',[subject,'_avgs.mat']));   
    avgs_grad{sub_idx}     = data.avg_grad;
    avgs_mag{sub_idx}      = data.avg_mag;
    n_trials_tot(sub_idx)  = sum(data.n_trials(1,:));
    n_trials_grad(sub_idx) = data.n_trials_grad(2);
    n_trials_mag(sub_idx)  = data.n_trials_mag(2);
    subjectnames{sub_idx}  = subject;
    fprintf('%s loaded.\n',subject)
 
    % Add combined gradiometers
    cfg               = [];
    cfg.method        = 'sum';
    avgs_cmb{sub_idx} = ft_combineplanar(cfg,avgs_grad{sub_idx});
    clear data

    % Convert Units of magnetometers T -> fT
    avgs_mag{sub_idx}.avg = 10^15*avgs_mag{sub_idx}.avg;

end

% Computa grand average over all subjects
%----------------------------------------
% cfg         = [];
% cfg.latency = 'all';
% gavg_grad   = ft_timelockgrandaverage(cfg,avgs_grad{:});
% gavg_mag    = ft_timelockgrandaverage(cfg,avgs_mag{:});

% Import grand average data
%--------------------------
subject   = 'grandaverage';
data      = importdata(fullfile(settings.path2project,'derivatives','grandaverage','chirp',[subject,'_avgs.mat']));   
gavg_grad = data.gavg_grad;
gavg_mag  = data.gavg_mag;
clear data
% Convert Units of magnetometers T -> fT
gavg_mag.avg = 10^15*gavg_mag.avg;

% Add combined gradiometers
% cfg        = [];
% cfg.method = 'sum';
% gavg_cmb   = ft_combineplanar(cfg,gavg_grad);

cfg         = [];
cfg.latency = 'all';
gavg_cmb    = ft_timelockgrandaverage(cfg,avgs_cmb{:});

%% Plot kept trials
%--------------------------------------------------------------------------

figure; 
b1 = bar(1:n_subj, [n_trials_tot(:), n_trials_tot(:)], 'grouped', 'FaceColor', [0.8 0.8 0.8]); % Light grey
hold on
b2 = bar(1:n_subj, [n_trials_grad(:), n_trials_mag(:)], 'grouped','FaceColor', 'flat');
ax                    = gca; 
ax.XTick              = 1:n_subj; 
ax.XTickLabel         = subjectnames;
ax.XTickLabelRotation = 45; 
title('Number of kept trials after preprocessing');
ylabel('Number of Trials');
xlabel('Subject');
grid on; 
b1(1).DisplayName      = 'Before Rejection';
b1(2).HandleVisibility = 'off';
b2(1).DisplayName      = 'Gradiometers';
b2(2).DisplayName      = 'Magnetometers';
b2(1).FaceColor        = 'red';  % Gradiometer bars
b2(2).FaceColor        = 'blue'; % Magnetometer bars
legend('Location', 'best');

%% Plot of event related fields for all sensors arranged topographically 
%--------------------------------------------------------------------------

% Magnetometer
%-------------
figure
cfg            = [];
cfg.showlabels = 'yes';
cfg.fontsize   = 6;
cfg.layout     = 'neuromag306mag_helmet.lay';
ft_multiplotER(cfg, gavg_mag);
title('magnetometer')

% Gradiometer
%------------
figure
cfg            = [];
cfg.showlabels = 'yes';
cfg.fontsize   = 6;
cfg.layout     = 'neuromag306planar_helmet.lay';
ft_multiplotER(cfg, gavg_grad);
title('gradiometer')

% Combined Planar Gradiometers
%-----------------------------
figure
cfg            = [];
cfg.showlabels = 'yes';
cfg.fontsize   = 6;
cfg.layout     = 'neuromag306cmb_helmet.lay';
ft_multiplotER(cfg, gavg_cmb);
title('combined gradiometer')

%% Plot the topographic distribution of the data
%--------------------------------------------------------------------------
timewin = [0.09, 0.12]; % time window to look at in more detail

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
ft_topoplotER(cfg, gavg_mag);
title('Magnetometer');

cfg.layout   = 'neuromag306planar_helmet.mat';
cfg.figure   = subplot(1,3,2);
ft_topoplotER(cfg, gavg_grad);
title('Gradiometer');

cfg.layout   = 'neuromag306cmb_helmet.mat';
cfg.figure   = subplot(1,3,3);
ft_topoplotER(cfg, gavg_cmb);
title('Combined Gradiometer');

sgtitle(sprintf('%s: [%d-%d] ms',subject, timewin(1)*1000,timewin(2)*1000))

%% Final Plot
%--------------------------------------------------------------------------
% Select which sensor to plot
sens2plot = 'mag'; % 'mag' 'grad' 'cmb'
chan2plot = cell(1,2);

switch sens2plot
    case 'mag'
        gavg2plot    = gavg_mag;
        avgs2plot    = avgs_mag;
        layout       = 'neuromag306mag_helmet.mat';
        chan2plot{1} = 'MEG1621'; % left
        chan2plot{2} = 'MEG2411'; % right

        b_unit       = 'fT';
    case 'grad'
        gavg2plot = gavg_grad;
        avgs2plot = avgs_grad;
        % chan2plot{1} = 'MEG0242'; % left
        % chan2plot{2} = 'MEG1332'; % right
        chan2plot{1} = 'MEG0243'; % left
        chan2plot{2} = 'MEG1333'; % right
        layout    = 'neuromag306planar_helmet.mat';

        b_unit       = 'T / cm';
    case 'cmb'
        gavg2plot    = gavg_cmb;
        avgs2plot    = avgs_cmb;
        layout       = 'neuromag306cmb_helmet.mat';
        chan2plot{1} = 'MEG0242+0243'; % left
        chan2plot{2} = 'MEG1332+1333'; % right

        b_unit       = 'T / cm';
end

aefs    = zeros(2,length(gavg2plot.time),n_subj); % auditory evoked fields for subjects
gaefs   = zeros(2,length(gavg2plot.time)); % auditory evoked fields for grand average
cmap    = distinguishable_colors(n_subj);
minvals = zeros(1,n_subj);
maxvals = zeros(1,n_subj);

% For each subject
%-----------------
for sub_idx = 1:n_subj
    for ch_idx = 1:2 % 1:left / 2:right
        idx                    = find(contains(avgs2plot{sub_idx}.label,chan2plot{ch_idx}));
        aefs(ch_idx,:,sub_idx) = avgs2plot{sub_idx}.avg(idx,:);
    end   
    minvals(sub_idx) = min(aefs(:,:,sub_idx),[],'all');
    maxvals(sub_idx) = max(aefs(:,:,sub_idx),[],'all');
end

% For grandaverage
%-----------------
for ch_idx = 1:2
    idx             = find(contains(gavg2plot.label,chan2plot{ch_idx}));
    gaefs(ch_idx,:) = gavg2plot.avg(idx,:);
end
% Compute standard error
stderror = std(aefs,0,3)./sqrt(n_subj);

% Plot aefs
%--------------------------------------------------------------------------

% Data Parameters
timevec = gavg2plot.time * 1000; % Convert time to ms
ylims   = [min(minvals), max(maxvals)]; % Centralized Y-axis limits
xlims   = [-100, 400];

% Style Parameters
subject_line_width = 2;
subject_alpha      = 0.3; % Transparency for individual subject lines
ga_line_width      = 2;
ga_color           = 'k'; % Grand average line color
error_color        = [0.8, 0.8, 0.8]; % Shaded error area color
error_alpha        = 0.8; % Shaded error area transparency

title_font_size    = 18; % Font size for subplot titles
label_font_size    = 16; % Font size for x and y axis labels
legend_font_size   = 10;
legend_num_columns = 7;

% Create Figure and Axes
% Pre-define axes handles for clarity and control.
fig = figure('Name', 'AEFs','Color', 'w');
ax_top1 = subplot(2, 2, 1);
ax_top2 = subplot(2, 2, 2);
ax_bottom = subplot(2, 2, [3, 4]);

top_axes = [ax_top1, ax_top2]; % Group axes for easy looping

% Plot Individual Subjects and Grand Average (Top Subplots)
for i = 1:2
    % Set the current axes for plotting
    axes(top_axes(i));
    hold on;
    
    % Plot all subject waveforms with transparency
    for sub_idx = 1:n_subj
        plot(timevec, aefs(i, :, sub_idx), ...
            'LineWidth', subject_line_width, ...
            'Color', [cmap(sub_idx, :), subject_alpha]);
    end
    
    % Overlay the grand average waveform
    plot(timevec, gaefs(i, :), ...
        'LineWidth', ga_line_width, ...
        'Color', ga_color, ...
        'HandleVisibility', 'off'); % Exclude from legend
        
    % Apply common formatting
    xlabel('t / ms', 'FontSize', label_font_size);
    ylabel(sprintf('sig / %s',b_unit), 'FontSize', label_font_size)
    ylim(ylims);
    xlim(xlims);
    grid on;
    grid minor;
end
title(ax_top1,sprintf('Left Hemisphere: %s', chan2plot{1}), 'FontSize', title_font_size);
title(ax_top2,sprintf('Right Hemisphere: %s', chan2plot{2}), 'FontSize', title_font_size);

% Link Y-axes of top plots for consistent scaling
linkaxes(top_axes, 'y');

% Add a single, comprehensive legend to the first top subplot
lgd            = legend(ax_top1, subjectnames);
lgd.FontSize   = legend_font_size;
lgd.NumColumns = legend_num_columns;
lgd.Box = 'off';

% Plot Grand Average Comparison with Standard Error (Bottom Subplot)
% Use the pre-defined handle for the bottom axes.
axes(ax_bottom);
hold on;

% Define line styles for the two channels
ga_linestyles = {'-', '--'};
ga_legends    = cell(1, 2);

% Plot grand averages with shaded error for both channels
for i = 1:2
    % Create shaded error patch
    patch([timevec, fliplr(timevec)], ...
          [gaefs(i, :) + stderror(i, :), fliplr(gaefs(i, :) - stderror(i, :))], ...
          error_color, ...
          'FaceAlpha', error_alpha, ...
          'LineStyle', 'none', ...
          'HandleVisibility', 'off'); % Exclude from legend
          
    % Plot grand average line
    plot(timevec, gaefs(i, :), ...
        'LineWidth', ga_line_width, ...
        'Color', ga_color, ...
        'LineStyle', ga_linestyles{i});
end

% Apply formatting
xlabel('t / ms', 'FontSize', label_font_size);
ylabel(sprintf('sig / %s',b_unit), 'FontSize', label_font_size)
xlim(xlims);
grid on;
grid minor;
legend(chan2plot, 'Location', 'northwest');
hold off;

% Add topoplopt
%--------------------------------------------------------------------------
timewin2plot = zeros(1,2);
timewin(1)   = 0.105;
timewin(2)   = 0.180;

figure('Color', 'w');
cfg                  = [];
cfg.xlim             =  [timewin(1),timewin(1)];
cfg.highlight        = 'on';
cfg.highlightchannel = chan2plot;
cfg.highlightcolor   = 'k';
cfg.highlightsize    = 90;
cfg.highlightsymbol  = '.';
cfg.layout           = layout;
cfg.interactive      = 'no';
cfg.comment          = 'no';
ft_topoplotER(cfg,gavg2plot);
title(sprintf('%s: %i ms',sens2plot,timewin(1)*1000))
set(gca,'fontsize', 12)
colormap(bluewhitered)

cfg.xlim             =  [timewin(2),timewin(2)];
figure('Color', 'w');
ft_topoplotER(cfg,gavg2plot);
title(sprintf('%s: %i ms',sens2plot,timewin(2)*1000))
set(gca,'fontsize', 12)
colormap(bluewhitered)