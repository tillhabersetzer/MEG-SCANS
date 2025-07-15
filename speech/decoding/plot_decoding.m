%--------------------------------------------------------------------------
% Till Habersetzer, 15.07.2025
% Communication Acoustics, CvO University Oldenburg
% till.habersetzer@uol.de 
%
% Description:
%   Analyzes and visualizes neural decoding results. This script evaluates
%   model performance on held-out test data from the training phase (audio-
%   book stimulus) and tests its generalization on the OLSA sentence task.
%
%   Key Steps:
%   1.  Loads pre-computed decoding results for all subjects.
%   2.  Plots model training performance, including cross-validation curves
%       and test accuracy (sorted vs. shuffled).
%   3.  Visualizes generalization performance by plotting decoding accuracy
%       against OLSA sentence SNR and intelligibility.
%   4.  Performs statistical tests to correlate decoding accuracy with OLSA
%       stimulus properties.
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
subjects   = 1:24;
n_subj     = length(subjects);
lambda_exp = settings.decoding.decoder.lambda_exp;
n_lambdas  = length(lambda_exp);
lambdas    = 10.^(lambda_exp );
n_snrs     = 6; % 6 different snrs per subject

% Take absolute value of correlation
% apply_abs = true;
apply_abs = false;

% accuracy metric
corr_metric = settings.decoding.decoder.corr_metric; 

%% Import and extract data
%--------------------------------------------------------------------------
subjectnames = cell(1,n_subj);

% Arrays - Audiobook
%-------------------
lambda_vals = zeros(1,n_subj); % lambda regularization values used during analysis
cv          = cell(1,n_subj);

% All trials and division into training and test trials
n_trials       = zeros(1,n_subj);
n_trials_train = zeros(1,n_subj);
n_trials_test  = zeros(1,n_subj);

model_acc          = zeros(1,n_subj); % correlation value of sorted epochs > shuffled epochs
model_r            = zeros(2,n_subj); % 2: mean correlation value of sorted epochs and shuffled epochs
model_r_std        = zeros(2,n_subj);

% Arrays - Olsa
%--------------
pred_acc          = zeros(120,n_subj); % correlation values
pred_acc_shuffled = zeros(120,100,n_subj); % correlation values, 100 permutations
pred_err          = zeros(120,n_subj); 
snr_data          = zeros(120,n_subj); % sentence snrs
intelligibilities = zeros(120,n_subj); % sentence intelligibilities
snrs              = zeros(n_snrs,n_subj); % unique snrs

for sub_idx = 1:n_subj
    subject               = sprintf('sub-%02d',subjects(sub_idx));
    subjectnames{sub_idx} = subject;

    data = importdata(fullfile(settings.path2derivatives,subject,'speech',sprintf('%s_decoding.mat',subject)));   
    % Extract prediction scores for olsa
    pred_acc(:,sub_idx) = data.stats_olsa.r; 
    pred_err(:,sub_idx) = data.stats_olsa.err;

    % shuffled data
    pred_acc_shuffled(:,:,sub_idx) = data.stats_olsa_shuffled.r;

    % Extract event information for olsa
    snr_data(:,sub_idx)          = data.event_description.snrs;
    intelligibilities(:,sub_idx) = data.event_description.intells;
    % Check unique SNRs
    snrs(:,sub_idx) = sort(unique(data.event_description.snrs),'ascend');

    % Extract information of trained model
    n_trials(sub_idx)       = data.eval_model.n_trials;
    n_trials_train(sub_idx) = data.eval_model.n_trials_train;
    n_trials_test(sub_idx)  = data.eval_model.n_trials_test;
    lambda_vals(sub_idx)    = data.eval_model.lambda;
    cv{sub_idx}             = data.eval_model.cv;

    % Compute statistics for trained model
    %-------------------------------------
    stats_sorted   = data.eval_model.stats_sorted;
    stats_shuffled = data.eval_model.stats_shuffled;

    % Metric: Absolute value of correlation
    if apply_abs
        model_acc(1,sub_idx)   = sum(abs(stats_sorted.r) > abs(stats_shuffled.r))*100/length(stats_sorted.r);
        model_r(:,sub_idx)     = [mean(abs(stats_sorted.r));mean(abs(stats_shuffled.r))]; % sorted / shuffled condition
        model_r_std(:,sub_idx) = [std(abs(stats_sorted.r),0,1);std(abs(stats_shuffled.r),0,1)] ./sqrt(data.eval_model.n_trials_test); % sorted / shuffled condition
    % Metric: Simple Correlation
    else
        model_acc(1,sub_idx)   = sum(stats_sorted.r > stats_shuffled.r)*100/length(stats_sorted.r);
        model_r(:,sub_idx)     = [mean(stats_sorted.r);mean(stats_shuffled.r)];
        model_r_std(:,sub_idx) = [std(stats_sorted.r,0,1);std(stats_shuffled.r,0,1)] ./sqrt(data.eval_model.n_trials_test);
    end

    clear data
    fprintf('%s loaded.\n',subject)
    clear stats_sorted stats_shuffled

end

%% Plot training and tes trials across subjects
%--------------------------------------------------------------------------

figure; 
bar(1:n_subj, [n_trials',n_trials_train',n_trials_test']); 
hold off
ax                    = gca; 
ax.XTick              = 1:n_subj; 
ax.XTickLabel         = subjectnames;
ax.XTickLabelRotation = 45; 
title('Number of trials (training, test overall)');
ylabel('Number of Trials');
xlabel('Subject');
ylim([0,max(n_trials)+5])
grid on; 
legend({'all trials','training','test'},'Location','southoutside', 'Orientation', 'horizontal')

%% Evaluate goodness of individual trained model 
%--------------------------------------------------------------------------

% decoding accuracies, correlation values and regularization parameters
figure;

% Plot correlation values and accuracy
%--------------------------------------------------------------------------
ymin        = min(model_r,[],'all') - 0.1;
ymax        = max(model_r,[],'all') + 0.1;
x_positions = 1:n_subj;  % X-axis locations for subjects
bar_width   = 0.3; % common bar width
green_color = [0.47 0.67 0.19]; % Define the green color once

subplot(2,3,1:3)

% Left Y-Axis: Correlation values grouped bars (sorted, shuffled) with errors 
%--------------------------------------------------------------------------
yyaxis left;
% Plot the grouped bars, slightly offset to the left
b = bar(x_positions - 0.15, model_r', bar_width); % Narrower bars
hold on;

% Add the error bars for the grouped plot
for i = 1:size(model_r, 1) % Loop over groups
    x_coords = b(i).XEndPoints;
    y_coords = b(i).YEndPoints;
    errors   = model_r_std(i, :);
    errorbar(x_coords, y_coords, errors, 'k', 'linestyle', 'none');
end
ylabel(sprintf('Correlation value (%s)',corr_metric));
ylim([ymin, ymax]);

% Right Y-Axis: Accuray metric (correctly decoded trials: sorted > shuffle) 
%--------------------------------------------------------------------------
yyaxis right;
% Plot the single accuracy bar, slightly offset to the right
b_acc = bar(x_positions + 0.15, model_acc, bar_width/2, 'FaceColor', green_color); % Green
ylabel('Acc / %');
ylim([0, 105]); 

hold off;

% Final Formatting
%-----------------
b(1).FaceColor = 'b'; 
b(2).FaceColor = 'r'; 
ax                    = gca;
ax.YAxis(1).Color     = 'k'; % Set left axis color to black
ax.YAxis(2).Color     = green_color; % Set right axis color to green
ax.XTick              = x_positions;
ax.XTickLabel         = subjectnames;
ax.XTickLabelRotation = 45;
xlabel('Subject');
title('Model Correlation Values and Accuracies across subjects');
grid on;

% Create a legend for all three data series
legend([b(1), b(2), b_acc], {'Correlation (sorted)', 'Correlation (shuffled)', 'Accuracy'}, 'Location', 'southoutside','Orientation', 'horizontal');

% Regularization parameters
%--------------------------------------------------------------------------
cmap = distinguishable_colors(n_subj);

subplot(2,3,4)
hold on;
for sub_idx = 1:n_subj
    cv_subj = cv{sub_idx};
    n_fold  = size(cv_subj.r,1); % number of folds for cross-validation
    errorbar(1:n_lambdas,mean(cv_subj.r,1),std(cv_subj.r,1)/sqrt(n_fold-1),'linewidth',1,'Color',cmap(sub_idx,:))
end
set(gca,'xtick',1:n_lambdas,'xticklabel',lambda_exp);
xlim([0,numel(lambdas)+1]);
grid on
title('Cross-Validation Accuracy');
xlabel('Regularization (1\times10^\lambda)');
ylabel('Correlation');

subplot(2,3,5)
hold on;
for sub_idx = 1:n_subj
    cv_subj = cv{sub_idx};
    n_fold  = size(cv_subj.err,1); % number of folds for cross-validation
    errorbar(1:n_lambdas,mean(cv_subj.err,1),std(cv_subj.err,1)/sqrt(n_fold-1),'linewidth',1,'Color',cmap(sub_idx,:))
end
set(gca,'xtick',1:n_lambdas,'xticklabel',lambda_exp);
xlim([0,numel(lambdas)+1]);
grid on
title('Cross-Validation Error')
xlabel('Regularization (1\times10^\lambda)')
ylabel('MSE')
set(gca, 'YScale', 'log')

subplot(2,3,6)
plot(1:n_subj,log10(lambda_vals),'x')
set(gca,'xtick',1:n_subj,'xticklabel',subjectnames,'xticklabelrotation',45);
grid on
title('Regularization value across subjects')
xlabel('subjects')
ylabel('Regularization (1\times10^\lambda)')
ylim([lambda_exp(1),lambda_exp(end)]);

% Shared Legend Creation 
%-----------------------
num_cols     = 12;
h_legend_ax = axes('Position', [0, 0, 1, 1], 'Visible', 'off');
% Plot invisible "dummy" lines on these new axes to generate handles
hold(h_legend_ax, 'on');
dummy_handles = gobjects(n_subj, 1);
for i = 1:n_subj
    dummy_handles(i) = plot(h_legend_ax, NaN, NaN, 'LineWidth', 2, 'Color', cmap(i,:));
end
hold(h_legend_ax, 'off');

% Create the legend from the dummy handles and place it
legend(dummy_handles, subjectnames, 'Location', 'southoutside', 'NumColumns', num_cols);

sgtitle('Evaluation of trained model (cross-validation results, accuracies, correlation values)')

%% Evaluation of olsa test data
%--------------------------------------------------------------------------

% Computation of olsa scores
%--------------------------------------------------------------------------

% subjects x snrs
olsa_corr              = zeros(n_snrs,n_subj); 
olsa_corr_std          = zeros(n_snrs,n_subj); 
olsa_corr_shuffled     = zeros(n_snrs,n_subj); 
olsa_corr_shuffled_std = zeros(n_snrs,n_subj); 
n_permutations         = size(pred_acc_shuffled,2);

for sub_idx = 1:n_subj % Loop over subjects
    for snr_idx = 1:n_snrs % Loop over snrs

        snr2check = snrs(snr_idx, sub_idx);
        idx       = (snr_data(:,sub_idx)==snr2check);
        if ~isequal(sum(idx),20)
            error('Unexpected number of olsa trials (%i)!',sum(idx))
        end
        corr_vals          = pred_acc(idx,sub_idx);
        corr_vals_shuffled = pred_acc_shuffled(idx,:,sub_idx);
        
        if apply_abs
            olsa_corr(snr_idx,sub_idx)     = mean(abs(corr_vals));
            olsa_corr_std(snr_idx,sub_idx) = std(abs(corr_vals))./sqrt(sum(idx));
            % shuffled
            olsa_corr_shuffled(snr_idx,sub_idx)     = mean(abs(corr_vals_shuffled),"all");
            olsa_corr_shuffled_std(snr_idx,sub_idx) = std(abs(corr_vals_shuffled),[],"all")./sqrt(sum(idx)*n_permutations);
        else
            olsa_corr(snr_idx,sub_idx)     = mean(corr_vals);
            olsa_corr_std(snr_idx,sub_idx) = std(corr_vals)./sqrt(sum(idx));
            % shuffled
            olsa_corr_shuffled(snr_idx,sub_idx)     = mean(corr_vals_shuffled,"all");
            olsa_corr_shuffled_std(snr_idx,sub_idx) = std(corr_vals_shuffled,[],"all")./sqrt(sum(idx)*n_permutations);
        end  
        
        clear snr2check corr_vals
    end
end

% Plot SNR-curve for each subject
%--------------------------------------------------------------------------
ymin = min(olsa_corr - olsa_corr_std,[],'all') - 0.05;
ymax = max(olsa_corr + olsa_corr_std,[],'all') + 0.05;
xmin = min(snr_data, [], 'all') - 1;
xmax = max(snr_data, [], 'all') + 1;

figure('Color', 'w')
for sub_idx = 1:n_subj

    subplot(4,6,sub_idx);
    hold on;  
    % plot(snrs(:,sub_idx), olsa_corr(:,sub_idx),'LineWidth',1,'LineStyle','--','Color', 'k','Marker', 'x','MarkerSize',14); 
    errorbar(snrs(:,sub_idx), olsa_corr(:,sub_idx), olsa_corr_std(:,sub_idx), 'LineWidth', 1, "LineStyle", "none", 'Color','k', 'Marker', 'x');
    % highlight first and last point
    plot([snrs(1,sub_idx),snrs(end,sub_idx)], [olsa_corr(1,sub_idx), olsa_corr(end,sub_idx)], 'LineStyle', 'none', 'LineWidth', 2, 'Color', 'r','Marker', 'x','MarkerSize',10,'HandleVisibility','off'); 
    % Add shuffled condition
    errorbar(snrs(:,sub_idx), olsa_corr_shuffled(:,sub_idx), olsa_corr_shuffled_std(:,sub_idx), 'LineWidth', 1, "LineStyle", "none", 'Color',[0.7, 0.7, 0.7, 0.4], 'Marker', 'x');
    title(subjectnames{sub_idx});
    grid on;
    grid minor;
    box on; 
    % Set the same axis limits for every subplot
    axis([xmin, xmax, ymin, ymax]);

    % xlabel
    if ismember(sub_idx,19:24)
       xlabel('SNR / dB')
    end

    % ylabel
    if ismember(sub_idx,[1,7,13,19])
        if apply_abs
            ylabel(sprintf("|r_{%s}|",corr_metric))
        else
            ylabel(sprintf("r_{%s}",corr_metric))
        end
    end

    % Legend
    if sub_idx == 1
        legend({'sorted','shuffled'},'location','northwest')
    end

end % loop over subjects
sgtitle('Correlation values against SNR')

% Plot data for all subjects
%--------------------------------------------------------------------------
intelli_vec = sort(unique(intelligibilities),'ascend');
% Compute grand average across subjects (given fixed intelligibilities
% across subjects)
gavg_olsa_corr              = mean(olsa_corr,2); % mean
gavg_olsa_corr_sem          = std(olsa_corr,0,2)./sqrt(n_subj); % standard error
gavg_olsa_corr_shuffled     = mean(olsa_corr_shuffled,2); % mean
gavg_olsa_corr_shuffled_sem = std(olsa_corr_shuffled,0,2)./sqrt(n_subj); % standard error

axis_font_size   = 18; 
legend_font_size = 16;
title_font_size  = 18; 
  
figure('Color', 'w')

subplot(1,3,1)
hold on;
plot(snrs(:),olsa_corr_shuffled(:),'LineWidth',1,'Color', [1, 0, 0, 0.4],'Marker', 'x','MarkerSize', 10, 'LineStyle', 'none', 'MarkerSize', 10)
plot(snrs(:),olsa_corr(:),'LineWidth',1,'Color', 'k','Marker', 'x','MarkerSize', 10, 'LineStyle', 'none', 'MarkerSize', 10)
hold off;
xlabel('SNR / dB')
if apply_abs
    ylabel(sprintf("|r_{%s}|",corr_metric))
else
    ylabel(sprintf("r_{%s}",corr_metric))
end
title('Correlation against SNR', 'FontSize', title_font_size)
axis square;
grid on;
grid minor;
box on; 
axis([xmin, xmax, ymin, ymax]);
legend({'shuffled','sorted'},'location','northwest')

subplot(1,3,2)
hold on;
for sub_idx = 1:n_subj

    plot([intelli_vec(1),intelli_vec(end)], [olsa_corr(1,sub_idx), olsa_corr(end,sub_idx)], 'LineStyle', '-', 'LineWidth', 2, 'Color', cmap(sub_idx,:),'Marker', 'x', 'MarkerSize', 15); 

end
hold off;
xlabel('Intelligibility / %')
if apply_abs
    ylabel(sprintf("|r_{%s}|",corr_metric))
else
    ylabel(sprintf("r_{%s}",corr_metric))
end
set(gca,'xtick',[intelli_vec(1),intelli_vec(end)],'xticklabel',[intelli_vec(1),intelli_vec(end)]);
title('Correlation against Intelligibility', 'FontSize', title_font_size)
axis square;
grid on;
grid minor;
box on; 
axis([intelli_vec(1)-5, intelli_vec(end)+5, ymin, ymax]);

subplot(1,3,3)
hold on;
errorbar(intelli_vec, gavg_olsa_corr_shuffled, gavg_olsa_corr_shuffled_sem, 'LineWidth', 1, 'LineStyle', 'none', 'Color',[1, 0, 0, 0.4], 'Marker', 'x', 'MarkerSize', 15);
errorbar(intelli_vec, gavg_olsa_corr, gavg_olsa_corr_sem, 'LineWidth', 1, 'LineStyle', 'none', 'Color','k', 'Marker', 'x', 'MarkerSize', 15);
hold off;
xlabel('Intelligibility / %')
if apply_abs
    ylabel(sprintf("|r_{%s}|",corr_metric))
else
    ylabel(sprintf("r_{%s}",corr_metric))
end
title('Correlation against Intelligibility', 'FontSize', title_font_size)
axis square;
grid on;
grid minor;
box on; 
xlim([intelli_vec(1)-5, intelli_vec(end)+5])
% axis([0, 100, ymin, ymax]);

sgtitle('Overview plot')

% Shared Legend Creation 
%-----------------------
num_cols    = 12;
h_legend_ax = axes('Position', [0, 0, 1, 1], 'Visible', 'off');
hold(h_legend_ax, 'on');
dummy_handles = gobjects(n_subj, 1);
for i = 1:n_subj
    dummy_handles(i) = plot(h_legend_ax, NaN, NaN, 'LineWidth', 2, 'Color', cmap(i,:));
end
hold(h_legend_ax, 'off');
legend(dummy_handles, subjectnames, 'Location', 'southoutside', 'NumColumns', num_cols, 'FontSize', legend_font_size);

% Formatting
%-----------
all_axes = findall(gcf, 'type', 'axes');
set(all_axes, 'FontSize', axis_font_size);

%% Statistics
%--------------------------------------------------------------------------

% Compute correlation between snr and correlation values
% (subplot(1,3,1)
[r_pearson1,p_pearson1]   = corr(snrs(:),olsa_corr(:),'Type','Pearson'); % Pearson's Linear Correlation Coefficient
[r_spearman1,p_spearman1] = corr(snrs(:),olsa_corr(:),'Type','Spearman'); % Spearman's rho 

fprintf('\nSNR - Correlation Coefficient: Pearson correlation %.3f (p=%.5f).\n', r_pearson1, p_pearson1);
fprintf('\nSNR - Correlation Coefficient: Spearman correlation %.3f (p=%.5f).\n', r_spearman1, p_spearman1);
    
% Compute correlation between intelligibilty and correlation values
% (subplot(1,3,3)
[r_pearson2,p_pearson2]   = corr(repmat(intelli_vec,n_subj,1),olsa_corr(:),'Type','Pearson'); % Pearson's Linear Correlation Coefficient
[r_spearman2,p_spearman2] = corr(repmat(intelli_vec,n_subj,1),olsa_corr(:),'Type','Spearman'); % Spearman's rho 

fprintf('\nIntelligibility - Correlation Coefficient: Pearson correlation %.3f (p=%.5f).\n', r_pearson2, p_pearson2);
fprintf('\nIntelligibility - Correlation Coefficient: Spearman correlation %.3f (p=%.5f).\n', r_spearman2, p_spearman2);
    
% Compute paired-samples t-test between intellibility 20% and 80%
% (subplot(1,3,3)
[httest,pttest] = ttest(olsa_corr(1,:),olsa_corr(end,:));
fprintf('\n Pairswise difference between Intelligibility %i%% and %i%% has mean equal to zero: %s (p=%.5f).\n', intelli_vec(1), intelli_vec(end), string(logical(~httest)), pttest);

% Compute paired-samples t-test between sorted and shuffled audiobooks 
% (subplot(1,3,3)
[httest,pttest] = ttest(model_r(1,:),model_r(2,:));
fprintf('\n Pairswise difference between sorted and shuffled audiobook has mean equal to zero: %s (p=%.5f).\n', string(logical(~httest)), pttest);

%% Final plot for paper
%--------------------------------------------------------------------------
axis_font_size   = 16; 
legend_font_size = 16;
title_font_size  = 18; 
text_font_size   = 16;
  
ymin        = min(olsa_corr,[],'all') - 0.01;
ymax        = max(olsa_corr,[],'all') + 0.01;
xmin        = min(snr_data, [], 'all') - 1;
xmax        = max(snr_data, [], 'all') + 1;
ymin1       = min(model_r,[],'all') - 0.1;
ymax1       = max(model_r,[],'all') + 0.1;
x_positions = 1:n_subj;  
bar_width   = 0.9;

figure('Color', 'w');
subplot(2,2,1:2)
hold on;
for sub_idx = 1:n_subj
    %     b1 = bar(x_positions, model_r(1,:), bar_width, 'FaceColor', [0 0.4470 0.7410],'FaceAlpha',0.9); 
    %     errorbar(x_positions, model_r(1,:),model_r_std(1,:), 'k', 'linestyle', 'none');
    %     b2 = bar(x_positions, model_r(2,:), bar_width, 'FaceColor', [0.8500 0.3250 0.0980],'FaceAlpha',0.9); 
    %     errorbar(x_positions, model_r(2,:),model_r_std(2,:), 'k', 'linestyle', 'none');
    b1 = bar(x_positions(sub_idx), model_r(1,sub_idx), bar_width, 'FaceColor', cmap(sub_idx,:),'FaceAlpha',0.8); 
    errorbar(x_positions(sub_idx), model_r(1,sub_idx), model_r_std(1,sub_idx), 'k', 'linestyle', 'none','CapSize',20,'LineWidth',1.5);
    b2 = bar(x_positions(sub_idx), model_r(2,sub_idx), bar_width, 'FaceColor', [0.7,0.7,0.7],'FaceAlpha',1); % always grey
    errorbar(x_positions(sub_idx), model_r(2,sub_idx), model_r_std(2,sub_idx), 'k', 'linestyle', 'none','CapSize',20,'LineWidth',1.5);

    h1 = hatchfill2(b1, 'single', 'HatchAngle', 45, 'HatchColor', 'k', 'HatchLineWidth', 1);
end
hold off;
if apply_abs
    ylabel(sprintf("|r_{%s}|",corr_metric))
else
    ylabel(sprintf("r_{%s}",corr_metric))
end
set(gca,'xtick',1:n_subj,'xticklabel',subjectnames,'XTickLabelRotation',45);
title('Audiobook: Correlation Values', 'FontSize', title_font_size)
grid on;
grid minor;
box on; 
ylim([ymin1, ymax1]);
%legend([b1,b2], {'sorted','shuffled'}, 'Location', 'northeast');
legend([h1,b2], {'sorted','shuffled'}, 'Location', 'northeast');

subplot(2,2,3)
hold on;
for sub_idx = 1:n_subj
    plot(snrs(:,sub_idx),olsa_corr(:,sub_idx),'LineWidth',1,'Color', cmap(sub_idx,:),'Marker', 'x','MarkerSize', 10, 'LineStyle', 'none')
end
xlabel('SNR / dB')
if apply_abs
    ylabel(sprintf("|r_{%s}|",corr_metric))
else
    ylabel(sprintf("r_{%s}",corr_metric))
end
title('OLSA: Correlation Values against SNR', 'FontSize', title_font_size)
grid on;
grid minor;
box on; 
axis([xmin, xmax, ymin, ymax]);

% Add p-val into plot  
%--------------------
ptext   = sprintf('$r \\approx %.2f \\left[p<10^{%d}\\right]$', r_pearson1, ceil(log10(p_pearson1)));
ax1     = gca; % Get current axes
xLimits = get(ax1, 'XLim');
yLimits = get(ax1, 'YLim');
text(xLimits(1) + .05*range(xLimits), ... % x-position with 5% margin
     yLimits(2) - .05*range(yLimits), ... % y-position with 5% margin
     ptext, ...
     'VerticalAlignment', 'top', ...
     'HorizontalAlignment', 'left', ...
     'BackgroundColor', 'white', ...
     'EdgeColor', 'black', ...
     'Interpreter', 'latex', ...
     'FontSize',text_font_size);

subplot(2,2,4)
hold on;
for sub_idx = 1:n_subj
    plot(intelli_vec,olsa_corr(:,sub_idx),'LineWidth',1, 'Color', [cmap(sub_idx,:),0.8],'Marker','x','MarkerSize',10,'LineStyle', 'none');
end
errorbar(intelli_vec, gavg_olsa_corr, gavg_olsa_corr_sem, 'LineWidth', 2, 'LineStyle', 'none', 'Color','k', 'Marker', 'x', 'MarkerSize', 20);
hold off;
xlabel('Intelligibility / %')
if apply_abs
    ylabel(sprintf("|r_{%s}|",corr_metric))
else
    ylabel(sprintf("r_{%s}",corr_metric))
end
title('OLSA: Correlation Values against Intelligibility', 'FontSize', title_font_size)
grid on;
grid minor;
box on; 
axis([intelli_vec(1)-5, intelli_vec(end)+5, ymin, ymax])

% Add p-val into plot  
%--------------------
ptext   = sprintf('$r \\approx %.2f \\left[p<10^{%d}\\right]$', r_pearson2, ceil(log10(p_pearson2)));
ax1     = gca; % Get current axes
xLimits = get(ax1, 'XLim');
yLimits = get(ax1, 'YLim');
text(xLimits(1) + .05*range(xLimits), ... % x-position with 5% margin
     yLimits(2) - .05*range(yLimits), ... % y-position with 5% margin
     ptext, ...
     'VerticalAlignment', 'top', ...
     'HorizontalAlignment', 'left', ...
     'BackgroundColor', 'white', ...
     'EdgeColor', 'black', ...
     'Interpreter', 'latex', ...
     'FontSize',text_font_size);

% Formatting
%-----------
all_axes = findall(gcf, 'type', 'axes');
set(all_axes, 'FontSize', axis_font_size);