%--------------------------------------------------------------------------
% Till Habersetzer, 24.07.2025
% Communication Acoustics, CvO University Oldenburg
% till.habersetzer@uol.de
%--------------------------------------------------------------------------

close all
clearvars
clc

%% Load data
%--------------------------------------------------------------------------
rawdata_path = fullfile('M:\meg_scans\bidsdata');

% Load audiograms
%----------------
subjects     = 1:24;
n_sub        = length(subjects);
audiograms   = zeros(12,2,n_sub); % freq x left/right x subjects
subjectnames = cell(1,n_sub);

for sidx = 1:n_sub

    sub_idx = subjects(sidx);

    subject              = sprintf('sub-%02d', sub_idx);
    subjectnames{sidx}   = subject;
    data                 = readtable(fullfile(rawdata_path,subject,'beh',sprintf('%s_task-audiogram_beh.tsv',subject)), 'FileType', 'text', 'Delimiter', '\t');
    
    frequency            = data.frequency;
    audiograms(:,:,sidx) = [data.hearingThreshold_left_,data.hearingThreshold_right_];
    clear data

end

%% Plot Audiograms
%--------------------------------------------------------------------------

audiogram_l = squeeze(audiograms(:,1,:))';
audiogram_r = squeeze(audiograms(:,2,:))';

figure('Name','Audiograms')

axs = cell(1,2);

axs{1} = subplot(2,1,1);
for i = 1:n_sub
    plot(frequency/1000,audiogram_r(i,:),'-o');
    hold on
end
% Add mean
errorbar(axs{1},frequency/1000,mean(audiogram_r,1),std(audiogram_r,1),'k-o','LineWidth',2)

legend([subjectnames,'mean + std'],'Location','eastoutside');
title('Right')
xlabel('f / lHz')
ylabel(' Level / dB HL')
grid on

axs{2} = subplot(2,1,2);
for i = 1:n_sub
    plot(frequency/1000,audiogram_l(i,:),'-o');
    hold on
end
% Add mean adn std
errorbar(axs{2},frequency/1000,mean(audiogram_l,1),std(audiogram_l,1),'k-o','LineWidth',2)

legend([subjectnames,'mean + std'],'Location','eastoutside');
title('Left')
xlabel('f / kHz')
ylabel(' Level / dB HL')
grid on


set([axs{:}],'xtick', frequency/1000, 'xticklabel',{'0.125','0.25','0.5','0.75','1','1.5','2','3','4','6','8','10'})
set([axs{:}],'XScale', 'log');
set([axs{:}], 'YDir','reverse','ylim',[-10,40])