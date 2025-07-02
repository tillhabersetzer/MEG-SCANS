%--------------------------------------------------------------------------
% Till Habersetzer, 27.06.2025
% Communication Acoustics, CvO University Oldenburg
% till.habersetzer@uol.de 
%
% This script preprocesses audiobook stimuli by calculating their speech
% envelopes and provides a visualization tool to inspect the raw audio
% alongside its derived envelope with a scrollable interface.
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
% For test purposes of eeg reference, only select single session and
% condition

% Select seubjects
% subjects = 2:24; % all files available for sub-02-sub-23
subjects = 4;
n_subj   = length(subjects);

% Filter settings
bpfreq = settings.crosscorr.meg.bpfreq; % meg
lpfreq = settings.crosscorr.audio.lpfreq.audio; % audio

% Epoch length
trialdur = settings.crosscorr.trialdur;

% Downsampling frequency
fs_down = settings.crosscorr.fs_down;

% check whether maxfiltered data should be analyzed
use_maxfilter = settings.use_maxfilter;

% Latency Correction
apply_latency_correction = settings.apply_latency_correction;
audio_latency            = settings.audio_latency;

bids_dir = settings.path2bids;

% Add fieldtrip
% addpath(settings.path2fieldtrip)
% ft_defaults

% Addpath for additional functions
addpath(fullfile(settings.path2project,'analysis','helper_functions'))


%% Proprocess audiobooks
%--------------------------------------------------------------------------

% check if filter is okay
fsamp = 44100;
[b,a] = butter(3,lpfreq/(fsamp/2),'low');
% h     = fvtool(b,a);
% h.Fs  = fsamp;

stim_dir         = fullfile(bids_dir,'stimuli','audiobooks');
envelopes        = cell(1,4);
raw_audiodata    = cell(1,4);
audiobook_labels = cell(1,4);

f_idx = 1;
for task_idx = 1:2     
    task = sprintf('audiobook%i',task_idx);
    for run_idx = 1:2
        run = sprintf('0%i',run_idx);

        % Import raw audio
        fname                     = sprintf('task-%s_run-%s_stim.wav',task,run);
        [raw_audiodata{f_idx},fs] = audioread(fullfile(stim_dir, fname)); 

        envelopes{f_idx}        = cal_envelope(raw_audiodata{f_idx},fs,b,a,'onset_envelope');
        audiobook_labels{f_idx} = fname(1:end-9);
        fprintf('%s processed.\n', fname)
        f_idx                   = f_idx + 1;
    end
end

% Optionally: Downsample data before plotting
%--------------------------------------------

envelopes_down     = cell(size(envelopes));
raw_audiodata_down = cell(size(raw_audiodata));

for f_idx = 1:4
    raw_audiodata_down{f_idx} = resample(raw_audiodata{f_idx},fs_down,fs);
    envelopes_down{f_idx}     = resample(envelopes{f_idx},fs_down,fs);
    fprintf('%s resampled.\n', audiobook_labels{f_idx})
end


%% Have a look at speech onset envelope
%--------------------------------------------------------------------------
% This section visualizes raw audio data and its corresponding speech onset
% envelope. It presents three subplots: the raw audio, the envelope, and
% a combined view. A slider control is included to allow for horizontal
% scrolling through the data, which is particularly useful for inspecting
% audio files.

% Choose which audiodata to show
downsample = 'yes';
switch downsample
    case 'yes'
        raw_audiodata2plot = raw_audiodata_down;
        envelopes2plot     = envelopes_down;
        fs2plot            = fs_down;
    case 'no'
        raw_audiodata2plot = raw_audiodata;
        envelopes2plot     = envelopes;
        fs2plot            = fs;
end

% Select the audio file to process by choosing an index (e.g., 1, 2, 3, or 4).
f_idx = 2;

% Define the width of the visible x-axis in seconds. This determines the
% duration of the audio segment shown at any one time.
time_win = 10; % in seconds

% Generate a time vector 't' in seconds for the x-axis of the plots.
t = (0:length(raw_audiodata2plot{f_idx})-1) / fs2plot;

% Figure and Plot Initialization
%-------------------------------
figure;
% Enable double buffering to reduce flickering during slider updates.
set(gcf, 'DoubleBuffer', 'on');

% Subplot 1: Raw Audio Data 
subplot(3, 1, 1);
plot(t, raw_audiodata2plot{f_idx});
title('Raw Audio Data');
ylabel('Amplitude');
grid on;

% Subplot 2: Speech Onset Envelope 
subplot(3, 1, 2);
% The envelope is one sample shorter (due to derivative), so we adjust the time vector
plot(t(1:length(envelopes2plot{f_idx})), envelopes2plot{f_idx});
title('Speech Onset Envelope');
ylabel('Envelope');
grid on;

% Subplot 3: Combined Plot 
subplot(3, 1, 3);
plot(t, raw_audiodata2plot{f_idx}, 'Color', [0.7 0.7 0.7]); % Plot audio in gray
hold on;
% The envelope might be one sample shorter.
plot(t(1:length(envelopes2plot{f_idx})), envelopes2plot{f_idx}/max(envelopes2plot{f_idx}), 'r', 'LineWidth', 1.5); % Plot envelope in red, adjust in height for plotting
hold off;
title('Combined View: Audio and Envelope');
xlabel('Time (s)');
ylabel('Amplitude / Envelope');
legend('Raw Audio', 'Envelope');
grid on;

sgtitle(sprintf('Analysis of: %s', audiobook_labels{f_idx}));

% Link the x-axes of all subplots so they scroll and zoom together.
linkaxes(findobj(gcf, 'type', 'axes'), 'x');

% Add Scrollability
%------------------
% Set the initial x-axis limits for the viewing window.
xlim([0, time_win]);

% Create a slider for horizontal scrolling.
uicontrol('Style', 'slider', ...
          'Units', 'normalized', ...
          'Position', [0.1, 0.01, 0.8, 0.04], ... % [left, bottom, width, height]
          'Min', 0, ...
          'Max', t(end) - time_win, ...
          'Value', 0, ... % Start at the beginning
          'Callback', @(slider, ~) scrollPlots(slider, time_win));

%% Functions
%--------------------------------------------------------------------------
function scrollPlots(slider, windowSize)
    % Callback function for the slider.
    % This function is executed whenever the slider's value changes.
    %
    % Inputs:
    %   slider     - The handle to the slider uicontrol object.
    %   windowSize - The width of the viewing window in seconds.

    % Get the current value of the slider (the start of the new view).
    scrollValue = get(slider, 'Value');

    % Set the x-axis limits of the current axes ('gca') to the new position.
    % Because the axes are linked, all subplots will update simultaneously.
    xlim(gca, [scrollValue, scrollValue + windowSize]);
end