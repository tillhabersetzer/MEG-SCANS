%--------------------------------------------------------------------------
% Till Habersetzer, 27.06.2025
% Communication Acoustics, CvO University Oldenburg
% till.habersetzer@uol.de 
%
% Description:
%   Calculates speech envelopes for audiobook stimuli and provides an
%   interactive, scrollable plot to visually compare the raw audio waveform
%   against its derived envelope. Supports multiple envelope types.
%--------------------------------------------------------------------------

close all
clearvars
clc 

%% Import main settings 
%--------------------------------------------------------------------------
settings_speech

% Addpath for additional functions
addpath(fullfile(settings.path2project,'analysis','helper_functions'))

%% Script settings 
%--------------------------------------------------------------------------

% Select audiobook & run
task = 'audiobook1'; % 'audiobook1' 'audiobook2'
run  = 1; % 1 2

% Select envelope type
type = 'auditory_envelope'; % type: 'envelope' 'onset_envelope' 'auditory_envelope'

if strcmpi(type, 'envelope') || strcmpi(type, 'onset_envelope')
    cfg_envelope              = [];
    cfg_envelope.type         = type;
    cfg_envelope.fs           = settings.fs_audio;
    cfg_envelope.lpfreq       = settings.crosscorr.audio.lpfreq;
    cfg_envelope.lpfiltord    = settings.crosscorr.audio.lpfiltord; 
    cfg_envelope.filtertype   = settings.crosscorr.filtertype;
    cfg_envelope.plotfiltresp = 'no';

    fs_down = settings.crosscorr.fs_down;
elseif strcmpi(type, 'auditory_envelope')
    cfg_envelope      = [];
    cfg_envelope.type = 'auditory_envelope';
    cfg_envelope.fs   = settings.fs_audio;

    fs_neuro          = 1000;
    fs_down           = settings.decoding.fs_down;
    filtertype        = settings.decoding.filtertype;
    plotfiltresp      = 'no';
    bpfreq            = settings.decoding.bpfreq;
else
    error('CAL_ENVELOPE: Unsupported Type (%s)!',type);
end

% Add fieldtrip
% addpath(settings.path2fieldtrip)
% ft_defaults

bids_dir = settings.path2bids;
stim_dir = fullfile(bids_dir,'stimuli','audiobooks');

%% Proprocess audiobook
%--------------------------------------------------------------------------

% Import raw audio
fname              = sprintf('task-%s_run-%s_stim.wav',task,sprintf('0%i',run));
[raw_audiodata,fs] = audioread(fullfile(stim_dir, fname)); 

% if strcmpi(type, 'auditory_envelope')
%     % only take first 2 minutes due to memory requirements
%     raw_audiodata = raw_audiodata(1:120*fs);
% end

envelope = cal_envelope(cfg_envelope, raw_audiodata);

if strcmpi(type, 'auditory_envelope')
    envelope      = resample(envelope, fs_neuro, settings.fs_audio);
    envelope      = ft_preproc_bandpassfilter(envelope, fs_neuro, bpfreq, [], filtertype, [], [], [], [], [], plotfiltresp, []);
    raw_audiodata = resample(raw_audiodata, fs_neuro, settings.fs_audio);
end

audiobook_label = fname(1:end-9);
fprintf('%s processed.\n', fname)
   
% Optionally: Downsample data before plotting
%--------------------------------------------

if strcmpi(type, 'auditory_envelope')
    raw_audiodata_down = resample(raw_audiodata,fs_down,fs_neuro);
    envelope_down      = resample(envelope,fs_down,fs_neuro);
else
    raw_audiodata_down = resample(raw_audiodata,fs_down,fs);
    envelope_down      = resample(envelope,fs_down,fs);
end
fprintf('%s resampled.\n', audiobook_label)

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
        envelope2plot      = envelope_down;
        fs2plot            = fs_down;
    case 'no'
        raw_audiodata2plot = raw_audiodata;
        envelope2plot      = envelope;

        if strcmpi(type, 'auditory_envelope')
            fs2plot = fs_neuro;
        else
            fs2plot = fs;
        end
end

% Define the width of the visible x-axis in seconds. This determines the
% duration of the audio segment shown at any one time.
time_win = 10; % in seconds

% Generate a time vector 't' in seconds for the x-axis of the plots.
t = (0:length(raw_audiodata2plot)-1) / fs2plot;

% Figure and Plot Initialization
%-------------------------------
figure;
% Enable double buffering to reduce flickering during slider updates.
set(gcf, 'DoubleBuffer', 'on');

% Subplot 1: Raw Audio Data 
subplot(3, 1, 1);
plot(t, raw_audiodata2plot);
title('Raw Audio Data');
ylabel('Amplitude');
grid on;

% Subplot 2: Speech Onset Envelope 
subplot(3, 1, 2);
% The envelope is one sample shorter (due to derivative), so we adjust the time vector
plot(t(1:length(envelope2plot)), envelope2plot);
title('Speech Envelope');
ylabel('Envelope');
grid on;

% Subplot 3: Combined Plot 
subplot(3, 1, 3);
plot(t, raw_audiodata2plot/max(abs(raw_audiodata2plot)), 'Color', [0.7 0.7 0.7]); % Plot audio in gray
hold on;
% The envelope might be one sample shorter.
plot(t(1:length(envelope2plot)), envelope2plot/max(abs(envelope2plot)), 'r', 'LineWidth', 1.5); % Plot envelope in red, adjust in height for plotting
hold off;
title('Combined View: Audio and Envelope');
xlabel('Time (s)');
ylabel('Amplitude / Envelope');
legend('Raw Audio', 'Envelope');
grid on;

sgtitle(sprintf('Envelope: %s | Analysis of %s', type, audiobook_label));

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