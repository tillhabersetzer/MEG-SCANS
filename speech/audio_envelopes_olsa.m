%--------------------------------------------------------------------------
% Till Habersetzer, 15.07.2025
% Communication Acoustics, CvO University Oldenburg
% till.habersetzer@uol.de 
%
% Description:
%   Pre-processes OLSA sentence stimuli by computing their auditory
%   envelopes. This script applies the same resampling and bandpass
%   filtering steps used for the neural data to ensure consistency.
%
%   The final output is an interactive figure that allows for a
%   sentence-by-sentence visual comparison of the processed audio
%   and its corresponding envelope.
%--------------------------------------------------------------------------

close all
clearvars
clc 

%% Import main settings 
%--------------------------------------------------------------------------
settings_speech

% Addpath for additional functions
addpath(fullfile(settings.path2project,'analysis','helper_functions'))

% Select which final fs to show: 
fs_down = 1000;
% fs_down = 64;

fs_neuro = settings.decoding.fs_neuro;
fs_audio = settings.fs_audio;
% Epoch length for olsa trials
prestim  = settings.decoding.olsa.prestim;
poststim = settings.decoding.olsa.poststim;

% Zeropadding before filtering
padding = settings.decoding.olsa.padding;

% Filter settings
filtertype = settings.decoding.filtertype; % windowed sinc type I linear phase FIR filter
bpfreq     = settings.decoding.bpfreq; % meg + audio

%% Import data
%--------------------------------------------------------------------------
% Import olsa envelopes
bids_dir        = settings.path2bids;
stim_dir        = fullfile(bids_dir,'stimuli','olsa','sentences'); 
contents        = dir(stim_dir);
fnames_audio    = {contents.name};
fnames_audio    = fnames_audio(endsWith(fnames_audio, ".wav", "IgnoreCase", true));
n_sentences     = length(fnames_audio);
padding_samples = round(padding*fs_neuro);

olsa_orig      = cell(1,n_sentences);
olsa_envelopes = cell(1,n_sentences);

% Load original data
for stn_idx = 1:n_sentences % loop over sentences

    % Compute olsa envelope
    %----------------------------------------------------------------------
    fname_audio = fnames_audio{stn_idx};
    [audio, fs] = audioread(fullfile(stim_dir,fname_audio));
    if ~isequal(fs,fs_audio)
        error('Unexpected sampling frequency (%i)!',fs)
    end
    % Audio is stereo signal; focus on left channel which was also used in
    % the original experiment for both left and right ear.
    audio = audio(:,1);

    % Add prestim/poststim interval
    audio = [zeros(round(prestim*fs_audio),1);audio;zeros(round(poststim*fs_audio),1)];

    % Compute auditory envelope
    cfg      = [];
    cfg.type = 'auditory_envelope';
    cfg.fs   = fs_audio;
    envelope = cal_envelope(cfg, audio);

    % figure
    % hold on
    % plot(audio)
    % plot(envelope)
    % hold off

    % Resample entire audio to the intermediate rate (fs_neuro)
    % The 'resample' function includes an anti-aliasing filter.
    envelope_resampled = resample(envelope, fs_neuro, fs_audio);
    % Downsample as well
    audio_resampled = resample(audio, fs_neuro, fs_audio);

    % figure
    % hold on
    % plot(audio_resampled)
    % plot(envelope_resampled)
    % hold off

    % Apply the same bandpass filter used on the neuro-data
    % This requires the FieldTrip toolbox.
    % Tha audio data will be padded cause it is rather short the filtere
    % frequences very low and the filter sharp
    envelope_padded   = [zeros(1,padding_samples),envelope_resampled,zeros(1,padding_samples)];
    plotfiltresp      = 'no';
    envelope_filtered = ft_preproc_bandpassfilter(envelope_padded, fs_neuro, bpfreq, [], filtertype, [], [], [], [], [], plotfiltresp, []);
    envelope_filtered = envelope_filtered(padding_samples+1:end-padding_samples);
    
    if fs_down == 64
        audio_resampled   = resample(audio_resampled,fs_down,fs_neuro);
        envelope_filtered = resample(envelope_filtered, fs_down, fs_neuro);
    end

    % Store data in cell
    olsa_orig{stn_idx}      = audio_resampled;
    olsa_envelopes{stn_idx} = envelope_filtered;
    fprintf('Sentence %i/%i (%s) processed.\n',stn_idx,n_sentences,fname_audio)
    clear envelope envelope_resampled envelope_padded envelope_filtered audio audio_resampled
    
end
clear olsa_data

%% Visualize sentences
%--------------------------------------------------------------------------

% --- Figure and UI Initialization ---

% Create a handles structure to store all data and UI handles
handles               = struct();
handles.n_files       = n_sentences;
handles.fs            = fs_down;
handles.audio_data    = olsa_orig;
handles.envelope_data = olsa_envelopes;
handles.fnames_audio  = fnames_audio; % Store filenames
handles.current_index = 1;

% Create the main figure window
handles.fig = figure('Name', 'Interactive Sentence Plotter', 'NumberTitle', 'off', 'Color', 'w');
set(handles.fig, 'DoubleBuffer', 'on'); % Reduce flickering

% Create the subplots
handles.ax1 = subplot(3, 1, 1);
handles.ax2 = subplot(3, 1, 2);
handles.ax3 = subplot(3, 1, 3);

% Initialize plot handles so we can update them later
handles.h_raw = plot(handles.ax1, NaN, NaN, 'b');
title(handles.ax1, 'Raw Audio Data');
ylabel(handles.ax1, 'Amplitude');
grid(handles.ax1, 'on');

handles.h_env = plot(handles.ax2, NaN, NaN, 'r');
title(handles.ax2, 'Estimated Envelope');
ylabel(handles.ax2, 'Envelope');
grid(handles.ax2, 'on');

hold(handles.ax3, 'on');
handles.h_comb_raw = plot(handles.ax3, NaN, NaN, 'Color', 'b');
handles.h_comb_env = plot(handles.ax3, NaN, NaN, 'r', 'LineWidth', 1.5);
hold(handles.ax3, 'off');
title(handles.ax3, 'Combined View: Audio and Envelope');
xlabel(handles.ax3, 'Time / s');
ylabel(handles.ax3, 'Normalized Amplitude');
legend(handles.ax3, {'Audio', 'Envelope'}, 'Location', 'northeast');
grid(handles.ax3, 'on');

% Link axes so they zoom and pan together
linkaxes([handles.ax1, handles.ax2, handles.ax3], 'x');

% --- Interactive UI Controls ---

% "Previous" Button
uicontrol('Style', 'pushbutton', 'String', '<< Previous', ...
          'Units', 'normalized', 'Position', [0.1, 0.01, 0.1, 0.05], ...
          'Callback', @previous_callback);
      
% "Next" Button
uicontrol('Style', 'pushbutton', 'String', 'Next >>', ...
          'Units', 'normalized', 'Position', [0.8, 0.01, 0.1, 0.05], ...
          'Callback', @next_callback);
      
% Slider for navigating through files
handles.h_slider = uicontrol('Style', 'slider', ...
                     'Units', 'normalized', 'Position', [0.22, 0.01, 0.56, 0.05], ...
                     'Min', 1, 'Max', handles.n_files, 'Value', handles.current_index, ...
                     'SliderStep', [1/(handles.n_files-1), 1/(handles.n_files-1)], ... % Step by 1
                     'Callback', @slider_callback);

% Save the handles structure to the figure's guidata
guidata(handles.fig, handles);

% Initial plot of the first sentence
update_plot(handles.fig);

% --- Local Callback Functions ---
% These functions are defined at the end of the script. They use guidata
% to retrieve and modify the shared handles structure.

function slider_callback(source, ~)
    % Get the figure handle and the handles structure
    fig = gcbf; 
    handles = guidata(fig);
    
    % Get integer value from slider and update the index
    handles.current_index = round(get(source, 'Value'));
    
    % Save the updated handles structure
    guidata(fig, handles);
    
    % Update the plot
    update_plot(fig);
end

function previous_callback(~, ~)
    fig = gcbf;
    handles = guidata(fig);
    
    % Decrease index, but not below 1
    if handles.current_index > 1
        handles.current_index = handles.current_index - 1;
    end
    
    guidata(fig, handles);
    update_plot(fig);
end

function next_callback(~, ~)
    fig = gcbf;
    handles = guidata(fig);
    
    % Increase index, but not beyond the number of files
    if handles.current_index < handles.n_files
        handles.current_index = handles.current_index + 1;
    end
    
    guidata(fig, handles);
    update_plot(fig);
end

function update_plot(fig)
    % This is the main function that redraws the plots for a given index
    handles = guidata(fig);
    
    % Get the data for the current sentence
    audio    = handles.audio_data{handles.current_index};
    envelope = handles.envelope_data{handles.current_index};

    % Get the current filename
    current_fname = handles.fnames_audio{handles.current_index};
    
    % Create the time vector for this specific sentence
    t = (0:length(audio)-1) / handles.fs;
    
    % --- Update Plot Data (more efficient than replotting) ---
    set(handles.h_raw, 'XData', t, 'YData', audio);
    set(handles.h_env, 'XData', t, 'YData', envelope);
    set(handles.h_comb_raw, 'XData', t, 'YData', audio / max(abs(audio)));
    set(handles.h_comb_env, 'XData', t, 'YData', envelope / max(abs(envelope)));
    
    % --- Update Titles ---
    % Update the main figure title with the current filename
    sgtitle(handles.fig, sprintf('Displaying: %s', current_fname), 'Interpreter', 'none');
    
    % Update subplot titles to show the current sentence number
    title(handles.ax1, 'Raw Audio Data');
    title(handles.ax2, 'Estimated Envelope');
    title(handles.ax3, 'Combined View');
    
    % Adjust x-axis limits to fit the current sentence length
    xlim(handles.ax1, [0, t(end)]);
    
    % Update the slider's position to stay in sync
    set(handles.h_slider, 'Value', handles.current_index);
    
    % Redraw the figure
    drawnow;
end