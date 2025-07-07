function preprocessing_olsa(settings)
%--------------------------------------------------------------------------
% Till Habersetzer, 05.07.2025
% Communication Acoustics, CvO University Oldenburg
% till.habersetzer@uol.de 
%
% Description:
%   Prepares individual OLSA sentence stimuli for a neural decoding analysis
%   by converting each raw sentence file into a processed auditory envelope.
%
%   Key Steps:
%   -   For each raw audio (.wav) sentence file, it computes the auditory
%       envelope and adds pre- and post-stimulus zero-padding.
%   -   Applies a multi-step filtering and resampling process to each sentence
%       envelope to align it with neural data processing parameters.
%   -   Saves a single structure containing all processed sentence envelopes
%       as separate fields to a final .mat file.
%
% Input:
%   settings (struct): Settings containing all necessary paths and
%                      processing parameters for the decoding analysis.
%--------------------------------------------------------------------------

%% Script settings 
%--------------------------------------------------------------------------

% Filter settings
filtertype = settings.decoding.filtertype; % windowed sinc type I linear phase FIR filter
bpfreq     = settings.decoding.bpfreq; % meg + audio

% Downsampling frequency
fs_down = settings.decoding.fs_down;

% Audio frequency
fs_audio = settings.fs_audio;

% MEG sampling frequency
fs_neuro = settings.decoding.fs_neuro;

% Epoch length for olsa trials
prestim  = settings.decoding.olsa.prestim;
poststim = settings.decoding.olsa.poststim;

% Zeropadding before filtering
padding = settings.decoding.olsa.padding;

bids_dir = settings.path2bids;
stim_dir = fullfile(bids_dir,'stimuli','olsa','sentences'); 

%% Process audio
%--------------------------------------------------------------------------
% Loop over all Olsa sentences
contents        = dir(stim_dir);
fnames_audio    = {contents.name};
fnames_audio    = fnames_audio(endsWith(fnames_audio, ".wav", "IgnoreCase", true));
n_sentences     = length(fnames_audio);
padding_samples = round(padding*fs_neuro);

audio_envelopes    = struct();
audio_envelopes.fs = fs_down;

for stn_idx = 1:n_sentences % loop over sentences

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

    % Resample entire audio to the intermediate rate (fs_neuro)
    % The 'resample' function includes an anti-aliasing filter.
    envelope_resampled = resample(envelope, fs_neuro, fs_audio);

    % Apply the same bandpass filter used on the neuro-data
    % This requires the FieldTrip toolbox.
    % Tha audio data will be padded cause it is rather short the filtere
    % frequences very low and the filter sharp
    envelope_padded   = [zeros(1,padding_samples),envelope_resampled,zeros(1,padding_samples)];
    plotfiltresp      = 'no';
    envelope_filtered = ft_preproc_bandpassfilter(envelope_padded, fs_neuro, bpfreq, [], filtertype, [], [], [], [], [], plotfiltresp, []);
    envelope_filtered = envelope_filtered(padding_samples+1:end-padding_samples);

    % figure
    % hold on
    % plot(envelope_resampled)
    % plot(envelope_filtered)
    % legend({'padded','filtered'})
    % -> with and without padding leads to same result

    % Resample the extracted epoch to the final target frequency (fs_down)
    [~, base_name, ~]                                  = fileparts(fname_audio);  
    audio_envelopes.(sprintf('envelope_%s',base_name)) = resample(envelope_filtered, fs_down, fs_neuro);
    fprintf('Sentence %i/%i (%s) processed.\n',stn_idx,n_sentences,fname_audio)
    clear envelope envelope_resampled envelope_padded envelope_filtered

end % sentences
    
% Save results
%--------------
dir2save = fullfile(settings.path2derivatives,'stimuli');
if ~exist(dir2save,'dir')
    mkdir(dir2save)
end
fname = 'preprocessed_olsa_envelopes_decoding.mat';

save(fullfile(dir2save,fname),'audio_envelopes','-v7.3'); 
fprintf("%s saved.\n",fname)

end % end of function