function [envelope] = cal_envelope(audiodata,fs,b,a,type)
%--------------------------------------------------------------------------
% Till Habersetzer, 26.06.2025
% Communication Acoustics, CvO University Oldenburg
% till.habersetzer@uol.de 
% 
% CAL_ENVELOPE Extracts different types of envelopes from audio data.
%
%   [envelope] = CAL_ENVELOPE(audiodata, fs, b, a, type) computes either
%   the instantaneous amplitude envelope (log-scaled) or an onset-specific
%   envelope (half-wave rectified derivative) after Hilbert transform and
%   zero-phase filtering.
%
% Input Arguments:
%   audiodata (vector): Input audio time-series.
%   fs (numeric): Sampling frequency in Hz.
%   b (vector): Numerator filter coefficients.
%   a (vector): Denominator filter coefficients.
%   type (char): 'envelope' for log-scaled amplitude, or 'onset_envelope' for derivative.
%
% Output Arguments:
%   envelope (vector): The calculated envelope time-series.
%--------------------------------------------------------------------------

% Calculate instantaneous amplitude using Hilbert transform
envelope_inst = abs(hilbert(audiodata));

% Apply zero-phase digital filtering to the instantaneous amplitude
envelope_filtered = filtfilt(b, a, envelope_inst);

switch type
    case 'onset_envelope'
        % Petersen, Eline Borch, et al. "Neural tracking of attended versus 
        % ignored speech is differentially affected by hearing loss." 
        % Journal of neurophysiology 117.1 (2017): 18-27.
        
        % Calculate the derivative 
        envelope_diff = diff(envelope_filtered) * fs;

        % Half-wave rectification: set negative values to zero
        envelope_diff(envelope_diff < 0) = 0;

        % Assign the result to the output variable
        envelope = envelope_diff;
        
    case 'envelope'
        % Ensure no negative values before log transformation
        envelope_filtered(envelope_filtered < 0) = 0;

        % Apply log transformation for dB-like scale (add epsilon to avoid log(0))
        envelope_log = 20 * log10(envelope_filtered + 1e-6);

        % Assign the result to the output variable
        envelope = envelope_log;

    otherwise
        % Throw an error for unsupported types
        error('CAL_ENVELOPE:UnsupportedType', 'Specified envelope type "%s" is not supported! Use "envelope" or "onset_envelope".', type);
end

end