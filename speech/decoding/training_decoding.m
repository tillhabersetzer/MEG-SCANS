function training_decoding(subject,settings)
%--------------------------------------------------------------------------
% Till Habersetzer, 27.01.2025
% Communication Acoustics, CvO University Oldenburg
% till.habersetzer@uol.de 
%
% Description:
%   Performs subject-level training and evaluation of a neural decoding
%   model (backward model) using the mTRF Toolbox.
%
%   Key Steps:
%   -   Trains a temporal response function (TRF) model to reconstruct a
%       speech envelope from MEG data, using a continuous audiobook dataset.
%   -   Finds the optimal regularization parameter (lambda) for the model
%       using leave-one-out cross-validation on the training data.
%   -   Evaluates the model on a held-out portion of the audiobook data,
%       comparing performance against a shuffled-trial baseline.
%   -   Tests the model's generalization by using the audiobook-trained
%       decoder to predict envelopes for a separate OLSA sentence task.
%   -   Saves the trained model, its evaluation statistics, and the OLSA
%       predictions to a single .mat file for the subject.
%
% Inputs:
%   subject (char/string): The subject identifier (e.g., 'sub-01').
%   settings (struct):     Settings structure with all necessary paths and 
%                          parameters.
%--------------------------------------------------------------------------

%% Script settings 
%--------------------------------------------------------------------------
% Apply zscoring of all trials
apply_zscoring = settings.decoding.zscore;

% Sanity check plots
plot_results = false;

%% Import data
%--------------------------------------------------------------------------

% Load audiobooks
data                   = importdata(fullfile(settings.path2project,'derivatives',subject,'speech',sprintf('%s_preprocessed_audiobooks_decoding.mat',subject)));   
epochs_audio_audiobook = data.epochs_audio;
epochs_neuro_audiobook = data.epochs_neuro;
n_trials               = length(epochs_neuro_audiobook.trial);
clear data
% Load olsa
data              = importdata(fullfile(settings.path2project,'derivatives',subject,'speech',sprintf('%s_preprocessed_olsa_decoding.mat',subject)));   
epochs_audio_olsa = data.epochs_audio;
epochs_neuro_olsa = data.epochs_neuro;
event_description = data.event_description;
clear data
fprintf("\nData from %s is loaded.\n",subject)

% Sanity check
n_trials_olsa = length(epochs_neuro_olsa.trial);
if ~isequal(n_trials_olsa,120)
    error('%s: Unexpected number of olsa trials (%i)!',subject, n_trials_olsa)
end

%% Apply zscoring
%--------------------------------------------------------------------------
% Z-scoring is applied using a single mean and standard deviation
% calculated from all data points across all trials combined. This ensures
% a consistent normalization across the entire dataset.

% For neuro data, this global approach preserves the relative amplitude
% differences (i.e., the topography) between channels, which would be
% lost if each channel were z-scored independently.
% For the audio data, the relative power across features is preserved.

if apply_zscoring   
    % Audiobooks
    %----------------------------------------------------------------------
    epochs_audio_audiobook_concat = horzcat(epochs_audio_audiobook{:});
    epochs_neuro_audiobook_concat = horzcat(epochs_neuro_audiobook.trial{:});
    
    mean_val_audio = mean(epochs_audio_audiobook_concat,2);
    std_val_audio  = std(epochs_audio_audiobook_concat,0,2);
    mean_val_neuro = mean(epochs_neuro_audiobook_concat,'all');
    std_val_neuro  = std(epochs_neuro_audiobook_concat,0,'all');
    
    for trl_idx = 1:n_trials
        % audio
        epochs_audio_audiobook{trl_idx} = (epochs_audio_audiobook{trl_idx}-mean_val_audio)./std_val_audio;
        % neuro
        epochs_neuro_audiobook.trial{trl_idx} = (epochs_neuro_audiobook.trial{trl_idx}-mean_val_neuro)./std_val_neuro;
    end
    clear epochs_audio_audiobook_concat epochs_neuro_audiobook_concat mean_val_audio std_val_audio mean_val_neuro std_val_neuro
    
    % Olsa
    %----------------------------------------------------------------------
    epochs_audio_olsa_concat = horzcat(epochs_audio_olsa{:});
    epochs_neuro_olsa_concat = horzcat(epochs_neuro_olsa.trial{:});
    
    mean_val_audio = mean(epochs_audio_olsa_concat,2);
    std_val_audio  = std(epochs_audio_olsa_concat,0,2);
    mean_val_neuro = mean(epochs_neuro_olsa_concat,'all');
    std_val_neuro  = std(epochs_neuro_olsa_concat,0,'all');

    for trl_idx = 1:n_trials_olsa
        % audio
        epochs_audio_olsa{trl_idx} = (epochs_audio_olsa{trl_idx}-mean_val_audio)./std_val_audio;
        % neuro
        epochs_neuro_olsa.trial{trl_idx} = (epochs_neuro_olsa.trial{trl_idx}-mean_val_neuro)./std_val_neuro;
    end
    clear epochs_audio_olsa_concat epochs_neuro_olsa_concat mean_val_audio std_val_audio mean_val_neuro std_val_neuro
    fprintf("\n %s: Zcoring finished.\n",subject)

end


%% Train decoder on audio books
%--------------------------------------------------------------------------

n_trials_train = ceil(0.8*n_trials); % number of training trials
n_trials_test  = n_trials - n_trials_train; % number of test trials


% Cross-validation for optimization
%----------------------------------
% To optimize the decoders ability to predict stimulus features from new 
% EEG data, we tune the regularization parameter using an efficient 
% leave-one-out cross-validation (CV) procedure.

% take random trials for training and test
rng("shuffle")
p = randperm(n_trials);

% lambda_exp = -6:2:6;
lambda_exp = settings.decoding.decoder.lambda_exp;
lambdas    = 10.^(lambda_exp );

% Run fast cross-validation
cv = mTRFcrossval(epochs_audio_audiobook(p(1:n_trials_train)), ... % stimulus
                  epochs_neuro_audiobook.trial(p(1:n_trials_train)), ... % response
                  epochs_neuro_audiobook.fsample, ... % sampling rate (Hz)
                  settings.decoding.decoder.drct, ... % direction
                  settings.decoding.decoder.tmin, ... % minimum time lag (ms)
                  settings.decoding.decoder.tmax, ... % maximum time lag (ms)
                  lambdas, ... % regularization values
                  'dim',settings.decoding.decoder.dim,... % work along the rows, observations in columns
                  'zeropad',settings.decoding.decoder.zeropad, ... % zero-pad the outer rows of the design matrix or delete them
                  'fast',settings.decoding.decoder.fast, ... % use the fast cross-validation method (requires more memory)
                  'corr', settings.decoding.decoder.corr_metric, ... % Spearman's rank correlation coefficient
                  'type', settings.decoding.decoder.type); % use all lags simultaneously to fit a multi-lag model

if plot_results
    % Plot CV accuracy 
    %-----------------
    nlambda = length(lambdas);
    nfold   = n_trials_train;
    figure
    subplot(1,2,1), errorbar(1:numel(lambdas),mean(cv.r),std(cv.r)/sqrt(nfold-1),'linewidth',2)
    set(gca,'xtick',1:nlambda,'xticklabel',lambda_exp), xlim([0,numel(lambdas)+1]), axis square, grid on
    title('CV Accuracy'), xlabel('Regularization (1\times10^\lambda)'), ylabel('Correlation')
    % Plot CV error
    %--------------
    subplot(1,2,2), errorbar(1:numel(lambdas),mean(cv.err),std(cv.err)/sqrt(nfold-1),'linewidth',2)
    set(gca,'xtick',1:nlambda,'xticklabel',lambda_exp), xlim([0,numel(lambdas)+1]), axis square, grid on
    title('CV Error'), xlabel('Regularization (1\times10^\lambda)'), ylabel('MSE')
end

% Get optimal  regularization hyperparameter
% (nfold-by-nlambda-by-yvar) == (20x7lambdas)
[rmax,idx] = max(mean(cv.r,1));
lambda     = lambdas(idx);

% Train model on training data
%-----------------------------
model = mTRFtrain(epochs_audio_audiobook(p(1:n_trials_train)), ...
                  epochs_neuro_audiobook.trial(p(1:n_trials_train)), ...
                  epochs_neuro_audiobook.fsample, ...
                  settings.decoding.decoder.drct, ... 
                  settings.decoding.decoder.tmin, ...
                  settings.decoding.decoder.tmax, ...
                  lambda, ...
                  'dim',settings.decoding.decoder.dim, ...
                  'zeropad',settings.decoding.decoder.zeropad, ...
                  'type', settings.decoding.decoder.type);

%% Test trained decoder model
%--------------------------------------------------------------------------

% Sorted data
%------------
epochs_audio_sorted = epochs_audio_audiobook(p(n_trials_train+1:end));
[pred,stats_sorted] = mTRFpredict(epochs_audio_sorted, ...
                                  epochs_neuro_audiobook.trial(p(n_trials_train+1:end)), ...
                                  model, ...
                                  'zeropad', settings.decoding.decoder.zeropad, ...
                                  'dim', settings.decoding.decoder.dim, ...
                                  'corr', settings.decoding.decoder.corr_metric);

% Shuffled data
%--------------
idx_shuffle = randperm(n_trials_test); % for random mapping of audiodata

% check if no elements are identical, so the difference should never be 0
while ~all(idx_shuffle-(1:n_trials_test)) % check for nonzero elements
    idx_shuffle = randperm(n_trials_test);
end
trials                = p(n_trials_train+1:end);
epochs_audio_shuffled = epochs_audio_audiobook(trials(idx_shuffle));

[~,stats_shuffled] = mTRFpredict(epochs_audio_shuffled, ...
                                 epochs_neuro_audiobook.trial(p(n_trials_train+1:end)), ...
                                 model, ...
                                 'zeropad', settings.decoding.decoder.zeropad, ...
                                 'dim', settings.decoding.decoder.dim, ...
                                 'corr', settings.decoding.decoder.corr_metric);

% Plot reconstruction 
%--------------------
if plot_results
    fs       = epochs_neuro_audiobook.fsample;
    xlim_max = settings.decoding.trialdur;
    figure
    subplot(2,1,1)
    plot((1:length(epochs_audio_sorted{1}))/fs,epochs_audio_sorted{1},'linewidth',2), hold on % original sorted audio
    plot((1:length(epochs_audio_shuffled{1}))/fs,epochs_audio_shuffled{1},'linewidth',2),  % original shuffled audio
    plot((1:length(pred{1}))/fs,pred{1},'linewidth',2), hold off, xlim([0,xlim_max]), grid on % predicted audio
    title('Reconstruction'), xlabel('Time (s)'), ylabel('Amplitude (a.u.)'), legend('Orig sorted','Orig shuffled','Pred sorted')
    
    % Plot test accuracy
    subplot(2,1,2), bar(1,rmax), hold on, bar(2,mean(stats_sorted.r)), bar(3,mean(stats_shuffled.r)), hold off
    set(gca,'xtick',1:2,'xticklabel',{'Val.','Test sorted','Test shuffled'}), grid on
    title('Model Performance'), xlabel('Dataset'), ylabel('Correlation')
end
% clear epochs_audio_sorted epochs_audio_shuffled

% Store data
%-----------
eval_model                = struct();
eval_model.cv             = cv; % cross-validation
eval_model.lambda         = lambda; % regularization value
eval_model.n_trials       = n_trials;
eval_model.n_trials_train = n_trials_train;
eval_model.n_trials_test  = n_trials_test;
eval_model.model          = model;
eval_model.stats_sorted   = stats_sorted;
eval_model.stats_shuffled = stats_shuffled;

%% Prediction on olsa sentences
%--------------------------------------------------------------------------

% Check if channel order is identical in olsa and audiobook meg recordings
%--------------------------------------------------------------------------
[chanidx1,chanidx2] = match_str(epochs_neuro_audiobook.label,epochs_neuro_olsa.label);
if ~isequal(chanidx1,chanidx2)
    error('Channel labels have different order!')
end
clear chanidx1 chanidx2

% Compute predictions
%--------------------
[sentences_pred,stats_olsa] = mTRFpredict(epochs_audio_olsa, ...
                                          epochs_neuro_olsa.trial, ...
                                          model,...
                                          'zeropad', settings.decoding.decoder.zeropad, ...
                                          'dim', settings.decoding.decoder.dim, ...
                                          'corr', settings.decoding.decoder.corr_metric);   

%% Save results
%--------------
dir2save = fullfile(settings.path2derivatives,subject,'speech');
if ~exist(dir2save,'dir')
    mkdir(dir2save)
end
fname = sprintf('%s_decoding.mat',subject);

results                   = strcut();
results.eval_model        = eval_model; % evaluation of trained model
results.stats_olsa        = stats_olsa;
results.sentences_pred    = sentences_pred;
results.sentences_orig    = epochs_audio_olsa;
results.event_description = event_description;

save(fullfile(dir2save,fname),'results','-v7.3'); 
fprintf("\n%s from %s saved.\n",fname,subject)

end % end of function