# MEG-SCANS: Code for Technical Validation

This repository contains the code used for the technical validation and generation of plots for the MEG-SCANS (Stories, Chirps, AND Noisy Sentences) dataset publication.

## Dataset and Publication

The MEG-SCANS (Stories, Chirps, And Noisy Sentences) dataset provides raw magnetoencephalography (MEG) recordings from 24 German-speaking participants. Each participant listened to approximately one hour of stimuli, including two audiobooks, sentences from the Oldenburger Matrix Sentence Test (OLSA) for Speech Reception Threshold (SRT) assessment, and short up-chirps used to assess MEG signal quality. The dataset comprises MEG data, corresponding audio material (audiobooks, OLSA envelopes, and chirp stimuli), and behavioral audiogram results from hearing screenings. Organized according to the Brain Imaging Data Structure (BIDS), this resource offers a robust benchmark for large-scale encoding/decoding analyses of temporally-resolved brain responses to speech. Comprehensive Matlab and Python code are included to replicate key data validations, ensuring transparency and reproducibility.

* **Dataset:** [DOI to OpenNeuro dataset]
* **Data Descriptor Paper:** [DOI to data descriptor paper]

## Technical Validation Overview

The technical validation covers the following aspects:

* Computation of auditory evoked fields (AEFs) based on chirp stimuli and a dipole fit of their N100m component.
* Computation of cross-correlation functions between presented audiobooks and MEG recordings, including a null distribution of permuted audio.
* Training of a decoder based on the audiobooks to infer speech intelligibility encoded in the MEG recordings.

MATLAB R2023a and the FieldTrip toolbox (https://www.fieldtriptoolbox.org/) were used for the analysis.

## Auditory Evoked Fields (AEFs)

The analysis of Auditory Evoked Fields (AEFs) is distinguished into sensor-level and source-level analyses. Settings for the AEF analysis are stored in `settings_chirps.m`. More detailed descriptions of each script's functionality can be found within the scripts themselves.

### Sensor Level Analysis 

| Script Name | Description |                                                                                                                                                                                             
| :--- | :--- |
| `compute_aefs.m` | Computes auditory evoked fields for each subject and the grand average across all subjects. It processes a list of subjects by defining trials around auditory chirp stimuli, applying a full preprocessing pipeline with artifact rejection and filtering, and then calculating an average evoked response for each individual. After processing all subjects, the script computes and saves a final grand average. |
| `plot_aefs_subjectlevel.m` | Loads a single subject's pre-computed Auditory Evoked Field (AEF) to generate a comprehensive set of visualizations. It produces multi-channel sensor layouts, single-channel waveforms, and 2D topographic maps. Each plot type is created for magnetometers, planar gradiometers, and a combined planar gradiometer representation. |
| `plot_aefs_grouplevel.m` | Aggregates pre-computed Auditory Evoked Field (AEF) data from all subjects to generate final summary figures. It produces overview plots of the grand average, including multi-channel layouts and topographies. This script generates Figure XX in the publication ([DOI to data descriptor paper]). For selected channels, it creates detailed visualizations that overlay every subject's response on the grand average and also plot the grand average with its standard error to illustrate variability. |

### Source Level Analysis

| Script Name | Description |                                                                                                                                                                                               
| :--- | :--- |
| `compute_headmodel_sourcemodel.m` | Automates the creation and validation of head and source models from anatomical MRIs for a batch of subjects for use in MEG source analysis. For each individual, the script sequentially co-registers the MRI with the MEG sensors, segments the brain volume, and computes subject-specific head and source models. The co-registration information is already provided by the BIDS dataset. For quality control, it generates a visualization of the alignment for each subject and then compiles all plots into a single HTML report for efficient visual inspection. |
| `check_coregistration.m` | Loads the subject's pre-computed files—including the head model, source model, segmented MRI, and sensor information—and generates a series of plots to verify their alignment. |
| `compute_dipolefit.m` | Performs dipole model fitting on individual subjects' MEG data to estimate the location and time course of neural sources. It uses a multi-step process, starting with a symmetric grid search to get an initial location, followed by a non-linear refinement to find the precise dipole positions. The dipole fit is performed separately on magnetometers and gradiometers. |                                                 | `plot_dipolefit_subjectlevel.m` | Visualizes and validates a single subject's dipole fitting results. It plots the final dipole locations on the subject's MRI and identifies their corresponding anatomical labels using an atlas. The script also plots the reconstructed source time courses. |                                                                                                                                                          | `plot_dipolefit_grouplevel.m` | This script generates Figure XX in the publication ([DOI to data descriptor paper]). It analyzes group dipole fitting results by first performing automated quality control to exclude subjects with poor fits. From the remaining valid data, it calculates the group-average dipole locations and source time courses. It then visualizes these results, creating a 3D plot of all individual and mean dipole locations on a template brain, and plots of the grand-average source activity over time with standard error. |

## Cross-correlation analysis
This section describes the cross-correlation analysis performed for the technical validation of the speech material. This analysis examines the relationship between audiobook onset envelopes and magnetoencephalography (MEG) recordings, drawing inspiration from Petersen et al. (2016) [[doi:10.1152/jn.00527.2016]](doi:10.1152/jn.00527.2016). The specific settings for this cross-correlation analysis are configured in `settings_speech.m`. The analysis workflow is divided into the following scripts:

| Script Name | Description |                                                                                                                                                                                               
| :--- | :--- | 
| `analysis_pipeline_crosscorr.m` | This master script orchestrates the "complete" MEG-audio cross-correlation analysis pipeline (without plotting). It preprocesses the shared audio stimuli, then processes each subject in parallel to prepare their MEG data and compute subject-level cross-correlations. The pipeline concludes by computing and saving the group-level grand average from all subject results. The entire process is logged to a text file. |
| `preprocessing_audio.m` | This script prepares audio stimuli for a cross-correlation analysis. For each raw audiobook file, it computes the onset envelope, segments it into fixed-length epochs, and downsamples them, saving the final data as .mat files. |
| `preprocessing_crosscorr.m` | This script prepares a single subject's data for a cross-correlation analysis. It processes continuous MEG data by filtering and epoching it, then loads corresponding audio envelope data. Finally, it rigorously synchronizes the two datasets by rejecting invalid trials and truncating all remaining pairs to a uniform length before saving the result. |
| `computation_crosscorr.m` | This script performs a subject-level temporal cross-correlation analysis between MEG data and an audio envelope. For a given subject, it loads preprocessed data and calculates the trial-by-trial correlation for two conditions: one with the correctly aligned audio and a second, shuffled-trial control condition for statistical comparison. It then saves the trial-averaged results for both conditions. |
| `computation_gavg_crosscorr.m` | This script computes a group-level grand average from individual subject data. It loads the pre-computed, trial-averaged cross-correlation results for each subject, then uses FieldTrip to calculate the grand average for both the actual and the shuffled control conditions, saving the final result to a single file. |
| `plot_crosscorr.m` | This script visualizes and statistically analyzes pre-computed MEG cross-correlation data. It loads group-level results and performs a cluster-based permutation test to find significant differences between the actual data and a shuffled-trial null condition. It then generates a comprehensive set of plots, including topographical maps of the results and detailed time-series plots for selected channels. This script is responsible for creating Figure XX in the publication ([DOI to data descriptor paper]). |

## Decoding analysis
This section describes the decoding-based analysis used for the technical validation of the speech material, for both audiobook and olsa recording. The analysis aimed to reproduce the results of Vanthornhout et al. (2018) [[doi:10.1007/s10162-018-0654-z]](doi:10.1007/s10162-018-0654-z) by training a backward model (decoder) on continuous audiobook data and testing its ability to reconstruct the speech envelopes of unseen OLSA sentences.
The primary objective was to reconstruct the psychometric function of speech intelligibility using an objective, neural-based metric. The temporal response function (TRF) framework used for training and evaluating the decoder is described in detail by Crosse et al. (2021) [[doi:10.3389/fnins.2021.705621]](doi:10.3389/fnins.2021.705621) and O'Sullivan et al. (2014) [[doi:10.1093/cercor/bht355]](doi:10.1093/cercor/bht355).






## helper functions

other functions
* maxfilter
* freesurfer, check trigger not shared
* + maybe onclude picture of coregistration report

| `audio_envelope.m` | This script calculates the speech envelope for various audiobook stimuli, primarily as a check. (Note: The envelopes used for the cross-correlation analysis are actually computed in `preprocessing_crosscorr.m`.) It generates an interactive figure displaying the raw audio waveform and its calculated envelope in three synchronized plots. A slider control allows users to scroll through the audio over time for easy signal inspection and comparison. |



