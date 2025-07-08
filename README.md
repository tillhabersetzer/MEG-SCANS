# MEG-SCANS: Code for Technical Validation

This repository contains the code used for the technical validation and generation of plots for the MEG-SCANS (Stories, Chirps, And Noisy Sentences) dataset publication.

## Dataset and Publication

The MEG-SCANS (Stories, Chirps, And Noisy Sentences) dataset provides raw magnetoencephalography (MEG) recordings from 24 German-speaking participants. Each participant listened to approximately one hour of stimuli, including two audiobooks, sentences from the Oldenburger Matrix Sentence Test (OLSA) for Speech Reception Threshold (SRT) assessment, and short up-chirps used to assess MEG signal quality. The dataset comprises MEG data, corresponding audio material (audiobooks, OLSA envelopes, and chirp stimuli), and behavioral audiogram results from hearing screenings. Organized according to the Brain Imaging Data Structure (BIDS), this resource offers a robust benchmark for large-scale encoding/decoding analyses of temporally-resolved brain responses to speech. Comprehensive Matlab and Python code are included to replicate key data validations, ensuring transparency and reproducibility.

* **Dataset:** [DOI to OpenNeuro dataset]
* **Data Descriptor Paper:** [DOI to data descriptor paper]

## Technical Validation Overview

The technical validation covers the following aspects:

* Computation of Auditory Evoked Fields (AEFs) from chirp stimuli and a dipole fit of their N100m component.
* A cross-correlation analysis between audiobook envelopes and MEG data.
* A decoding analysis where models are trained on audiobooks and tested on OLSA sentences.

MATLAB R2023a and the FieldTrip toolbox (https://www.fieldtriptoolbox.org/) were used for the analysis.

## Auditory Evoked Fields (AEFs)

The AEF analysis is divided into sensor-level and source-level analyses, with settings configured in `settings_chirps.m`. More detailed descriptions of each script's functionality can be found within the scripts themselves.

### Sensor Level Analysis 

| Script Name | Description |                                                                                                                                                                                             
| :--- | :--- |
| `compute_aefs.m` | Processes chirp stimuli to compute subject-level AEFs. It then calculates and saves the grand average across all subjects. |
| `plot_aefs_subjectlevel.m` | Generates a comprehensive set of visualizations for a single subject's AEF data, including sensor layouts, single-channel visualizations, and topographies. |
| `plot_aefs_grouplevel.m` | Creates summary figures for the group-level AEF analysis. It plots the grand-average and overlays individual subject data for comparison. It generates Figure XX in the publication ([DOI to data descriptor paper]). |

### Source Level Analysis

| Script Name | Description |                                                                                                                                                                                               
| :--- | :--- |
| `compute_headmodel_sourcemodel.m` | Creates subject-specific head and source models from anatomical MRIs. It also generates an HTML report for quality control. An excerpt from this report is shown below. The co-registration information is already provided by the BIDS dataset. |
| `check_coregistration.m` | Loads pre-computed model files to generate plots that verify the coregistration alignment between sensors, head model, source model and mri. |
| `compute_dipolefit.m` | Performs dipole model fitting using a two dipole model on a subject's AEF data (N100m peak) to estimate the location of the underlying neural sources. |                                             | `plot_dipolefit_subjectlevel.m` | Visualizes a single subject's dipole fit results, including dipole locations on the MRI and the reconstructed source time courses. |                                                       | `plot_dipolefit_grouplevel.m` | Analyzes group dipole data by excluding poor fits and computing the group-average location and source time course. It then visualizes these group-level results. It generates Figure XX in the publication ([DOI to data descriptor paper]). |

## Cross-correlation analysis
This section describes the cross-correlation analysis performed for the technical validation of the speech material. This analysis examines the relationship between audiobook onset envelopes and magnetoencephalography (MEG) recordings, drawing inspiration from Petersen et al. (2016) [[doi:10.1152/jn.00527.2016]](https://doi.org/10.1152/jn.00527.2016). The specific settings for this cross-correlation analysis are configured in `settings_speech.m`. The analysis workflow is divided into the following scripts:

| Script Name | Description |                                                                                                                                                                                               
| :--- | :--- | 
| `analysis_pipeline_crosscorr.m` | This master script orchestrates the "complete" MEG-audio cross-correlation analysis pipeline (without plotting). It preprocesses the shared audio stimuli, then processes each subject in parallel to prepare their MEG data and compute subject-level cross-correlations. The pipeline concludes by computing and saving the group-level grand average from all subject results. The entire process is logged to a text file. |
| `preprocessing_audio.m` | This script prepares audio stimuli for a cross-correlation analysis. For each raw audiobook file, it computes the onset envelope, segments it into fixed-length epochs, and downsamples them, saving the final data as .mat files. |
| `preprocessing_crosscorr.m` | This script prepares a single subject's data for a cross-correlation analysis. It processes continuous MEG data by filtering and epoching it, then loads corresponding audio envelope data. Finally, it rigorously synchronizes the two datasets by rejecting invalid trials and truncating all remaining pairs to a uniform length before saving the result. |
| `computation_crosscorr.m` | This script performs a subject-level temporal cross-correlation analysis between MEG data and an audio envelope. For a given subject, it loads preprocessed data and calculates the trial-by-trial correlation for two conditions: one with the correctly aligned audio and a second, shuffled-trial control condition for statistical comparison. It then saves the trial-averaged results for both conditions. |
| `computation_gavg_crosscorr.m` | This script computes a group-level grand average from individual subject data. It loads the pre-computed, trial-averaged cross-correlation results for each subject, then uses FieldTrip to calculate the grand average for both the actual and the shuffled control conditions, saving the final result to a single file. |
| `plot_crosscorr.m` | This script visualizes and statistically analyzes pre-computed MEG cross-correlation data. It loads group-level results and performs a cluster-based permutation test to find significant differences between the actual data and a shuffled-trial null condition. It then generates a comprehensive set of plots, including topographical maps of the results and detailed time-series plots for selected channels. This script is responsible for creating Figure XX in the publication ([DOI to data descriptor paper]). |

## Decoding analysis
This section describes the decoding-based analysis used for the technical validation of the speech material, for both audiobook and olsa recording. The analysis aimed to reproduce the results of Vanthornhout et al. (2018) [[doi:10.1007/s10162-018-0654-z]](https://doi.org/10.1007/s10162-018-0654-z) by training a backward model (decoder) on continuous audiobook data and testing its ability to reconstruct the speech envelopes of unseen OLSA sentences.
The primary objective was to reconstruct the psychometric function of speech intelligibility using an objective, neural-based metric. The temporal response function (TRF) framework used for training and evaluating the decoder is described in detail by Crosse et al. (2021) [[doi:10.3389/fnins.2021.705621]](https://doi.org/10.3389/fnins.2021.705621) and O'Sullivan et al. (2014) [[doi:10.1093/cercor/bht355]](https://doi.org/10.1093/cercor/bht355).
The decoding analysis pipeline closely follows the structure of the cross-correlation pipeline. The analysis workflow is divided into the following scripts:

| Script Name | Description |                                                                                                                                                                                               
| :--- | :--- | 
| `analysis_pipeline_decoding.m` | This master script orchestrates the complete neural decoding analysis pipeline. It begins by preprocessing all audio stimuli (both OLSA sentences and audiobooks). It then uses a parallel loop to process each subject, which involves preparing their corresponding MEG data and then executing the main script to train and evaluate a subject-specific decoding model. The entire workflow and total execution time are logged to a text file. |
| `preprocessing_olsa.m` | This script prepares OLSA sentence stimuli for a neural decoding analysis. It individually processes each raw .wav file by computing its auditory envelope, applying padding, and then performing a multi-step filtering and resampling process to align with neural data parameters. All processed envelopes are saved together in a single .mat file. |
| `preprocessing_audiobooks.m` | This script prepares audiobook stimuli for a neural decoding analysis. For each raw .wav file, it computes the auditory envelope, then applies a multi-step process of filtering, resampling, and segmenting the envelope into fixed-length epochs. The final data is saved into separate .mat files for the pilot subject and main subject group. |
| `preprocessing_audiobooks_decoding.m` | This script prepares a single subject's audiobook data for a decoding analysis. It processes the continuous MEG recording by filtering, epoching, and downsampling it. The script then loads the corresponding pre-processed audio envelopes and rigorously synchronizes the two datasets by rejecting invalid trials and truncating all remaining pairs to a uniform length before saving the final data.  |
| `preprocessing_olsa_decoding.m` | This script prepares OLSA sentence stimuli for a neural decoding analysis. It individually processes each raw .wav file by computing its auditory envelope, applying padding, and then performing a multi-step filtering and resampling process to align with neural data parameters. All processed envelopes are saved together in a single .mat file. |
| `training_decoding.m` | This script trains and evaluates a subject-specific neural decoding model using the mTRF Toolbox. It first trains a Temporal Response Function (TRF) model to reconstruct speech envelopes from MEG data using a continuous audiobook dataset. It then evaluates the model's performance on a held-out portion of the audiobook data and tests its ability to generalize by using it to reconstruct speech envelopes from a separate OLSA sentence task. |
| `plot_decoding.m` | This script is responsible for creating Figure XX in the publication ([DOI to data descriptor paper]). |

## helper functions
The `audio_envelopes.m` script in the `speech` directory provides a tool for detailed visualization of the envelopes used in both the cross-correlation and decoding analyses.
| Script Name | Description |                                                                                                                                                                                               
| :--- | :--- | 
| `audio_envelopes.m` | This script is an interactive tool for visualizing speech envelopes. It processes a selected audiobook file to compute a specific type of envelope (e.g., auditory, onset) and then generates a multi-panel plot with a slider that allows for easy scrolling and comparison of the raw audio waveform and its derived envelope. |

An example of the coregistration report is shown below. This report is created for each subject to allow for visual inspection of the alignment between the sensor and anatomical data.
![Coregistration example](./images/coregistration_report_example.png)





other functions
* maxfilter
* freesurfer, check trigger not shared
* + maybe onclude picture of coregistration report

| `audio_envelope.m` | This script calculates the speech envelope for various audiobook stimuli, primarily as a check. (Note: The envelopes used for the cross-correlation analysis are actually computed in `preprocessing_crosscorr.m`.) It generates an interactive figure displaying the raw audio waveform and its calculated envelope in three synchronized plots. A slider control allows users to scroll through the audio over time for easy signal inspection and comparison. |



