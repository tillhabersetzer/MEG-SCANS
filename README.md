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
* Training of a decoder based on the audiobooks to infer speech intelligibility encoded in the MEG recordings (see: https://doi.org/10.1007/s10162-018-0654-z).

## Auditory Evoked Fields (AEFs)

The analysis of Auditory Evoked Fields (AEFs) is distinguished into sensor-level and source-level analyses. Settings for the AEF analysis are stored in `settings_chirps.m`. More detailed descriptions of each script's functionality can be found within the scripts themselves.

### Sensor Level Analysis 

| Script Name | Description                                                                                                                                                                                                                               
| --- | --- |
| `compute_aefs.m` | Computes auditory evoked fields for each subject and the grand average across all subjects. It processes a list of subjects by defining trials around auditory chirp stimuli, applying a full preprocessing pipeline with artifact rejection and filtering, and then calculating an average evoked response for each individual. After processing all subjects, the script computes and saves a final grand average. |
| `plot_aefs_subjectlevel.m` | Loads a single subject's pre-computed Auditory Evoked Field (AEF) to generate a comprehensive set of visualizations. It produces multi-channel sensor layouts, single-channel waveforms, and 2D topographic maps. Each plot type is created for magnetometers, planar gradiometers, and a combined planar gradiometer representation. |
| `plot_aefs_grouplevel.m` | Aggregates pre-computed Auditory Evoked Field (AEF) data from all subjects to generate final summary figures. It produces overview plots of the grand average, including multi-channel layouts and topographies. This script generates Figure XX in the publication ([DOI to data descriptor paper]). For selected channels, it creates detailed visualizations that overlay every subject's response on the grand average and also plot the grand average with its standard error to illustrate variability. |

### Source Level Analysis

| Script Name | Description |    
| --- | --- |
| `compute_headmodel_sourcemodel.m` | Automates the creation and validation of head and source models from anatomical MRIs for a batch of subjects for use in MEG source analysis. For each individual, the script sequentially co-registers the MRI with the MEG sensors, segments the brain volume, and computes subject-specific head and source models. The co-registration information is already provided by the BIDS dataset. For quality control, it generates a visualization of the alignment for each subject and then compiles all plots into a single HTML report for efficient visual inspection. |
| `check_coregistration.m` | Loads the subject's pre-computed files—including the head model, source model, segmented MRI, and sensor information—and generates a series of plots to verify their alignment. |
| `compute_dipolefit.m` | Performs dipole model fitting on individual subjects' MEG data to estimate the location and time course of neural sources. It uses a multi-step process, starting with a symmetric grid search to get an initial location, followed by a non-linear refinement to find the precise dipole positions. The dipole fit is performed separately on magnetometers and gradiometers. |                                                                                                                                                 
| `plot_dipolefit_subjectlevel.m` | Visualizes and validates a single subject's dipole fitting results. It plots the final dipole locations on the subject's MRI and identifies their corresponding anatomical labels using an atlas. The script also plots the reconstructed source time courses. |                                                                                                                                                                                                                                                        
| `plot_dipolefit_grouplevel.m` | This script generates Figure XX in the publication ([DOI to data descriptor paper]). It analyzes group dipole fitting results by first performing automated quality control to exclude subjects with poor fits. From the remaining valid data, it calculates the group-average dipole locations and source time courses. It then visualizes these results, creating a 3D plot of all individual and mean dipole locations on a template brain, and plots of the grand-average source activity over time with standard error. |                                                                                                                          
