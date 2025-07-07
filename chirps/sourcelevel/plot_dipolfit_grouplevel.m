%--------------------------------------------------------------------------
% Till Habersetzer, 24.06.2025
% Communication Acoustics, CvO University Oldenburg
% till.habersetzer@uol.de 
% 
% Group Analysis and Visualization of Dipole Fits
%
% This script performs a group-level analysis of MEG dipole fitting 
% results obtained from a chirp stimulus experiment.
%
% Key steps:
% 1. Loads pre-computed dipole data for a specified range of subjects.
% 2. Performs quality control using the `check_dipolefitting` subfunction 
%    to exclude subjects with poor fits (e.g., unilateral or outlier dipoles).
% 3. Calculates the group-average dipole positions, moments, and time courses
%    for the left and right hemispheres, along with the standard error of the mean.
% 4. Generates visualizations:
%    - A 3D plot of individual and mean dipoles on a template MRI.
%    - Plots of the grand-average dipole moment time courses (both vector 
%      components and scalar magnitude) for each hemisphere.
%--------------------------------------------------------------------------

close all
clearvars
clc 

%% Import main settings 
%--------------------------------------------------------------------------
current_dir = pwd;
cd(fullfile('..'))
settings_chirps
cd(current_dir)

% Addpath for additional functions
addpath(fullfile(settings.path2project,'analysis','helper_functions'))

%% Script settings
%--------------------------------------------------------------------------
% Choose subject for plotting
subjects  = 2:24;
n_subj    = length(subjects);

chan2plot         = 'megplanar'; % 'megmag' 'megplanar'
time2plot         = [-100 400]; % for dipole timecourse
timewin_dipolefit = settings.timewindow_dipfit*1000; % in ms

%% Import data
%--------------------------------------------------------------------------
template_mri  = ft_read_mri(fullfile(settings.path2fieldtrip,'template','anatomy','single_subj_T1_1mm.nii'));
template_grid = importdata(settings.path2template_grid);
template_grid = ft_convert_units(template_grid,'mm'); 
atlas         = ft_read_atlas(fullfile(settings.path2fieldtrip,'template','atlas','aal','ROI_MNI_V4.nii')); % mm

sources_sym      = cell(1,n_subj);
sources_nosym    = cell(1,n_subj);
sources_vec      = cell(1,n_subj);
sources_sca_mean = cell(1,n_subj);
subjectnames     = cell(1,n_subj);
% In mni coordinates
sources_pos      = zeros(2,3,n_subj); % 2 dipoles x 3 coordinates x subjects

for sub_idx = 1:n_subj
    subject               = sprintf('sub-%02d',subjects(sub_idx));
    subjectnames{sub_idx} = subject;

    data = importdata(fullfile(settings.path2project,'derivatives',subject,'chirp',sprintf('%s_dipolfits.mat',subject))); 

    % Select correct data based on channel selection
    switch chan2plot
        case 'megmag'
            sources_sym{sub_idx}      = data.source_sym_mag;
            sources_nosym{sub_idx}    = data.source_nosym_mag;
            sources_vec{sub_idx}      = data.source_vec_mag;
            sources_sca_mean{sub_idx} = data.source_sca_mean_mag;
        case 'megplanar'
            sources_sym{sub_idx}      = data.source_sym_grad;
            sources_nosym{sub_idx}    = data.source_nosym_grad;
            sources_vec{sub_idx}      = data.source_vec_grad;
            sources_sca_mean{sub_idx} = data.source_sca_mean_grad;
    end
    clear data

    % Rescale units - optional for plotting (Am -> nAm)
    %--------------------------------------------------
    sources_vec{sub_idx}.dip.mom = 10^9*sources_vec{sub_idx}.dip.mom;
    sources_sca_mean{sub_idx}    = 10^9*sources_sca_mean{sub_idx};

    % Add source locations in mni space
    %----------------------------------
    sourcemodel = importdata(fullfile(settings.path2derivatives,subject,'chirp',sprintf('%s_sourcemodel-volumetric.mat',subject)));
    sourcemodel = ft_convert_units(sourcemodel, 'm'); 

    idx                      = dsearchn(sourcemodel.pos,sources_nosym{sub_idx}.dip.pos);
    sources_pos(:,:,sub_idx) = template_grid.pos(idx,:);
    clear sourcemodel

    fprintf('%s loaded.\n',subject)
end

%% Compute mean dipole positions and locations
%--------------------------------------------------------------------------
% Detect bad subjects for dipole fitting
[bad_subjects, mapping, idx_good_subject] = check_dipolefitting(sources_pos,subjectnames);
n_subj_remaining = sum(idx_good_subject);

dipole_positions_left        = zeros(n_subj_remaining,3);
dipole_positions_right       = zeros(n_subj_remaining,3);
dipole_moments_left          = zeros(n_subj_remaining,3,size(sources_nosym{1}.dip.mom,2));
dipole_moments_right         = zeros(n_subj_remaining,3,size(sources_nosym{1}.dip.mom,2));
dipole_timecourses_vec_left  = zeros(n_subj_remaining,3,length(sources_vec{1}.time));
dipole_timecourses_vec_right = zeros(n_subj_remaining,3,length(sources_vec{1}.time));
dipole_timecourses_sca_left  = zeros(n_subj_remaining,length(sources_vec{1}.time));
dipole_timecourses_sca_right = zeros(n_subj_remaining,length(sources_vec{1}.time));

counter = 1;
for sub_idx = find(idx_good_subject)
 
    if mapping(sub_idx, 1)==1 % first index is left side
        idx_left1  = 1;
        idx_right1 = 2;
        idx_left2  = 1:3;
        idx_right2 = 4:6;
    else 
        idx_left1  = 2;
        idx_right1 = 1;
        idx_left2  = 4:6;
        idx_right2 = 1:3;
    end

    % Append dipole positions
    dipole_positions_left(counter,:)  = sources_pos(idx_left1,:,sub_idx);
    dipole_positions_right(counter,:) = sources_pos(idx_right1,:,sub_idx);
    % Append dipole moments
    dipole_moments_left(counter,:,:)  = sources_nosym{sub_idx}.dip.mom(idx_left2,:);
    dipole_moments_right(counter,:,:) = sources_nosym{sub_idx}.dip.mom(idx_right2,:);
    % Append dipole timecourse - vector
    dipole_timecourses_vec_left(counter,:,:)  = sources_vec{sub_idx}.dip.mom(idx_left2,:);
    dipole_timecourses_vec_right(counter,:,:) = sources_vec{sub_idx}.dip.mom(idx_right2,:);
    % Append dipole timecourse - scalar
    dipole_timecourses_sca_left(counter,:)  = sources_sca_mean{sub_idx}(idx_left1,:);
    dipole_timecourses_sca_right(counter,:) = sources_sca_mean{sub_idx}(idx_right1,:);

    counter = counter + 1;
end

% Compute mean values
%--------------------
avg_dipole_position       = [squeeze(mean(dipole_positions_left,1));squeeze(mean(dipole_positions_right,1))];
avg_dipole_moment         = [squeeze(mean(dipole_moments_left,1));squeeze(mean(dipole_moments_right,1))];     
avg_dipole_timecourse_vec = [squeeze(mean(dipole_timecourses_vec_left,1));squeeze(mean(dipole_timecourses_vec_right,1))]; 
avg_dipole_timecourse_sca = [mean(dipole_timecourses_sca_left,1);mean(dipole_timecourses_sca_right,1)]; 
% Standard error of mean
avg_dipole_timecourse_vec_sem = [squeeze(std(dipole_timecourses_vec_left,1));squeeze(std(dipole_timecourses_vec_right,1))]./sqrt(n_subj_remaining);
avg_dipole_timecourse_sca_sem = [std(dipole_timecourses_sca_left,1);std(dipole_timecourses_sca_right,1)]./sqrt(n_subj_remaining);

%% Plot
%--------------------------------------------------------------------------

% plot 3D dipoles
cmap = distinguishable_colors(n_subj_remaining);
figure
hold on
counter = 1;
for sub_idx = find(idx_good_subject)
    ft_plot_dipole(sources_pos(1,:,sub_idx), mean(sources_nosym{sub_idx}.dip.mom(1:3,:),2), 'color', cmap(counter,:),  'unit', 'mm','alpha',0.5)
    ft_plot_dipole(sources_pos(2,:,sub_idx), mean(sources_nosym{sub_idx}.dip.mom(4:6,:),2), 'color', cmap(counter,:),  'unit', 'mm','alpha',0.5)
    counter = counter + 1;
end
ft_plot_dipole(avg_dipole_position(1,:), mean(avg_dipole_moment(1:3,:),2), 'color', [0.6350, 0.0780, 0.1840], 'unit', 'mm','thickness', 5, 'length', 20, 'diameter', 10)
ft_plot_dipole(avg_dipole_position(2,:), mean(avg_dipole_moment(4:6,:),2), 'color', [0.6350, 0.0780, 0.1840], 'unit', 'mm','thickness', 5, 'length', 20, 'diameter', 10)

% Position of crosshair - Heschl_R or mean dipole position (center x-axis)
%--------------------------------------------------------------------------
cfg         = [];
cfg.atlas   = atlas;
cfg.roi     = {'Heschl_R'};
mask_heschl = ft_volumelookup(cfg, atlas);
% Compute mean
linear_indices                = find(mask_heschl);
[row_idx, col_idx, slice_idx] = ind2sub(size(mask_heschl), linear_indices); % Convert linear indices to 3D subscripts (i, j, k voxel coordinates)
voxel_coords                  = [row_idx, col_idx, slice_idx, ones(length(linear_indices), 1)]'; % Combine these into a matrix where each row is a voxel coordinate (i,j,k), add 1 for homogeneous coordinates
coords_mni_homogeneous        = atlas.transform * voxel_coords; % 4. Transform voxel coordinates to mni coordinates
coords_mni                    = coords_mni_homogeneous(1:3, :)'; % Transpose to N x 3 matrix
pos                           = mean(coords_mni, 1); % Calculate the mean position

% or use centered dipole position
% pos    = mean(avg_dipole_position,1);
% pos(1) = 0; % center x-axis
%--------------------------------------------------------------------------

ft_plot_slice(template_mri.anatomy, 'transform', template_mri.transform, 'location', pos, 'orientation', [1 0 0], 'resolution', 0.1)
ft_plot_slice(template_mri.anatomy, 'transform', template_mri.transform, 'location', pos, 'orientation', [0 1 0], 'resolution', 0.1)
ft_plot_slice(template_mri.anatomy, 'transform', template_mri.transform, 'location', pos, 'orientation', [0 0 1], 'resolution', 0.1)

ft_plot_crosshair(pos, 'color', [1 1 1]/2);
axis tight
axis off
% view(0, 90) % axial
% view(0, 0) % coronal
view(90, 0) % sagittal

cfg         = [];
cfg.atlas   = atlas;
cfg.output  = 'multiple';
cfg.roi     = avg_dipole_position;
labels      = ft_volumelookup(cfg, atlas);
[~, indx]   = max(labels(1).count);
label_left  = labels(1).name(indx);
[~, indx]   = max(labels(2).count);
label_right = labels(2).name(indx);

fprintf('Mean dipole position left: %s\n',label_left{1})
fprintf('Mean dipole position right: %s\n',label_right{1})

%% (3) Visualize dipole time courses
%--------------------------------------------------------------------------

%% (3.1) All moments, xyz-directions and both hemispheres
%--------------------------------------------------------------------------

title_font_size  = 24; % Font size for subplot titles
label_font_size  = 24; % Font size for x and y axis labels
axis_font_size   = 24; % title_font_size  
legend_font_size = 20;

minval = min(avg_dipole_timecourse_vec-avg_dipole_timecourse_vec_sem,[],'all');
maxval = max(avg_dipole_timecourse_vec+avg_dipole_timecourse_vec_sem,[],'all');

name    = {'Left Hemisphere','Right Hemisphere'};
timevec = sources_vec{1}.time*1000;  
cmap    = {"r","g","b"}; %'x' 'y' 'z'
% Loop over both hemispheres
figure('Name','Freely oriented dipole (fixed positions, loose orientation','Color','w')
axisvec = horzcat(time2plot,[minval,maxval]);

for lr_idx = 1:2
    idxs2plot = (lr_idx-1)*3+(1:3); % 1:3 or 4:6
    
    sp(lr_idx) = subplot(2,1,lr_idx); 
    hold on
    % add patch for dipole fit timewindow
    patch([timewin_dipolefit(1), timewin_dipolefit(2), timewin_dipolefit(2), timewin_dipolefit(1)], ...
          [maxval, maxval, minval, minval], ...
                  [0.8,0.8,0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    
    % Plot dipole time courses and standard errors
    patch([timevec,fliplr(timevec)],[avg_dipole_timecourse_vec(idxs2plot(1),:)+avg_dipole_timecourse_vec_sem(idxs2plot(1),:),fliplr(avg_dipole_timecourse_vec(idxs2plot(1),:)-avg_dipole_timecourse_vec_sem(idxs2plot(1),:))],cmap{1},'FaceAlpha',0.1,'LineStyle','none','HandleVisibility','off');
    patch([timevec,fliplr(timevec)],[avg_dipole_timecourse_vec(idxs2plot(2),:)+avg_dipole_timecourse_vec_sem(idxs2plot(2),:),fliplr(avg_dipole_timecourse_vec(idxs2plot(2),:)-avg_dipole_timecourse_vec_sem(idxs2plot(2),:))],cmap{2},'FaceAlpha',0.1,'LineStyle','none','HandleVisibility','off');
    patch([timevec,fliplr(timevec)],[avg_dipole_timecourse_vec(idxs2plot(3),:)+avg_dipole_timecourse_vec_sem(idxs2plot(3),:),fliplr(avg_dipole_timecourse_vec(idxs2plot(3),:)-avg_dipole_timecourse_vec_sem(idxs2plot(3),:))],cmap{3},'FaceAlpha',0.1,'LineStyle','none','HandleVisibility','off');
    arrayfun(@(i) plot(timevec, avg_dipole_timecourse_vec(idxs2plot(i),:), '-', 'color',cmap{i},'LineWidth',2), 1:3) 

    if lr_idx ==2
        xlabel('t / ms')
    end
    ylabel('dipole moment / nAm'); 
    legend(sp(1),{'timewindow dipolfit','x', 'y', 'z'},'Location','southwest', 'FontSize', legend_font_size);
    axis(axisvec) 
    grid on;
    grid minor;
    box on; 
    title(name{lr_idx}, 'FontSize', title_font_size)

    xticks(time2plot(1):100:time2plot(2)); % Set x-axis ticks in steps of 100 ms
    set(gca, 'FontSize', axis_font_size); % Set font size for axis values
    
end
sgtitle(sprintf('%s',chan2plot),'FontSize', title_font_size+4)

% Link Y-axes of top plots for consistent scaling
linkaxes(sp, 'y');


%% (3.2) Fixed dipole, mean dipolmoment orientation as orientation constraint
%--------------------------------------------------------------------------

minval = min(avg_dipole_timecourse_sca-avg_dipole_timecourse_sca_sem,[],'all');
maxval = max(avg_dipole_timecourse_sca+avg_dipole_timecourse_sca_sem,[],'all');

name    = {'left hemisphere','right hemisphere'};
timevec = sources_vec{1}.time*1000;  
% Loop over both hemispheres
figure('Name','Freely oriented dipole (fixed positions, loose orientation','Color','w')
axisvec   = horzcat(time2plot,[minval,maxval]);
% line styles
ls = {'-','--'};

hold on
% add patch for dipole fit timewindow
patch([timewin_dipolefit(1), timewin_dipolefit(2), timewin_dipolefit(2), timewin_dipolefit(1)], ...
      [maxval, maxval, minval, minval], ...
       'y', 'EdgeColor', 'none', 'FaceAlpha', 0.2);
for lr_idx = 1:2
   
    % sp(lr_idx) = subplot(2,1,lr_idx);     
    % Plot dipole time courses and standard errors
    patch([timevec,fliplr(timevec)],[avg_dipole_timecourse_sca(lr_idx,:)+avg_dipole_timecourse_sca_sem(lr_idx,:),fliplr(avg_dipole_timecourse_sca(lr_idx,:)-avg_dipole_timecourse_sca_sem(lr_idx,:))],[0.8,0.8,0.8],'FaceAlpha',0.4,'LineStyle','none','HandleVisibility','off');
    plot(timevec, avg_dipole_timecourse_sca(lr_idx, :), '-', 'color','k','LineWidth',2, 'LineStyle', ls{lr_idx})
    
end
xlabel('t / ms')
ylabel('dipole moment / nAm'); 
legend([{'dipolfit timewindow'},name],'Location','northeast');
axis(axisvec) 
grid on
grid minor
sgtitle(sprintf('%s',chan2plot))

%% Additional functions
%------------------------------------------------------------------------

function [bad_subjects, hemisphere_map, is_good_subject] = check_dipolefitting(dipole_positions,subject_ids)
%-----------------------------------------------------------------------
% Till Habersetzer, 24.06.2025
% Communication Acoustics, CvO University Oldenburg
% till.habersetzer@uol.de 
%
% CHECK_DIPOLEFITTING Performs quality control on group dipole fitting results.
%
%   This function assesses the quality of bilateral dipole fits from a group
%   of subjects based on two criteria:
%   1.  Anatomical Plausibility: It rejects subjects where both dipoles are
%       fitted in the same hemisphere.
%   2.  Group Consistency: It identifies subjects whose dipole locations are
%       spatial outliers compared to the rest of the group.
%
%   SYNTAX:
%   [bad_subjects, hemisphere_map, is_good_subject] = check_dipolefitting(dipole_positions, subject_ids)
%
%   INPUTS:
%   dipole_positions : A 3D numeric array of dipole positions with dimensions
%                      [2 x 3 x n_subjects], where the 2nd dimension
%                      represents the [X, Y, Z] coordinates. The x-axis is
%                      assumed to separate the hemispheres (left: x<0, right: x>0).
%   subject_ids      : A cell array of strings containing the identifier for
%                      each subject, matching the 3rd dimension of 'dipole_positions'.
%
%   OUTPUTS:
%   bad_subjects     : A cell array of identifiers for subjects who failed
%                      one or more quality control checks.
%   hemisphere_map   : An [n_subjects x 2] array that maps the original
%                      dipole order to hemispheres. For a given subject (row),
%                      [1, 2] means dipole 1 is Left/dipole 2 is Right.
%                      [2, 1] means dipole 1 is Right/dipole 2 is Left.
%                      Rows for bad subjects contain NaN.
%   is_good_subject  : A logical row vector indicating which subjects passed
%                      all checks (1 = good, 0 = bad).
%-----------------------------------------------------------------------

% 0.) Initialization 
%-----------------------------------------------------------------------
% dipole_positions = sources_pos;
% subject_ids      = subjectnames;

n_subjects = size(dipole_positions, 3);
if n_subjects ~= length(subject_ids)
    error('Number of subjects in dipole_positions and subject_ids must match.');
end

hemisphere_map = nan(n_subjects, 2);

% 1.) Check for Bilateral Dipoles 
%-----------------------------------------------------------------------
% Dipoles are considered unilateral if their x-coordinates have the same sign
x_coords          = squeeze(dipole_positions(:, 1, :));
is_bad_unilateral = (prod(x_coords, 1) > 0); 
    
% Initialize a logical index of subjects who are currently considered "good"
is_good_subject = ~is_bad_unilateral;

% 2.) Map Hemispheres and Sort Positions for good subjects 
%-----------------------------------------------------------------------
% Pre-allocate arrays to hold sorted dipole positions.
left_dipole_pos  = nan(n_subjects, 3);
right_dipole_pos = nan(n_subjects, 3);

% Only loop over subjects who passed the first check.
for sub_idx = find(is_good_subject)
    % Check the sign of the first dipole's x-coordinate to map it.
    if sign(dipole_positions(1, 1, sub_idx)) < 0 % Dipole 1 is in the left hemisphere (x<0)
       hemisphere_map(sub_idx, :) = [1, 2]; % [1,2] -> [L,R]
    else
       hemisphere_map(sub_idx, :) = [2, 1]; % [2,1] -> [R,L]
    end

    % Use the mapping to sort positions into left and right arrays.
    is_left_dipole  = (hemisphere_map(sub_idx, :) == 1);
    is_right_dipole = (hemisphere_map(sub_idx, :) == 2);
    
    left_dipole_pos(sub_idx, :)  = dipole_positions(is_left_dipole, :, sub_idx);
    right_dipole_pos(sub_idx, :) = dipole_positions(is_right_dipole, :, sub_idx);
end

% 3.) Identify Spatial Outliers using a Group Distance Metric 
%-----------------------------------------------------------------------
% Calculate the pairwise Euclidean distance between all subjects' dipoles.
distance_matrix_L = pdist2(left_dipole_pos, left_dipole_pos, 'euclidean');
distance_matrix_R = pdist2(right_dipole_pos, right_dipole_pos, 'euclidean');

% For each subject, sum their distances to all other subjects. A large sum
% indicates that the subject's dipole is far from the group cluster.
group_distance_L = sum(distance_matrix_L, 2);
group_distance_R = sum(distance_matrix_R, 2);

% Use the Median Absolute Deviation (MAD) method to find outliers in the
% group distance metric. This flags subjects whose dipoles are unusually
% far from the rest of the group. isoutlier ignores NaNs by default.
is_bad_outlier_L = isoutlier(group_distance_L);
is_bad_outlier_R = isoutlier(group_distance_R);

is_bad_outlier = is_bad_outlier_L | is_bad_outlier_R;

% 4.) Finalize bad subject list
%-----------------------------------------------------------------------
% A subject is bad if they failed the unilateral check OR the outlier check.
is_good_subject(is_bad_outlier) = false;

% Get the names of all subjects who failed any check.
bad_subjects = subject_ids(is_good_subject == false);

end