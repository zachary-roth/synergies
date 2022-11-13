function [EMG_max_all] = find_EMG_max(EMG_dir,EMG_files)
% Find the maximum Activation values across all trials

% Import the first EMG file
EMG_all = readmatrix(fullfile(EMG_dir,EMG_files{1}));

% Vertically concatenate all the EMG files into a single array
for f = 2:length(EMG_files)
    trial = readmatrix(fullfile(EMG_dir,EMG_files{f})); % Import the next trial
    EMG_all = [EMG_all;trial]; % Vert Cat
end

EMG_all = abs(EMG_all); % Rectify the EMG data
EMG_max_all = max(EMG_all); % Find the max activation for each column (muscle)