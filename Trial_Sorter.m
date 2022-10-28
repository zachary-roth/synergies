% Author: Zach Roth <zachary.roth@student.kuleuven.be>
% Created: 25-October-2022

% This script uses the meta data provided in the 'trials.xlsx' to
% reorganize the data into movement directories with data-type
% sub-directories. The files are copied into the new folders, leaving the
% original data intact.

close all
clear
clc

% HARD-CODE the data types
data_types = {'EMG','IK','ID'}';

% Select the file containing the session metadata
[meta_file, meta_path] = uigetfile('*.xlsx');

% Move to the folder containing the metadata
cd(meta_path)
mkdir('Movements')

% Import the metadata
meta = readtable(strcat(meta_path,meta_file));

% Get a list of the movements in the session
movements = unique(meta.movement);

% Create the Movement directory with data-type sub-directories
for m = 1:length(movements)
    dir_path = fullfile(meta_path,'Movements',movements{m});
    mkdir(dir_path)
    for s = 1:length(data_types)
        sub_path = fullfile(meta_path,'Movements',movements{m},data_types{s});
        mkdir(sub_path)
    end
end

% Sort the data files into their respective movement and dataype
% sub-directories

for m = 1:length(movements)
    trial_i = find(strcmp(movements(m),meta.movement)); %Get the trial INDICES of a specific movement
    trial_n = meta.ID(trial_i); %Get the trial NUMBERS of a specific movement
    % Get a list of the files in a data folder
    for d = 1:length(data_types)
        % Create a list of file names
        file_names = dir(fullfile(meta_path,data_types{d}));
        % Find the files that match the relevant trial number
        for f = 1:length(file_names)
            % Copy the relevant files into thier respective sub-directories
            for n = 1:length(trial_n)
                if contains(file_names(f).name,trial_n(n)) == 1
                    source = fullfile(meta_path,data_types{d},file_names(f).name);
                    destination = fullfile(meta_path,'Movements',movements{m},data_types{d},file_names(f).name);
                    copyfile(source,destination)
                else
                    continue
                end
            end
        end
    end
end

