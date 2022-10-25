%                                   _            __  
%                                  (_)          /  | 
%    ___ _   _ _ __   ___ _ __ __ _ _  ___  ___ `| | 
%   / __| | | | '_ \ / _ \ '__/ _` | |/ _ \/ __| | | 
%   \__ \ |_| | | | |  __/ | | (_| | |  __/\__ \_| |_
%   |___/\__, |_| |_|\___|_|  \__, |_|\___||___/\___/
%         __/ |                __/ |                 
%        |___/                |___/                           
% 
% Author: Zach Roth <zachary.roth@student.kuleuven.be>
% Created: 25-October-2022

clc
clear
close all

tic

%% EMG DATA
%--------------------------------------------------------------------------

% 1. Import Data

emg_dir = uigetdir();

tic

cd(emg_dir)

% Get list of EMG files
emg_files = dir(fullfile(emg_dir, '*.xlsx'));
emg_files = {emg_files.name}';

% Collate the EMG files into a single structure
for ii = 1:length(emg_files)
    trial_num = sscanf(emg_files{ii},strcat("EMG_raw_","%f",".xlsx"));
    data = readtable(emg_files{ii});
    if isempty(data) == 1
        continue
    else
        emg_data.(['trial' num2str(trial_num)]) = data;
    end
end

toc