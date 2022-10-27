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
% Point to the folder containing the EMG data
emg_dir = uigetdir();

tic

cd(emg_dir)

% Get a list of EMG files
emg_files = dir(fullfile(emg_dir, '*.xlsx'));
emg_files = {emg_files.name}';

% Collate the EMG files into a single structure
for file = 1:length(emg_files)
    trial_num = sscanf(emg_files{file},strcat("EMG_raw_","%f",".xlsx"));
    data = readtable(emg_files{file});
    if isempty(data) == 1
       disp(['File ',num2str(file),' empty'])
        continue
    else
        emg_raw.(['trial' num2str(trial_num)]) = data;
        disp(['File ',num2str(file),' imported'])
    end
end

% 2. Filtering

% Create a HIGH-pass 4-order Butterworth Filter
fc = 10; % cutoff frequency
fs = 1000; % sampling rate
n  = 4; % nth-order
Wn = fc/(fs/2); % Cut off frequency
ftype = 'high'; % Filter type

[b,a] = butter(n,Wn,ftype);

% HIGH-pass filter EMG data
trials = fieldnames(emg_raw);
for trial = 1:length(trials)
    trial_num = sscanf(trials{trial},strcat("trial","%f"));
    for muscle = 2:17
        x = table2array(emg_raw.(['trial' num2str(trial_num)])(:,muscle));
        emg_raw.(['trial' num2str(trial_num)])(:,muscle) = num2cell(filtfilt(b,a,x));
    end
    disp(['Trial',num2str(trial),' high-pass filtered'])
end

% Create a LOW-pass 4-order Butterworth Filter
fc = 50; % cutoff frequency
fs = 1000; % sampling rate
n  = 4; % nth-order
Wn = fc/(fs/2); % Cut off frequency
ftype = 'low'; % Filter type

[b,a] = butter(n,Wn,ftype);

% LOW-pass filter EMG data
trials = fieldnames(emg_raw);
for trial = 1:length(trials)
    trial_num = sscanf(trials{trial},strcat("trial","%f"));
    for muscle = 2:17
        x = table2array(emg_raw.(['trial' num2str(trial_num)])(:,muscle));
        emg_raw.(['trial' num2str(trial_num)])(:,muscle) = num2cell(filtfilt(b,a,x));
    end
    disp(['Trial',num2str(trial_num),' low-pass filtered'])
end

% Rectify EMG data
for trial = 1:length(trials)
    trial_num = sscanf(trials{trial},strcat("trial","%f"));
    for muscle = 2:17
        filt_emg = table2array(emg_raw.(['trial' num2str(trial_num)])(:,muscle));
        emg_raw.(['trial' num2str(trial_num)])(:,muscle) = num2cell(abs(filt_emg));
    end
    disp(['Trial',num2str(trial_num),' rectified'])
end

% Smooth EMG data


% Resample EMG data
proc_cols = {'Percent','RREF','RVAL','RBIF','RMEH','RTIA','RGAS','RSOL','RGLU','LREF','LVAL','LBIF','LMEH','LTIA','LGAS','LSOL','LGLU'};

for trial = 1:length(trials)
    trial_num = sscanf(trials{trial},strcat("trial","%f"));
    emg_proc.(['trial' num2str(trial_num)]) = array2table(zeros(101,17), 'VariableNames',proc_cols);
    emg_proc.(['trial' num2str(trial_num)])(:,1) = num2cell((0:100)');
    for muscle = 2:17
        rect_emg = table2array(emg_raw.(['trial' num2str(trial_num)])(:,muscle));
        % sample values
        v = table2array(emg_raw.(['trial' num2str(trial_num)])(:,muscle));
        % sample points
        x = (1:1:length(rect_emg))';
        % query point
        xq = (1:(length(rect_emg)/101):length(rect_emg))';
        % resample the kinetic data
        emg_proc.(['trial' num2str(trial_num)])(:,muscle) = num2cell(interp1(x,v,xq));
    end
    disp(['Trial',num2str(trial_num),' resampled'])
end


toc