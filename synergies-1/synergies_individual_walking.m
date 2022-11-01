
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
walking_dir = uigetdir();
tic
cd(walking_dir)

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
% Remove the 50Hz noise with Notch filter
fs = 1000;
fo = 50;  
q = 35; 
bw = (fo/(fs/2))/q;
[b,a] = iircomb(round(fs/fo),bw,'notch');

trials = fieldnames(emg_raw);
for trial = 1:length(trials)
    trial_num = sscanf(trials{trial},strcat("trial","%f"));
    for muscle = 2:17
        raw_EMG = table2array(emg_raw.(['trial' num2str(trial_num)])(:,muscle));
        emg_raw.(['trial' num2str(trial_num)])(:,muscle) = num2cell(filter(b,a,raw_EMG));
    end
    disp(['Trial',num2str(trial),' notch filtered'])
end

% Bandpass filter (20-400 Hz)
order = 4; 
cutoff = [20 400];
[b,a] = butter(order/2, cutoff/(0.5*fs));

trials = fieldnames(emg_raw);
for trial = 1:length(trials)
    trial_num = sscanf(trials{trial},strcat("trial","%f"));
    for muscle = 2:17
        rawEMG_notch = table2array(emg_raw.(['trial' num2str(trial_num)])(:,muscle));
        emg_raw.(['trial' num2str(trial_num)])(:,muscle) = num2cell(filtfilt(b,a,rawEMG_notch));
    end
    disp(['Trial',num2str(trial),' band-pass filtered'])
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

% Smooth EMG data (15 Hz low pass filter)
order = 4;
cutoff = 15;
[b, a] = butter(order, cutoff/(0.5*fs));

trials = fieldnames(emg_raw);
for trial = 1:length(trials)
    trial_num = sscanf(trials{trial},strcat("trial","%f"));
    for muscle = 2:17
        rectEMG = table2array(emg_raw.(['trial' num2str(trial_num)])(:,muscle));
        emg_raw.(['trial' num2str(trial_num)])(:,muscle) = num2cell(filtfilt(b,a,rectEMG));
    end
    disp(['Trial',num2str(trial),' smoothed'])
end

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

% Visualize EMG

muscles = {'RREF','RVAL','RBIF','RMEH','RTIA','RGAS','RSOL','RGLU',...
    'LREF','LVAL','LBIF','LMEH','LTIA','LGAS','LSOL','LGLU'};

t = tiledlayout(2,4);

X = (0:100)';

for muscle = 1:8
    nexttile
    for trial = 1:length(trials)
        trial_num = sscanf(trials{trial},strcat("trial","%f"));
        hold on
        Y1 = emg_proc.(['trial' num2str(trial_num)]){:,muscles{muscle}};
        Y2 = emg_proc.(['trial' num2str(trial_num)]){:,muscles{muscle+8}};
        plot(X,Y1,'LineWidth',0.5,'Color',[0 1 0 0.5])
        plot(X,Y2,'LineWidth',0.5,'Color',[1 0 0 0.5])
    end
    ylim([0 0.3])
    hold off
end

toc
