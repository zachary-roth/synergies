%{
   _____                   _       
  /  ___|                 | |      
  \ `--.  __ _ _   _  __ _| |_ ___ 
   `--. \/ _` | | | |/ _` | __/ __|
  /\__/ / (_| | |_| | (_| | |_\__ \
  \____/ \__, |\__,_|\__,_|\__|___/
            | |                    
            |_|                    
%}

%{ 
Author: Zach Roth <zachary.roth@student.kuleuven.be>
Created: 28-October-2022

Part 1: Movement Events
Part 2: Process EMG Data

%}

close all
clear
clc

%% 1. Movement Events

% 1. Import Data
% Select the folder containing the squat data
squat_dir = uigetdir();
cd(squat_dir)

% Get a list of IK files
ik_dir = fullfile(squat_dir,'IK');
ik_files = dir(fullfile(ik_dir, '*.mot'));
ik_files = {ik_files.name}';

% Initialize a structure to hold all the data called squats
squats = struct();

% Import the IK files into the squats structure
for f = 1:length(ik_files)
    squats.IK.(strcat('squat_',num2str(f))) = ReadMotFile(fullfile(squat_dir,'IK',ik_files{f}));
    disp(['IK ',num2str(f),' imported'])
end

% Visualize the Knee Flexion Data
% Get the time column index
time_col = find(strcmp('time',squats.IK.squat_1.names));

% Get the column index of the right knee angle
r_knee_col = find(strcmp('knee_angle_r',squats.IK.squat_1.names));

% Get the column index of the left knee angle
l_knee_col = find(strcmp('knee_angle_l',squats.IK.squat_1.names));

max_flexion = zeros(length(ik_files),2);

t = tiledlayout('flow');
t.Title.String = 'Knee Flexion in the Squat Movement';
t.Title.FontWeight = 'bold';

for tile = 1:length(ik_files)
    nexttile
    
    X = squats.IK.(strcat('squat_',num2str(tile))).data(:,time_col);
    Y1 = squats.IK.(strcat('squat_',num2str(tile))).data(:,r_knee_col);
    Y2 = squats.IK.(strcat('squat_',num2str(tile))).data(:,l_knee_col);
    
    % Get the indices where max knee flexion occurs
    r_max_i = find(Y1 == max(Y1));
    l_max_i = find(Y2 == max(Y2));

    % Get the time where max knee flexion occurs
    r_max_t = squats.IK.(strcat('squat_',num2str(tile))).data(r_max_i,time_col);
    l_max_t = squats.IK.(strcat('squat_',num2str(tile))).data(l_max_i,time_col);

    % Save the max flexion indices
    squats.IK.(strcat('squat_',num2str(tile))).end_index = max(l_max_i,r_max_i);
    squats.IK.(strcat('squat_',num2str(tile))).end_time = max(l_max_t,r_max_t);

    hold on
    plot(X,Y1,'LineWidth',1.5,'color',[0 1 0 .25])
    plot(X,Y2,'LineWidth',1.5,'color',[1 0 0 .25])
    xline(r_max_t,'LineStyle','--','color',[0 1 0 .25])
    xline(l_max_t,'LineStyle','--','color',[1 0 0 .25])
    grid on
    xlim([0 6])
    ylim([-120 0])
    xlabel('Time [s]')
    ylabel(['Knee Flexion [' char(176) ']'])
    hold off
end
leg = legend('Right Knee Flexion','Left Knee Flexion', 'Max R Knee Flexion','Max L Knee Flexion','Orientation', 'Horizontal');
leg.Layout.Tile = 'north';

savefig(fullfile(squat_dir,'Figures','Knee_Flexion.fig'))

%% Process EMG Data

% Import EMG Data
% Get a list of IK files
EMG_dir = fullfile(squat_dir,'EMG');
EMG_files = dir(fullfile(EMG_dir, '*.xlsx'));
EMG_files = {EMG_files.name}';

% Import the EMG files into the squats structure
for f = 1:length(EMG_files)
    data = readtable(fullfile(squat_dir,'EMG',EMG_files{f}));
    squats.EMG.names = data.Properties.VariableNames;
    squats.EMG.EMG_raw.(strcat('squat_',num2str(f))) = table2array(data);
    disp(['EMG ',num2str(f),' imported'])
end

% Trim the trials
% The duration of one squat is defined as the beginning of the trial until
% the instant of maximum knee flexion in EITHER the right or the left knee.
% Data collected after the instant of max flexion is discarded.

for f = 1:length(EMG_files)
    i_first = find(squats.EMG.EMG_raw.(strcat('squat_',num2str(f)))(:,2),1,'first');
    i_last = find(squats.EMG.EMG_raw.(strcat('squat_',num2str(f)))(:,1) == squats.IK.(strcat('squat_',num2str(f))).end_time);
    squats.EMG.EMG_trim.(strcat('squat_',num2str(f))) = squats.EMG.EMG_raw.(strcat('squat_',num2str(f)))(i_first:i_last,:);
    disp(['EMG squat_',num2str(f),' trimmed'])
end

% Filtering
% Notch filter: Remove the 50Hz noise
fs = 1000;
fo = 50;  
q = 35; 
bw = (fo/(fs/2))/q;
[b,a] = iircomb(round(fs/fo),bw,'notch');


for f = 1:length(EMG_files)
    squats.EMG.EMG_notch.(strcat('squat_',num2str(f))) = zeros(size(squats.EMG.EMG_trim.(strcat('squat_',num2str(f)))));
    squats.EMG.EMG_notch.(strcat('squat_',num2str(f)))(:,1) = squats.EMG.EMG_trim.(strcat('squat_',num2str(f)))(:,1);
    for m = 2:length(squats.EMG.names)
        EMG_trim = squats.EMG.EMG_trim.(strcat('squat_',num2str(f)));
        squats.EMG.EMG_notch.(strcat('squat_',num2str(f)))(:,m) = filter(b,a,EMG_trim(:,m));
    end
    disp(['EMG squat_',num2str(f),' notch filtered'])
end

% Bandpass filter (20-400 Hz)
order = 4; 
cutoff = [20 400];
[b,a] = butter(order/2, cutoff/(0.5*fs));

for f = 1:length(EMG_files)
    squats.EMG.EMG_bandpass.(strcat('squat_',num2str(f))) = zeros(size(squats.EMG.EMG_notch.(strcat('squat_',num2str(f)))));
    squats.EMG.EMG_bandpass.(strcat('squat_',num2str(f)))(:,1) = squats.EMG.EMG_notch.(strcat('squat_',num2str(f)))(:,1);
    for m = 2:length(squats.EMG.names)
        EMG_notch = squats.EMG.EMG_notch.(strcat('squat_',num2str(f)));
        squats.EMG.EMG_bandpass.(strcat('squat_',num2str(f)))(:,m) = filtfilt(b,a,EMG_notch(:,m));
    end
    disp(['EMG squat_',num2str(f),' band-pass filtered'])
end

% Rectification
for f = 1:length(EMG_files)
    squats.EMG.EMG_rect.(strcat('squat_',num2str(f))) = abs(squats.EMG.EMG_bandpass.(strcat('squat_',num2str(f))));
    disp(['EMG squat_',num2str(f),' rectified'])
end

% 15 Hz low pass filter: Smoothing
order = 4;
cutoff = 15;
[b, a] = butter(order, cutoff/(0.5*fs));

for f = 1:length(EMG_files)
    squats.EMG.EMG_smooth.(strcat('squat_',num2str(f))) = zeros(size(squats.EMG.EMG_rect.(strcat('squat_',num2str(f)))));
    squats.EMG.EMG_smooth.(strcat('squat_',num2str(f)))(:,1) = squats.EMG.EMG_rect.(strcat('squat_',num2str(f)))(:,1);
    for m = 2:length(squats.EMG.names)
        EMG_rect = squats.EMG.EMG_rect.(strcat('squat_',num2str(f)));
        squats.EMG.EMG_smooth.(strcat('squat_',num2str(f)))(:,m) = filtfilt(b,a,EMG_rect(:,m));
    end
    disp(['EMG squat_',num2str(f),' smoothed'])
end

% Resample
for f = 1:length(EMG_files)
    squats.EMG.EMG_resamp.(strcat('squat_',num2str(f))) = zeros(101,length(squats.EMG.names));
    squats.EMG.EMG_resamp.(strcat('squat_',num2str(f)))(:,1) = (0:100)';
    for m = 2:length(squats.EMG.names)
        EMG_smooth = squats.EMG.EMG_smooth.(strcat('squat_',num2str(f)));
        % sample values
        v = EMG_smooth(:,m);
        % sample points
        x = (1:1:length(EMG_smooth))';
        % query point
        xq = (1:(length(EMG_smooth)/101):length(EMG_smooth))';
        % resample the kinetic data
        squats.EMG.EMG_resamp.(strcat('squat_',num2str(f)))(:,m) = interp1(x,v,xq);
    end
    disp(['EMG squat_',num2str(f),' resampled'])
end
