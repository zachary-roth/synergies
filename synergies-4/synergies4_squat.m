%   _____                   _       
%  /  ___|                 | |      
%  \ `--.  __ _ _   _  __ _| |_ ___ 
%   `--. \/ _` | | | |/ _` | __/ __|
%  /\__/ / (_| | |_| | (_| | |_\__ \
%  \____/ \__, |\__,_|\__,_|\__|___/
%            | |                    
%            |_|                    

% Author: Zach Roth <zachary.roth@student.kuleuven.be>
% Created: 28-October-2022

%% Import Data

% 1. Import Data
% Select the folder containing the squat data
squat_dir = uigetdir();
cd(squat_dir)

% Get a list of IK files
ik_dir = fullfile(squat_dir,'IK');
ik_files = dir(fullfile(ik_dir, '*.mot'));
ik_files = {ik_files.name}';

squats = struct();

for f = 1:length(ik_files)
    squats.IK.(strcat('squat_',num2str(f))) = ReadMotFile(fullfile(squat_dir,'IK',ik_files{f}));
end

%% Visualize the Knee Flexion Data
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