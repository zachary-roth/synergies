%% Default example solve muscle redundancy 
% (as in DeGroote2016)

clear 
close all
clc

%% Input information
% Add here the paths of IK, ID , US and EMG data trials you want to work with
Misc.IKfile = {'C:\Users\zacharyroth\Data\synergies-data\Subject1\IK\IK_11.mot'};
Misc.IDfile = {'C:\Users\zacharyroth\Data\synergies-data\Subject1\ID\ID_IK_11.sto'};
Misc.model_path  = 'C:\Users\zacharyroth\Data\synergies-data\Subject1\Models\GEN_scaled.osim';
Misc.OutPath     = 'C:\Users\zacharyroth\Data\synergies-data\Subject1\Movements\gait_normal\Calculated-Activations\Gait2';                    % folder to store results

% Get start and end time of the different files
time=[2.66 3.48]; % Right stance phase (+50ms beginning and end of time interval, more details see manual and publication)

% Settings
Misc.DofNames_Input={'ankle_angle_r','knee_angle_r','hip_flexion_r','hip_adduction_r','hip_rotation_r'};    % select the DOFs you want to include in the optimization

% Plotter Bool: Boolean to select if you want to plot lots of output information of intermediate steps in the script
Misc.PlotBool = false;

% MRS Bool: Select if you want to run the generic muscle redundancy solver
Misc.MRSBool = 1;

% Validation Bool: Select if you want to run the muscle redundancy solver with the optimized parameters
Misc.ValidationBool = 1; 	% TO DO: we should report results of EMG driven simulation as well

% name output
Misc.OutName = 'Normal_Gait_';

%% Run muscle tendon estimator:
[Results,DatStore] = solveMuscleRedundancy(time,Misc);
