
clear
close all
clc

% Import Opensim
import org.opensim.modeling.*

% Import CasADi
addpath('/Users/zacharyroth/Documents/MATLAB/casadi-osx-matlabR2015a-v3.5.5')
import casadi.*

% Add the muscle redundancy solver to the path
addpath(genpath('/Users/zacharyroth/Repos/MuscleRedundancySolver'))

DataPath = '/Users/zacharyroth/Data/synergies-data/Subject1/';

Misc.model_path  = {fullfile(DataPath,'Models/GEN_scaled.osim')};
Misc.IKfile = {fullfile(DataPath,'Movements/gait_normal/IK/IK_06.mot')};
Misc.IDfile = {fullfile(DataPath,'Movements/gait_normal/ID/ID_IK_06.sto')};

% Right stance phase (+50ms beginning and end of time interval, more details see manual and publication)
time = [2.835 3.675];

% select the DOFs you want to include in the optimization
Misc.DofNames_Input = {'ankle_angle_r','knee_angle_r','hip_flexion_r','hip_adduction_r','hip_rotation_r'};

% folder to store results
Misc.OutPath = {fullfile(DataPath,'Movements/gait_normal/Calculated-Activations')};

% Plotter Bool: Boolean to select if you want to plot lots of output information of intermediate steps in the script
Misc.PlotBool = 1;
% name output
Misc.OutName = 'Gait_';

[Results,DatStore] = solveMuscleRedundancy(time,Misc);
