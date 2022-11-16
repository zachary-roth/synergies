%{
_____/\\\\\\\\\\\\____________________________________        
 ___/\\\//////////_____________________________________       
  __/\\\_____________________________/\\\_____/\\\______      
   _\/\\\____/\\\\\\\__/\\\\\\\\\____\///___/\\\\\\\\\\\_     
    _\/\\\___\/////\\\_\////////\\\____/\\\_\////\\\////__    
     _\/\\\_______\/\\\___/\\\\\\\\\\__\/\\\____\/\\\______   
      _\/\\\_______\/\\\__/\\\/////\\\__\/\\\____\/\\\_/\\__  
       _\//\\\\\\\\\\\\/__\//\\\\\\\\/\\_\/\\\____\//\\\\\___ 
        __\////////////_____\////////\//__\///______\/////____

Author: Zach Roth <zachary.roth@student.kuleuven.be>
Created: 31-October-2022

This script is functionally identical to the calculated activations section
synergies1_GaitData.m. It does not include the EMG section (which was
already calculated) and differs in the way the data is imported (the file
naming conventions and metavdata structure was different than in
synergies1).

%}

close all
clear
clc

%%   Select the Subject1 data folder
disp('Select the Subject1 data folder')
dataDir = uigetdir('','Select the Subject1 data folder');
cd(dataDir)

% Delete the old Results Folder
if exist("Results","dir") == 7
    rmdir("Results",'s')
end

% Create a New Results folder and subfolders
mkdir("Results")
mkdir(fullfile("Results","Calculated_Activations"))
mkdir(fullfile("Results","Calculated_Activations","L"))
mkdir(fullfile("Results","Calculated_Activations","R"))
mkdir(fullfile("Results","Figures"))


%% Meta Data
%Import the meta data
meta = readtable(fullfile(dataDir,'timepoints_functional_shorter_2.xlsx'));

% Get a list of indices for all the normal gait trials
normal_gaits = find(strcmp(meta.movement,'gait_normal'));

% Get a list of indices for all the valid trials (recorded ICs)
valid_ICs = find(~isnan(meta.IC1));

% Filter out the invalid normal gait trials
gait_indices = intersect(normal_gaits,valid_ICs);

% Write the meta data to the trials structure
GaitData.meta.ID = meta.Trial(gait_indices);
GaitData.meta.side = meta.side(gait_indices);
GaitData.meta.IC1 = meta.IC1(gait_indices);
GaitData.meta.IC2 = meta.IC2(gait_indices);

% %% EMG
% % Get a list of EMG files
% EMG_dir = fullfile(dataDir,'EMG');
% EMG_files = dir(fullfile(EMG_dir, '*.xlsx'));
% EMG_files = {EMG_files.name};
% 
% % Find the max EMG activation values across all trials
% disp('Finding Max EMG Values')
% EMG_max_all = find_EMG_max(EMG_dir,EMG_files);
% 
% % Get the full list of EMG muscle names
% EMG_names = readtable(fullfile(EMG_dir,EMG_files{1}));
% EMG_names = EMG_names.Properties.VariableNames;
% 
% % Get the column names and indices for Left and Right muscles
% EMG_muscles_right_idx = find(strncmpi('R',EMG_names,1));
% EMG_muscles_right = EMG_names(EMG_muscles_right_idx);
% 
% EMG_muscles_left_idx = find(strncmpi('L',EMG_names,1));
% EMG_muscles_left = EMG_names(EMG_muscles_left_idx);
% 
% % Combine the muscles into a 1x2 array
% muscles = [{EMG_muscles_left}, {EMG_muscles_right}];
% 
% % Import and Process the EMG files into the GaitData Structure
% 
% % Initialize the left and right indices
% li = 1;
% ri = 1;
% 
% % Loop over the valid gait trials
% for f = 1:length(gait_indices)
%     % Read the gait side from the meta file
%     side = GaitData.meta.side{f};
%     % Assign an s index variable for slicing muscle names and writing to
%     % output structures
%     if side == 'L'
%         s = 1;
%     else
%         s = 2;
%     end
%     
%     % Read in the raw EMG files
%     rawFile = strcat('EMG_raw',GaitData.meta.ID{f},'.xlsx');
%     raw = readmatrix(fullfile(EMG_dir,rawFile));
% 
%     % Get the indices of the initial contacts (heel strikes)
%     IC1 = GaitData.meta.IC1(f);
%     IC1_idx = find(raw(:,1)==IC1);
%     IC2 = GaitData.meta.IC2(f);
%     IC2_idx = find(raw(:,1)==IC2);
%     
%     % Take the right muscles for right gaits and vice-versa
%     if GaitData.meta.side{f} == 'R'
%         trim = raw(IC1_idx:IC2_idx,EMG_muscles_right_idx);
%     elseif GaitData.meta.side{f} == 'L'
%         trim = raw(IC1_idx:IC2_idx,EMG_muscles_left_idx);
%     end
%     
%     % Notch filter: Remove the 50Hz noise
%     % Create the notch filter
%     fs = 1000;
%     fo = 50;
%     q = 35;
%     bw = (fo/(fs/2))/q;
%     [b,a] = iircomb(round(fs/fo),bw,'notch');
%     
%     % Initialize an output array
%     notch = zeros(size(trim));
%     
%     % Filter each column
%     for m = 1:length(muscles{s})
%         notch(:,m) = filter(b,a,trim(:,m));
%     end
% 
%     % Bandpass filter (20-400 Hz)
%     % Create the bandpass filter
%     order = 4;
%     cutoff = [20 400];
%     [b,a] = butter(order/2, cutoff/(0.5*fs));
%     
%     % Intialize an output array
%     bandPass = zeros(size(notch));
%     
%     % Filter each column
%     for m = 1:length(muscles{s})
%         bandPass(:,m) = filtfilt(b,a,notch(:,m));
%     end
% 
%     % Rectification
%     rect = abs(bandPass);
% 
%     % Smoothing (15 Hz low pass filter)
%     % Create the smoothing filter
%     order = 4;
%     cutoff = 15;
%     [b, a] = butter(order, cutoff/(0.5*fs));
%     
%     % Initialize an output array
%     lowPass = zeros(size(rect));
%     
%     % Filter each column
%     for m = 1:length(muscles{s})
%         lowPass(:,m) = filtfilt(b,a,rect(:,m));
%     end
% 
%     % Normalize
%     % Initialize an output array
%     norm = zeros(size(lowPass));
%     
%     % Get the correct left or right EMG max values 
%     if GaitData.meta.side{f} == 'R'
%         EMG_max = EMG_max_all(:,EMG_muscles_right_idx);
%     elseif GaitData.meta.side{f} == 'L'
%         EMG_max = EMG_max_all(:,EMG_muscles_left_idx);
%     end
%     
%     % Normalize each column
%     for m = 1:length(muscles{s})
%         norm(:,m) = lowPass(:,m)/EMG_max(m);
%     end
% 
%     % Resample
%     % Initialize an output array
%     resample = zeros(101,length(muscles{s}));
%     
%     % Specifiy the input arguments for the interp1 function
%     x = (1:1:length(norm))'; % sample points
%     xq = linspace(1,length(norm),101); % query point
%     
%     % Resample each column
%     for m = 1:length(muscles{s})
%         v = norm(:,m); % sample values
%         % resample the kinetic data
%         resample(:,m) = interp1(x,v,xq);
%     end
% 
%     % Store the processing steps in the GaitData structure    
%     if side == 'R'
%         GaitData.EMG.(side).raw.(strcat('raw',num2str(ri))) = raw;
%         GaitData.EMG.(side).trim.(strcat('trim',num2str(ri))) = trim;
%         GaitData.EMG.(side).notch.(strcat('notch',num2str(ri))) = notch;
%         GaitData.EMG.(side).bandPass.(strcat('bandPass',num2str(ri))) = bandPass;
%         GaitData.EMG.(side).rect.(strcat('rect',num2str(ri))) = rect;
%         GaitData.EMG.(side).lowPass.(strcat('lowPass',num2str(ri))) = lowPass;
%         GaitData.EMG.(side).norm.(strcat('norm',num2str(ri))) = norm;
%         GaitData.EMG.(side).resample.(strcat('resample',num2str(ri))) = resample;
%         ri = ri +1;
%     else
%         GaitData.EMG.(side).raw.(strcat('raw',num2str(li))) = raw;
%         GaitData.EMG.(side).trim.(strcat('trim',num2str(li))) = trim;
%         GaitData.EMG.(side).notch.(strcat('notch',num2str(li))) = notch;
%         GaitData.EMG.(side).bandPass.(strcat('bandPass',num2str(li))) = bandPass;
%         GaitData.EMG.(side).rect.(strcat('rect',num2str(li))) = rect;
%         GaitData.EMG.(side).lowPass.(strcat('lowPass',num2str(li))) = lowPass;
%         GaitData.EMG.(side).norm.(strcat('norm',num2str(li))) = norm;
%         GaitData.EMG.(side).resample.(strcat('resample',num2str(li))) = resample;
%         li = li +1;
%     end
%     
%     disp(strcat('EMG_',num2str(f),' processed'))
% end
% 
% % Store the muscle names in the GaitData Structure
% GaitData.EMG.R.muscles = EMG_muscles_right;
% GaitData.EMG.L.muscles = EMG_muscles_left;
% 
% % Concatenate the resampled EMG trials
% % Assign an s index variable for slicing muscle names and writing to output
% % structures
% for s = 1:2
%     if s == 1
%         side = 'L';
%     else
%         side = 'R';
%     end
%     
%     % Get a list of the resampled trials
%     fields = fieldnames(GaitData.EMG.(side).resample);
%     
%     % Vertically concatenate the trials
%     for f = 1:length(fields)
%         if f == 1
%             concat = GaitData.EMG.(side).resample.(strcat('resample',num2str(f)));
%         else
%             nextTrial = GaitData.EMG.(side).resample.(strcat('resample',num2str(f)));
%             concat = vertcat(concat,nextTrial);
%         end
%     end
%     GaitData.EMG.(side).concat = concat;
% end

%% Calculated Activations

% Select the OpenSim Model
disp('Select the OpenSim Model')
[modelFile, modelPath] = uigetfile('*.osim','Select the OpenSim model');
model_path = fullfile(modelPath, modelFile);
Misc.model_path = model_path;

% Initialize the left and right indices
li = 1;
ri = 1;

% Loop over the valid gait trials
for f = 1:length(gait_indices)
    % Read the gait side from the meta file
    side = GaitData.meta.side{f};
    % Assign the OutName variable based on the side of the gait
    if side == 'L'
        OutName = strcat('gait_',num2str(li));
    else
        OutName = strcat('gait_',num2str(ri));
    end
    
    % Adjust the Cost Function weights
    Misc.wTres = 10000;
    Misc.wVm= 0.001;

    % Make an output folder
    OutPath = fullfile(dataDir,'Results','Calculated_Activations',side,OutName);
    mkdir(OutPath)
    
    % Set the paths of the output folder, IK file and ID file
    Misc.OutPath = OutPath;
    Misc.IKfile = {fullfile(dataDir,'IK',strcat(GaitData.meta.ID{f},'_MRI_IK.mot'))};
    Misc.IDfile = {fullfile(dataDir,'ID',strcat(GaitData.meta.ID{f},'_MRI_IK_ID.sto'))};

    % Get start and end time of the different files
    % (+50ms beginning and end of time interval, more details see manual and publication)
    IC1 = GaitData.meta.IC1(f);
    IC2 = GaitData.meta.IC2(f);
    time = [IC1 IC2];

    % Settings
    if side == 'L'
        Misc.DofNames_Input={'ankle_angle_l','knee_angle_l','hip_flexion_l','hip_adduction_l','hip_rotation_l'};
    else
        Misc.DofNames_Input={'ankle_angle_r','knee_angle_r','hip_flexion_r','hip_adduction_r','hip_rotation_r'};
    end

    % Plotter Bool: Boolean to select if you want to plot lots of output information of intermediate steps in the script
    Misc.PlotBool = false;

    % MRS Bool: Select if you want to run the generic muscle redundancy solver
    Misc.MRSBool = 1;

    % Validation Bool: Select if you want to run the muscle redundancy solver with the optimized parameters
    Misc.ValidationBool = 0;

    % name output
    Misc.OutName = 'Normal_Gait_';

    % Run muscle tendon estimator:
    [Results,DatStore] = solveMuscleRedundancy(time,Misc);
    
    raw = Results.MActivation.genericMRS'; 

    % Resample
    resample = zeros(101,width(raw));
    
    % Specifiy the input arguments for the interp1 function
    x = (1:1:length(raw))'; % sample points
    xq = linspace(1,length(raw),101); % query point
    
    % Resample each column
    for m = 1:width(raw)
        v = raw(:,m); % sample values
        % resample the kinetic data
        resample(:,m) = interp1(x,v,xq);
    end

    % Store the calculated activations in the GaitData structure
    % Store the processing steps in the GaitData structure    
    if side == 'R'
        GaitData.calc.(side).raw.(strcat('raw',num2str(ri))) = raw;
        GaitData.calc.(side).resample.(strcat('resample',num2str(ri))) = resample;
        if ri == 1
            GaitData.calc.R.muscles = DatStore.MuscleNames;
        end
        ri = ri +1;
    else
        GaitData.calc.(side).raw.(strcat('raw',num2str(li))) = raw;
        GaitData.calc.(side).resample.(strcat('resample',num2str(li))) = resample;
        if li == 1
            GaitData.calc.L.muscles = DatStore.MuscleNames;
        end
        li = li +1;
    end
end

% Concatenate the resampled Calculated Activations
% Assign an s index variable for slicing muscle names and writing to output
% structures
for s = 1:2
    if s == 1
        side = 'L';
    else
        side = 'R';
    end

    % Get a list of the resampled trials
    fields = fieldnames(GaitData.calc.(side).resample);

    % Vertically concatenate the trials
    for f = 1:length(fields)
        if f == 1
            concat = GaitData.calc.(side).resample.(strcat('resample',num2str(f)));
        else
            nextTrial = GaitData.calc.(side).resample.(strcat('resample',num2str(f)));
            concat = vertcat(concat,nextTrial);
        end
    end
    GaitData.calc.(side).concat = concat;
end

%% Get the Reduced Set of Muscle Activations

% Hard Code the names of the calc muscles that correspond to the EMG
% muscles in a side x muscle array
reducedMuscles = ["rect_fem_l","vas_lat_l","bifemlh_l","semiten_l","tib_ant_l","lat_gas_l","soleus_l","glut_med1_l";
    "rect_fem_r","vas_lat_r","bifemlh_r","semiten_r","tib_ant_r","lat_gas_r","soleus_r","glut_med1_r"];

% Loop over the left and right side
% Assign an s index variable for slicing muscle names and writing to output
% structures
for s = 1:2
    if s == 1
        side = 'L';
    else
        side = 'R';
    end

    % Initialize the reduced set array
    calcReduced = zeros(length(GaitData.calc.(side).concat), width(reducedMuscles));

    % Loop over each reduced set muscle
    for m = 1:width(reducedMuscles)
        % Get the full set column that corresponds to the EMG muscle
        full_set_col = find(strcmp(reducedMuscles(s,m),GaitData.calc.(side).muscles));
        % Write each full set column to the reduced set array
        calcReduced(:,m) = GaitData.calc.(side).concat(:,full_set_col);
    end

    % Store the muscle names and reduced set data in the GaitData structure
    GaitData.calcReduced.(side).muscles = reducedMuscles(s,:);
    GaitData.calcReduced.(side).concat = calcReduced;
end

%% Save the GaitData structure
% Create a save path to the results folder
saveName = fullfile(dataDir,'Results','GaitData.mat');
% Save the GaitData Structure
save(saveName,"GaitData","-mat");
