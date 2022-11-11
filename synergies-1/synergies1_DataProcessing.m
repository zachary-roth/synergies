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
%}

close all
clear
clc

%%   Select the Subject1 data folder
disp('Select the Subject1 data folder')
dataDir = uigetdir('','Select the Subject1 data folder');
cd(dataDir)

if exist("Results","dir") == 7
    rmdir("Results",'s')
end

% Create the results folder and subfolders
mkdir("Results")
mkdir(fullfile("Results","Calculated_Activations"))
mkdir(fullfile("Results","Calculated_Activations","L"))
mkdir(fullfile("Results","Calculated_Activations","R"))
mkdir(fullfile("Results","Figures"))


%% Meta Data
%Import the meta data
meta = readtable(fullfile(dataDir,'trials.xlsx'));

% Get a list of indices for all the normal gait trials
normal_gaits = find(strcmp(meta.movement,'gait_normal'));

% Get a list of indices for all the valid trials (recorded ICs)
valid_ICs = find(~isnan(meta.IC1));

% Filter out the invalid normal gait trials
gait_indices = intersect(normal_gaits,valid_ICs);

% Write the meta data to the trials structure
GaitData.meta.ID = meta.ID(gait_indices);
GaitData.meta.side = meta.side(gait_indices);
GaitData.meta.IC1 = meta.IC1(gait_indices);
GaitData.meta.IC2 = meta.IC2(gait_indices);

%% EMG
EMG_max_path = fullfile(dataDir,"EMG_max.xlsx");
EMG_max_all = readmatrix(EMG_max_path);

% Get a list of EMG files
EMG_dir = fullfile(dataDir,'EMG');
EMG_files = dir(fullfile(EMG_dir, '*.xlsx'));
EMG_files = {EMG_files.name};


% Import and Process the EMG files into the gaits structure

li = 1;
ri = 1;

for f = 1:length(gait_indices)
    
    side = GaitData.meta.side{f};
    
    % Read in the raw EMG files
    rawFile = strcat('EMG_raw',GaitData.meta.ID{f},'.xlsx');
    raw = readmatrix(fullfile(EMG_dir,rawFile));

    % Get the indices of the initial contacts (heel strikes)
    IC1 = GaitData.meta.IC1(f);
    IC1_idx = find(raw(:,1)==IC1);
    IC2 = GaitData.meta.IC2(f);
    IC2_idx = find(raw(:,1)==IC2);
    
    % Take the right muscles for right gaits and vice-versa
    if GaitData.meta.side{f} == 'R'
        trim = raw(IC1_idx:IC2_idx,2:9);
    elseif GaitData.meta.side{f} == 'L'
        trim = raw(IC1_idx:IC2_idx,10:17);
    end
    
    % Notch filter: Remove the 50Hz noise
    fs = 1000;
    fo = 50;
    q = 35;
    bw = (fo/(fs/2))/q;
    [b,a] = iircomb(round(fs/fo),bw,'notch');
    
    notch = zeros(size(trim));
    
    for m = 1:8
        notch(:,m) = filter(b,a,trim(:,m));
    end

    % Bandpass filter (20-400 Hz)
    order = 4;
    cutoff = [20 400];
    [b,a] = butter(order/2, cutoff/(0.5*fs));

    bandPass = zeros(size(notch));

    for m = 1:8
        bandPass(:,m) = filtfilt(b,a,notch(:,m));
    end

    % Rectification
    rect = abs(bandPass);

    % Smoothing (15 Hz low pass filter)
    order = 4;
    cutoff = 15;
    [b, a] = butter(order, cutoff/(0.5*fs));

    lowPass = zeros(size(rect));

    for m = 1:8
        lowPass(:,m) = filtfilt(b,a,rect(:,m));
    end

    % Normalize
    
    norm = zeros(size(lowPass));

    if GaitData.meta.side{f} == 'R'
        EMG_max = EMG_max_all(:,2:9);
    elseif GaitData.meta.side{f} == 'L'
        EMG_max = EMG_max_all(:,10:17);
    end

    for m = 1:8
        norm(:,m) = lowPass(:,m)/EMG_max(m);
    end

    % Resample 
    
    resample = zeros(101,8);

    x = (1:1:length(norm))'; % sample points
    xq = linspace(1,length(norm),101); % query point

    for m = 1:8
        v = norm(:,m); % sample values
        % resample the kinetic data
        resample(:,m) = interp1(x,v,xq);
    end

    % Store the processing steps in the GaitData structure    
    if side == 'R'
        GaitData.EMG.(side).raw.(strcat('raw',num2str(ri))) = raw;
        GaitData.EMG.(side).trim.(strcat('trim',num2str(ri))) = trim;
        GaitData.EMG.(side).notch.(strcat('notch',num2str(ri))) = notch;
        GaitData.EMG.(side).bandPass.(strcat('bandPass',num2str(ri))) = bandPass;
        GaitData.EMG.(side).rect.(strcat('rect',num2str(ri))) = rect;
        GaitData.EMG.(side).lowPass.(strcat('lowPass',num2str(ri))) = lowPass;
        GaitData.EMG.(side).norm.(strcat('norm',num2str(ri))) = norm;
        GaitData.EMG.(side).resample.(strcat('resample',num2str(ri))) = resample;
        ri = ri +1;
    else
        GaitData.EMG.(side).raw.(strcat('raw',num2str(li))) = raw;
        GaitData.EMG.(side).trim.(strcat('trim',num2str(li))) = trim;
        GaitData.EMG.(side).notch.(strcat('notch',num2str(li))) = notch;
        GaitData.EMG.(side).bandPass.(strcat('bandPass',num2str(li))) = bandPass;
        GaitData.EMG.(side).rect.(strcat('rect',num2str(li))) = rect;
        GaitData.EMG.(side).lowPass.(strcat('lowPass',num2str(li))) = lowPass;
        GaitData.EMG.(side).norm.(strcat('norm',num2str(li))) = norm;
        GaitData.EMG.(side).resample.(strcat('resample',num2str(li))) = resample;
        li = li +1;
    end
    
    disp(strcat('EMG_',num2str(f),' processed'))
end

% Concatenate the resampled EMG trials
for s = 1:2
    if s == 1
        side = 'L';
    else
        side = 'R';
    end

    fields = fieldnames(GaitData.EMG.(side).resample);
    
    for f = 1:length(fields)
        if f == 1
            concat = GaitData.EMG.(side).resample.(strcat('resample',num2str(f)));
        else
            nextTrial = GaitData.EMG.(side).resample.(strcat('resample',num2str(f)));
            concat = vertcat(concat,nextTrial);
        end
    end
    GaitData.EMG.(side).concat = concat;
end

%% Calculated Activations

% HARD CODED
model_path = fullfile(dataDir,'Models','GEN_scaled.osim');
Misc.model_path = model_path;

li = 1;
ri = 1;

for f = 1:length(gait_indices)
    
    side = GaitData.meta.side{f};
    
    if side == 'L'
        OutName = strcat('gait_',num2str(li));
    else
        OutName = strcat('gait_',num2str(ri));
    end
    
    % Make an output folder
    OutPath = fullfile(dataDir,'Results','Calculated_Activations',side,OutName);
    mkdir(OutPath)
    
    % Set the paths of the output folder, IK file and ID file
    Misc.OutPath = OutPath;
    Misc.IKfile = {fullfile(dataDir,'IK',strcat('IK',GaitData.meta.ID{f},'.mot'))};
    Misc.IDfile = {fullfile(dataDir,'ID',strcat('ID_IK',GaitData.meta.ID{f},'.sto'))};

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
    resample = zeros(101,43);

    x = (1:1:length(raw))'; % sample points
    xq = linspace(1,length(raw),101); % query point

    for m = 1:43
        v = raw(:,m); % sample values
        % resample the kinetic data
        resample(:,m) = interp1(x,v,xq);
    end

    % Store the calculated activations in the GaitData structure
    % Store the processing steps in the GaitData structure    
    if side == 'R'
        GaitData.calc.(side).raw.(strcat('raw',num2str(ri))) = raw;
        GaitData.calc.(side).resample.(strcat('resample',num2str(ri))) = resample;
        ri = ri +1;
    else
        GaitData.calc.(side).raw.(strcat('raw',num2str(li))) = raw;
        GaitData.calc.(side).resample.(strcat('resample',num2str(li))) = resample;
        li = li +1;
    end
end

saveName = fullfile(dataDir,'Results','GaitData');
save(saveName,'GaitData','-mat');