%{
__/\\\\____________/\\\\_________________________________________________/\\\_        
 _\/\\\\\\________/\\\\\\_____________________________________________/\\\\\\\_       
  _\/\\\//\\\____/\\\//\\\____________________________________________\/////\\\_      
   _\/\\\\///\\\/\\\/_\/\\\_____/\\\\\_____/\\\____/\\\_____/\\\\\\\\______\/\\\_     
    _\/\\\__\///\\\/___\/\\\___/\\\///\\\__\//\\\__/\\\____/\\\/////\\\_____\/\\\_    
     _\/\\\____\///_____\/\\\__/\\\__\//\\\__\//\\\/\\\____/\\\\\\\\\\\______\/\\\_   
      _\/\\\_____________\/\\\_\//\\\__/\\\____\//\\\\\____\//\\///////_______\/\\\_  
       _\/\\\_____________\/\\\__\///\\\\\/______\//\\\______\//\\\\\\\\\\_____\/\\\_ 
        _\///______________\///_____\/////_________\///________\//////////______\///_ 

Author: Zach Roth <zachary.roth@student.kuleuven.be>

Synergies1 Workflow:
  ->MovementData1.m
    nmf1.m
    Visualization1.m

Description:
    This script is the first step in the synergies1 workflow. 

Ensure that the data conforms to the following organization and naming
conventions:
    
    Synergies1 (Dir)
        CP* (Dir for each participant)
            EMG (Dir containing EMG data)
                EMG_raw_(trial_ID).xlsx
            ID (Dir containing ID data)
                IK_ID_(trial_ID).sto
            IK
                IK_(trial_ID).mot
            Models
                *.osim
            CP*_trials.xlsx
Inputs:
    - CP*_trials.xlsx: Metadata
    - EMG Data: folder containing a .xlsx file for each movement trial
    - IK Data: folder containing a .mot file for each movement trial
    - ID Data: folder containing a .sto file for each movement trial
    - OpenSim Model: a generic scaled model

Outputs:
    - Results: Dir containing:
        - A folder with the calculated muscle activations from the muscle
        redundancy solver
        - MoveData.mat: a MATLAB structure containing the metadata and all
        intermediate processing steps for the muscle activations
%}

close all; clear; clc

% IMPORT DATA
% Select the synergies1 data folder
disp('Select the synergies1 data folder')
synergiesPath = uigetdir('','Select the synergies1');
cd(synergiesPath)

tic
% Delete the old Results Folder
resultsDir = fullfile(synergiesPath,"Results");
if exist(resultsDir,"dir") == 7
    rmdir(resultsDir,'s')
end

mkdir(resultsDir)

% Get a list of subjects
subjects = dir(fullfile(synergiesPath,"CP*"));
subjects = {subjects.name};

% Loop over the subjects
for subj = 1:length(subjects)
    % Get the path to the Subject data folder
    subjDir = fullfile(synergiesPath,subjects{subj});
    cd(subjDir)
    
    subjResults = fullfile(resultsDir,subjects{subj});
    mkdir(subjResults)
    mkdir(fullfile(subjResults,"Calculated_Activations"))

    % Import the meta data
    meta = readtable(fullfile(subjDir,strcat(subjects{subj},'_trials.xlsx')));
    
    % Get a list of indices for all the valid trials (recorded ICs)
    valid_trials = find(~isnan(meta.IC1));
    
    % Write the meta data to the trials structure
    MoveData.(subjects{subj}).meta.trials = meta(valid_trials,1:6);
    
    % Change the movement names to include gait sides (normal_gait_L,
    % normal_gait_R, etc)
    for f = 1:height(MoveData.(subjects{subj}).meta.trials)
        if MoveData.(subjects{subj}).meta.trials.side{f} == "both"
            continue
        elseif MoveData.(subjects{subj}).meta.trials.side{f} == 'L'
            MoveData.(subjects{subj}).meta.trials.movement{f} = strcat(MoveData.(subjects{subj}).meta.trials.movement{f},'_L');
        elseif MoveData.(subjects{subj}).meta.trials.side{f} == 'R'
            MoveData.(subjects{subj}).meta.trials.movement{f} = strcat(MoveData.(subjects{subj}).meta.trials.movement{f},'_R');
        end
    end
    
    % Sort the meta data based on the new movement names
    MoveData.(subjects{subj}).meta.trials = sortrows(MoveData.(subjects{subj}).meta.trials,"movement","ascend");
    
    % Create new trial names to incorporate numbers (cmj_1, cmj_2, etc.)
    i = 1;
    for f = 1:height(MoveData.(subjects{subj}).meta.trials)
        movement = string(MoveData.(subjects{subj}).meta.trials.movement{f});
        if f == 1
            trialName(f,1) = strcat(movement,"_",num2str(i));
            i = i+1;
        else
            old_movement = string(MoveData.(subjects{subj}).meta.trials.movement{f-1});
            if old_movement == movement
                trialName(f,1) = strcat(movement,"_",num2str(i));
                i = i+1;
            else
                i = 1;
                trialName(f,1) = strcat(movement,"_",num2str(i));
                i = i+1;
            end
        end
    end
    
    % Write the new trial names to the meta data
    MoveData.(subjects{subj}).meta.trials.trial = trialName;
    clear trialName
    % Resort the meta data into the original order
    MoveData.(subjects{subj}).meta.trials = sortrows(MoveData.(subjects{subj}).meta.trials,"ID","ascend");
    
    % Reassign the meta variable to the cleaned trial data
    meta = MoveData.(subjects{subj}).meta.trials;

    % FIND THE MAX EMG VALUES FOR NORMALIZATION
    % Get a list of EMG files
    EMG_dir = fullfile(subjDir,'EMG');
    EMG_files = dir(fullfile(EMG_dir, '*.xlsx'));
    EMG_files = {EMG_files.name};

    % Get the full list of EMG muscle names
    EMG_muscles = readtable(fullfile(EMG_dir,EMG_files{1}));
    EMG_muscles = EMG_muscles.Properties.VariableNames;
    MoveData.(subjects{subj}).meta.EMG_muscles = EMG_muscles(2:end);

    % Get the column names and indices for Left and Right muscles
    EMG_muscles_right_idx = find(strncmpi('R',EMG_muscles,1));
    EMG_muscles_right = EMG_muscles(EMG_muscles_right_idx);

    EMG_muscles_left_idx = find(strncmpi('L',EMG_muscles,1));
    EMG_muscles_left = EMG_muscles(EMG_muscles_left_idx);

    % Combine the muscles into a 1x2 array
    %muscles = [{EMG_muscles_left}, {EMG_muscles_right}];
    
    for f = 1:height(meta)
        % Import the raw EMG files
        rawFile = strcat('EMG_raw',meta.ID{f},'.xlsx');
        raw = readmatrix(fullfile(EMG_dir,rawFile));
        
        % Get the indices of the initial contacts (heel strikes)
        IC1 = meta.IC1(f);
        IC1_idx = find(raw(:,1)==IC1);
        IC2 = meta.IC2(f);
        IC2_idx = find(raw(:,1)==IC2);

        if meta.side{f} == 'L'
            raw = raw(:,EMG_muscles_left_idx);
        elseif meta.side{f} == 'R'
            raw = raw(:,EMG_muscles_right_idx);
        elseif meta.side{f} == "both"
            raw = raw(:,2:end);
        end
        
        % Trim the trials from the first to the second IC
        trim = raw(IC1_idx:IC2_idx,:);

        % Set the mean of the trials to zero
        demean = trim - mean(trim,1);

        % Notch filter: Remove the 50Hz noise
        % Create the notch filter
        fs = 1000;
        fo = 50;
        q = 35;
        bw = (fo/(fs/2))/q;
        [b,a] = iircomb(round(fs/fo),bw,'notch');

        % Initialize an output array
        notch = zeros(size(demean));

        % Filter each column
        for musc = 1:width(notch)
            notch(:,musc) = filter(b,a,trim(:,musc));
        end
        
        % Bandpass filter (20-400 Hz)
        % Create the bandpass filter
        order = 4;
        cutoff = [20 400];
        [b,a] = butter(order/2, cutoff/(0.5*fs));

        % Intialize an output array
        bandPass = zeros(size(notch));

        % Filter each column
        for musc = 1:width(bandPass)
            bandPass(:,musc) = filtfilt(b,a,notch(:,musc));
        end
        
        % Rectification
        rect = abs(bandPass);

        % Smoothing (15 Hz low pass filter)
        % Create the smoothing filter
        order = 4;
        cutoff = 15;
        [b, a] = butter(order, cutoff/(0.5*fs));

        % Initialize an output array
        lowPass = zeros(size(rect));

        % Filter each column
        for musc = 1:width(lowPass)
            lowPass(:,musc) = filtfilt(b,a,rect(:,musc));
        end

        % Write the processed EMG data to the MoveData structure
        MoveData.(subjects{subj}).movements.(meta.movement{f}).EMG.raw.(meta.trial{f}) = raw;
        MoveData.(subjects{subj}).movements.(meta.movement{f}).EMG.trim.(meta.trial{f}) = trim;
        MoveData.(subjects{subj}).movements.(meta.movement{f}).EMG.notch.(meta.trial{f}) = notch;
        MoveData.(subjects{subj}).movements.(meta.movement{f}).EMG.bandPass.(meta.trial{f}) = bandPass;
        MoveData.(subjects{subj}).movements.(meta.movement{f}).EMG.rect.(meta.trial{f}) = rect;
        MoveData.(subjects{subj}).movements.(meta.movement{f}).EMG.lowPass.(meta.trial{f}) = lowPass;
        MoveData.(subjects{subj}).meta.maxEMG.(meta.trial{f}) = lowPass;
    end
    
    % MAX EMG
    % Initialize a output array for the max EMG values
    maxEMG_all = zeros(1,width(MoveData.(subjects{subj}).meta.EMG_muscles));
    
    % Find the max EMG values from each trial
    for f = 1:height(meta)
        maxEMG = max(MoveData.(subjects{subj}).meta.maxEMG.(meta.trial{f}));
        if meta.side{f} == 'L'
            maxCols = maxEMG_all(EMG_muscles_left_idx-1);
            maxEMG_all(EMG_muscles_left_idx-1) = max(maxEMG, maxCols);
        elseif meta.side{f} == 'R'
            maxCols = maxEMG_all(EMG_muscles_right_idx-1);
            maxEMG_all(EMG_muscles_right_idx-1) = max(maxEMG, maxCols);
        elseif meta.side{f} == 'both'
            maxCols = maxEMG_all;
            maxEMG_all = max(maxEMG, maxCols);
        end
    end
    MoveData.(subjects{subj}).meta.maxEMG.maxEMG = maxEMG_all;

    % EMG Processing pt2: Normalize, Resample
    
    for f = 1:height(meta)
        % Get the lowPass data
        lowPass = MoveData.(subjects{subj}).movements.(meta.movement{f}).EMG.lowPass.(meta.trial{f});
        % Select the correct EMG max values
        if meta.side{f} == 'L'
            emgMax = MoveData.(subjects{subj}).meta.maxEMG.maxEMG(EMG_muscles_left_idx-1);
        elseif meta.side{f} == 'R'
            emgMax = MoveData.(subjects{subj}).meta.maxEMG.maxEMG(EMG_muscles_right_idx-1);
        elseif meta.side{f} == "both"
            emgMax = MoveData.(subjects{subj}).meta.maxEMG.maxEMG;
        end
        
        % Normalize
        norm = lowPass./emgMax;

        % Resample
        % Initialize an output array
        resample = zeros(101,width(norm));

        % Specifiy the input arguments for the interp1 function
        x = (1:1:length(norm))'; % sample points
        xq = linspace(1,length(norm),101); % query point

        % Resample each column
        for musc = 1:width(norm)
            v = norm(:,musc); % sample values
            % resample the kinetic data
            resample(:,musc) = interp1(x,v,xq);
        end
        
        % Write the processed EMG data to the MoveData structure
        MoveData.(subjects{subj}).movements.(meta.movement{f}).EMG.norm.(meta.trial{f}) = norm;
        MoveData.(subjects{subj}).movements.(meta.movement{f}).EMG.resample.(meta.trial{f}) = resample;
    end
    
    %% Calculated Activations
    % Select the OpenSim Model
    modelDir = fullfile(subjDir,"Models");
    modelFile = dir(fullfile(modelDir,'*.osim')).name;
    model_path = fullfile(modelDir, modelFile);
    Misc.model_path = model_path;

    for f = 1:height(meta)
        if meta.Calc_trial(f) == 0
            continue
        else
            % Check for the output folder
            MovementOutPath = fullfile(subjResults,"Calculated_Activations",meta.movement{f});
            if exist(MovementOutPath,"dir") == 0
                mkdir(MovementOutPath)
            end

            OutPath = fullfile(subjResults,"Calculated_Activations",meta.movement{f},meta.trial{f});
            mkdir(OutPath)

            Misc.OutPath = OutPath;
            Misc.IKfile = {fullfile(subjDir,'IK',strcat('IK',meta.ID{f},'.mot'))};
            Misc.IDfile = {fullfile(subjDir,'ID',strcat('ID_IK',meta.ID{f},'.sto'))};

            % Get start and end time of the different files
            % (+50ms beginning and end of time interval, more details see manual and publication)
            IC1 = meta.IC1(f);
            IC2 = meta.IC2(f);
            time = [IC1 IC2];

            if meta.side{f} == 'L'
                Misc.DofNames_Input={'ankle_angle_l','knee_angle_l','hip_flexion_l','hip_adduction_l','hip_rotation_l'};
            elseif meta.side{f} == 'R'
                Misc.DofNames_Input={'ankle_angle_r','knee_angle_r','hip_flexion_r','hip_adduction_r','hip_rotation_r'};
            elseif meta.side{f} == "both"
                Misc.DofNames_Input={'ankle_angle_l','knee_angle_l','hip_flexion_l','hip_adduction_l','hip_rotation_l','ankle_angle_r','knee_angle_r','hip_flexion_r','hip_adduction_r','hip_rotation_r'};
            end

            % Plotter Bool: Boolean to select if you want to plot lots of output information of intermediate steps in the script
            Misc.PlotBool = false;

            % MRS Bool: Select if you want to run the generic muscle redundancy solver
            Misc.MRSBool = 1;

            % Validation Bool: Select if you want to run the muscle redundancy solver with the optimized parameters
            Misc.ValidationBool = 0;

            % name output
            Misc.OutName = meta.trial{f};

            % Run muscle tendon estimator:
            [Results,DatStore] = solveMuscleRedundancy(time,Misc);

            MRS = Results.MActivation.genericMRS';

            % Resample
            resample = zeros(101,width(MRS));

            % Specifiy the input arguments for the interp1 function
            x = (1:1:length(MRS))'; % sample points
            xq = linspace(1,length(MRS),101); % query point

            % Resample each column
            for musc = 1:width(MRS)
                v = MRS(:,musc); % sample values
                % resample the MRS data
                resample(:,musc) = interp1(x,v,xq);
            end

            % Store the MRS muscle activations and the resampled data
            MoveData.(subjects{subj}).movements.(meta.movement{f}).calc.MRS.(meta.trial{f}) = MRS;
            MoveData.(subjects{subj}).movements.(meta.movement{f}).calc.resample.(meta.trial{f}) = resample;
            MoveData.(subjects{subj}).meta.calc_Muscles.((meta.movement{f})) = Results.MuscleNames;
        end
    end
 end

%% Save the GaitData structure
% Create a save path to the results folder
saveName = fullfile(synergiesPath,'Results','MoveData.mat');
% Save the GaitData Structure
save(saveName,"MoveData","-mat");

toc