close all; clear; clc;

% Select the synergies Data folder
disp('Select the synergies Data folder')
synergiesPath = uigetdir('','Select the synergies Data folder');
cd(synergiesPath)

dataDir = fullfile(synergiesPath,"Data");
synergies2Dir = fullfile(synergiesPath,"Synergies2");

subjectsAll = dir(fullfile(dataDir,"CP*"));
subjectsAll = {subjectsAll.name};

subjects = dir(fullfile(synergies2Dir,"CP*"));
subjects = {subjects.name};

for subj = 1:length(subjects)
    % Get the index of the Synergies2 subject from the list of all subjects
    iAll = find(strcmp(subjects{subj},subjectsAll));

    % Import the meta data
    meta = readtable(fullfile(synergies2Dir,subjects{subj},strcat(subjects{subj},"_trials")));
    
    files = unique(meta.ID);

    % Create an output folder for the EMG data
    EMG_dir = fullfile(synergies2Dir,subjects{subj},"EMG");
    if exist(EMG_dir,"dir") == 7
        rmdir(EMG_dir,"s")
    end
    mkdir(EMG_dir)

    for f = 1:length(files)
        % Import the c3d file
        c3dName = strcat(subjects{subj},"_T0",files{f},".c3d");
        c3dfile = fullfile(dataDir,subjects{subj},"Vicon",c3dName);
        [Markers,MLabels,VideoFrameRate,AnalogSignals,ALabels,AUnits,AnalogFrameRate,Event,ParameterGroup,CameraInfo]= readC3D(c3dfile, []);

        % Extract EMG Data
        EMG_pat = 'Voltage.' + lettersPattern(4);
        iEMG = find(contains(ALabels,EMG_pat));
        EMG = AnalogSignals(:,iEMG);

        EMG_names = extractAfter(ALabels(iEMG),"Voltage.");
        EMG = array2table(EMG,'VariableNames',EMG_names);
        Time = ((1:length(AnalogSignals))/1000)';
        Time = array2table(Time,'VariableNames',"Time");

        EMG = [Time,EMG];

        filename = fullfile(EMG_dir, strcat("EMG",files{f},".xlsx"));
        writetable(EMG,filename)
        disp(strcat("Created"," ",fullfile(EMG_dir, strcat("EMG",files{f},".xlsx"))))
    end % File Loop
end