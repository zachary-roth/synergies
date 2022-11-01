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

% Hard Code Valid Gait Trials
valid_trials = [6,7,8,11];

%% Import EMG Data

% Select the metadata file
[metafile, metapath] = uigetfile('*.xlsx','Select the meta data file')
meta = readtable(fullfile(metapath,metafile));

% Select the folder containing the normal gait data
gait_dir = uigetdir('','Select the folder containing the normal gait data')
tic
cd(gait_dir)

% Get a list of EMG files
EMG_dir = fullfile(gait_dir,'EMG');
EMG_files = dir(fullfile(EMG_dir, '*.xlsx'));
EMG_files = {EMG_files.name}';

% Import the EMG files into the squats structure
for f = 1:length(EMG_files)
    trial_num = sscanf(EMG_files{f},strcat("EMG_raw_","%f",".xlsx"));
    if ismember(trial_num,valid_trials)
        data = readtable(fullfile(gait_dir,'EMG',EMG_files{f}));
        gaits.EMG.names = data.Properties.VariableNames;
        gaits.EMG.EMG_raw.(strcat('gait_',num2str(f))) = table2array(data);
        disp(['EMG ',num2str(f),' imported'])
    else
        disp(['EMG ',num2str(f),' skipped (Invalid Trial)'])
        continue
    end
end

%% Process EMG Data

% Get an array of the heel strikes

% Initialize the array
gaits.ICs = zeros(length(valid_trials),2);

for f = 1:height(meta)
    if contains(meta.ID{f},string(valid_trials))
        gaits.ICs(f,1) = meta.IC1(f);
        gaits.ICs(f,2) = meta.IC2(f);
    else
        continue
    end
end

% Trim the EMG files



%% NNMF

%% VAF
