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
subjectDir = uigetdir('','Select the Subject1 data folder');
cd(subjectDir)

%% Import the meta data
meta = readtable(fullfile(subjectDir,'trials.xlsx'));

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

EMG_max_path = fullfile(subjectDir,"EMG_max.xlsx");
EMG_max_all = readmatrix(EMG_max_path);

% Get a list of EMG files
EMG_dir = fullfile(subjectDir,'EMG');
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