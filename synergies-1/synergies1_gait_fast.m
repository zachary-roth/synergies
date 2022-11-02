%{
__/\\\\\\\\\\\\\\\___________________________________________        
 _\/\\\///////////____________________________________________       
  _\/\\\_____________________________________________/\\\______      
   _\/\\\\\\\\\\\______/\\\\\\\\\_____/\\\\\\\\\\__/\\\\\\\\\\\_     
    _\/\\\///////______\////////\\\___\/\\\//////__\////\\\////__    
     _\/\\\_______________/\\\\\\\\\\__\/\\\\\\\\\\____\/\\\______   
      _\/\\\______________/\\\/////\\\__\////////\\\____\/\\\_/\\__  
       _\/\\\_____________\//\\\\\\\\/\\__/\\\\\\\\\\____\//\\\\\___ 
        _\///_______________\////////\//__\//////////______\/////____

Author: Zach Roth <zachary.roth@student.kuleuven.be>
Created: 2-November-2022
%}

close all
clear
clc

% Hard Code Valid Gait Trials
valid_trials = [12,13];
valid_IDs = ["_12" "_13"];

%% Import EMG Data

% Select the metadata file
[metafile, metapath] = uigetfile('*.xlsx','Select the meta data file')
meta = readtable(fullfile(metapath,metafile));

% Select the folder containing the fast gait data
fast_dir = uigetdir('','Select the folder containing the fast gait data')
tic
cd(fast_dir)

% Get a list of EMG files
EMG_dir = fullfile(fast_dir,'EMG');
EMG_files = dir(fullfile(EMG_dir, '*.xlsx'));
EMG_files = {EMG_files.name}';

% Import the EMG files into the fast structure
i = 1;
for f = 1:length(EMG_files)
    trial_num = sscanf(EMG_files{f},strcat("EMG_raw_","%f",".xlsx"));
    if ismember(trial_num,valid_trials)
        data = readtable(fullfile(fast_dir,'EMG',EMG_files{f}));
        fast.EMG.names = data.Properties.VariableNames;
        fast.EMG.EMG_raw.(strcat('gait_',num2str(i))) = table2array(data);
        disp(['EMG ',num2str(f),' imported'])
        i = i+1;
    else
        disp(['EMG ',num2str(f),' skipped (Invalid Trial)'])
        continue
    end
end

%% Process EMG Data

% Get an array of the heel strikes

% Initialize the array
fast.ICs = zeros(length(valid_trials),2);

i = 1;
for f = 1:height(meta)
    if contains(meta.ID{f},valid_IDs)
        fast.ICs(i,1) = meta.IC1(f);
        fast.ICs(i,2) = meta.IC2(f);
        i = i+1;
    else
        continue
    end
end

% Trim the EMG files
for f = 1:length(valid_trials)
    IC1 = fast.ICs(f,1);
    i_first = find(fast.EMG.EMG_raw.(strcat('gait_',num2str(f)))(:,1) == fast.ICs(f,1));
    IC2 = fast.ICs(f,2);
    i_last = find(fast.EMG.EMG_raw.(strcat('gait_',num2str(f)))(:,1) == fast.ICs(f,2));
    fast.EMG.EMG_trim.(strcat('gait_',num2str(f))) = fast.EMG.EMG_raw.(strcat('gait_',num2str(f)))(i_first:i_last,:);
    disp(['EMG gait_',num2str(f),' trimmed'])
end    

% Notch filter: Remove the 50Hz noise
fs = 1000;
fo = 50;  
q = 35; 
bw = (fo/(fs/2))/q;
[b,a] = iircomb(round(fs/fo),bw,'notch');


for f = 1:length(valid_trials)
    fast.EMG.EMG_notch.(strcat('gait_',num2str(f))) = zeros(size(fast.EMG.EMG_trim.(strcat('gait_',num2str(f)))));
    fast.EMG.EMG_notch.(strcat('gait_',num2str(f)))(:,1) = fast.EMG.EMG_trim.(strcat('gait_',num2str(f)))(:,1);
    for m = 2:length(fast.EMG.names)
        EMG_trim = fast.EMG.EMG_trim.(strcat('gait_',num2str(f)));
        fast.EMG.EMG_notch.(strcat('gait_',num2str(f)))(:,m) = filter(b,a,EMG_trim(:,m));
    end
    disp(['EMG gait_',num2str(f),' notch filtered'])
end

% Bandpass filter (20-400 Hz)
order = 4; 
cutoff = [20 400];
[b,a] = butter(order/2, cutoff/(0.5*fs));

for f = 1:length(valid_trials)
    fast.EMG.EMG_bandpass.(strcat('gait_',num2str(f))) = zeros(size(fast.EMG.EMG_notch.(strcat('gait_',num2str(f)))));
    fast.EMG.EMG_bandpass.(strcat('gait_',num2str(f)))(:,1) = fast.EMG.EMG_notch.(strcat('gait_',num2str(f)))(:,1);
    for m = 2:length(fast.EMG.names)
        EMG_notch = fast.EMG.EMG_notch.(strcat('gait_',num2str(f)));
        fast.EMG.EMG_bandpass.(strcat('gait_',num2str(f)))(:,m) = filtfilt(b,a,EMG_notch(:,m));
    end
    disp(['EMG gait_',num2str(f),' band-pass filtered'])
end

% Rectification
for f = 1:length(valid_trials)
    fast.EMG.EMG_rect.(strcat('gait_',num2str(f))) = abs(fast.EMG.EMG_bandpass.(strcat('gait_',num2str(f))));
    disp(['EMG gait_',num2str(f),' rectified'])
end

% 15 Hz low pass filter: Smoothing
order = 4;
cutoff = 15;
[b, a] = butter(order, cutoff/(0.5*fs));

for f = 1:length(valid_trials)
    fast.EMG.EMG_smooth.(strcat('gait_',num2str(f))) = zeros(size(fast.EMG.EMG_rect.(strcat('gait_',num2str(f)))));
    fast.EMG.EMG_smooth.(strcat('gait_',num2str(f)))(:,1) = fast.EMG.EMG_rect.(strcat('gait_',num2str(f)))(:,1);
    for m = 2:length(fast.EMG.names)
        EMG_rect = fast.EMG.EMG_rect.(strcat('gait_',num2str(f)));
        fast.EMG.EMG_smooth.(strcat('gait_',num2str(f)))(:,m) = filtfilt(b,a,EMG_rect(:,m));
    end
    disp(['EMG gait_',num2str(f),' smoothed'])
end

% Normalize
fast.EMG.EMG_max = zeros(length(valid_trials),length(fast.EMG.names));
for f = 1:length(valid_trials)
    EMG_smooth = fast.EMG.EMG_smooth.(strcat('gait_',num2str(f)));
    for m = 2:length(fast.EMG.names)
        fast.EMG.EMG_max(f,m) = max(EMG_smooth(:,m));
    end
end

for f = 1:length(valid_trials)
    fast.EMG.EMG_norm.(strcat('gait_',num2str(f))) = zeros(size(fast.EMG.EMG_smooth.(strcat('gait_',num2str(f)))));
    fast.EMG.EMG_norm.(strcat('gait_',num2str(f)))(:,1) = fast.EMG.EMG_smooth.(strcat('gait_',num2str(f)))(:,1);
    for m = 2:length(fast.EMG.names)
        EMG_smooth = fast.EMG.EMG_smooth.(strcat('gait_',num2str(f)));
        fast.EMG.EMG_norm.(strcat('gait_',num2str(f)))(:,m) = EMG_smooth(:,m)/max(fast.EMG.EMG_max(:,m));
    end
    disp(['EMG gait_',num2str(f),' normalized'])
end

% Resample
for f = 1:length(valid_trials)
    fast.EMG.EMG_resamp.(strcat('gait_',num2str(f))) = zeros(101,length(fast.EMG.names));
    fast.EMG.EMG_resamp.(strcat('gait_',num2str(f)))(:,1) = (0:100)';
    for m = 2:length(fast.EMG.names)
        EMG_norm = fast.EMG.EMG_norm.(strcat('gait_',num2str(f)));
        % sample values
        v = EMG_norm(:,m);
        % sample points
        x = (1:1:length(EMG_norm))';
        % query point
        xq = (1:(length(EMG_norm)/101):length(EMG_norm))';
        % resample the kinetic data
        fast.EMG.EMG_resamp.(strcat('gait_',num2str(f)))(:,m) = interp1(x,v,xq);
    end
    disp(['EMG gait',num2str(f),' resampled'])
end

%% Visualize EMG
X = (0:100)';
muscles = fast.EMG.names(2:end);

for f = 1:length(valid_trials)
    t = tiledlayout('flow');
    for m = 1:(length(muscles))
        nexttile
        EMG_resamp = fast.EMG.EMG_resamp.(strcat('gait_',num2str(f)));
        Y = EMG_resamp(:,m+1) * 100;
        plot(X,Y,'LineWidth',1)
        title(muscles{m})
        ylim([0 105])
        grid on
    end
    
    t.Title.String = 'Relative EMG Activation';
    t.Title.FontWeight = 'bold';

    subtitle = strcat('Fast Gait-',num2str(f));
    t.Subtitle.String = subtitle;

    t.XLabel.String = 'Movement Cycle [%]';
    t.XLabel.FontSize = 14;

    t.YLabel.String = 'Muscle Activation [%]';
    t.YLabel.FontSize = 14;

    savefig(fullfile(fast_dir,'Figures',strcat('fast_gait_',num2str(f),'_EMG.fig')))
end

%% NNMF

max_n_synergies = 6;

for f = 1:length(valid_trials)
    A = fast.EMG.EMG_resamp.(strcat('gait_',num2str(f)))(:,2:9);
    for s = 1:max_n_synergies
        k = s;
        [W,H,D] = nnmf(A,k,"algorithm","als","replicates",10,'Options',statset('Display','final','MaxIter',50))
        fast.EMG.nnmf.(strcat('gait_',num2str(f))).(strcat('k',num2str(s))).W = W;
        fast.EMG.nnmf.(strcat('gait_',num2str(f))).(strcat('k',num2str(s))).H = H;
        fast.EMG.nnmf.(strcat('gait_',num2str(f))).(strcat('k',num2str(s))).D = D;
        fast.EMG.nnmf.(strcat('gait_',num2str(f))).(strcat('k',num2str(s))).recon = W*H;
    end
end

%% VAF
for f = 1:length(valid_trials)
    data = fast.EMG.EMG_resamp.(strcat('gait_',num2str(f)))(:,2:9);
    dim_data = size(data);
    for k = 1:max_n_synergies
        data_rec = fast.EMG.nnmf.(strcat('gait_',num2str(f))).(strcat('k',num2str(k))).recon;
        ursqr = zeros(1,dim_data(2));
        VAF = zeros(1,dim_data(2));
        for i = 1:dim_data(2)
            X = [data(:,i) data_rec(:,i)];
            ursqr(i) = sum(prod(X,2))^2 / (sum(data(:,i).^2)*sum(data_rec(:,i).^2));
            VAF(i) = 1 - sum((data(:,i)-data_rec(:,i)).^2)/sum(data(:,i).^2);
        end
        fast.EMG.nnmf.(strcat('gait_',num2str(f))).ursqr(k,:) = ursqr;
        fast.EMG.nnmf.(strcat('gait_',num2str(f))).VAF(k,:) = VAF;
    end
end

%% Import IK Data

% Get a list of IK files
IK_dir = fullfile(fast_dir,'IK');
IK_files = dir(fullfile(IK_dir,'*.mot'));
IK_files = {IK_files.name}';

% Import the IK files into the gaits structure
i = 1;
for f = 1:length(IK_files)
    trial_num = sscanf(IK_files{f},strcat("IK_","%f",".mot"));
    if ismember(trial_num,valid_trials)
        Data = ReadMotFile(fullfile(fast_dir,'IK',IK_files{f}));
        fast.IK.names = Data.names;
        fast.IK.IK_full.(strcat('gait_',num2str(i))) = Data.data;
        disp(['IK ',num2str(f),' imported'])
        i = i+1;
    else
        disp(['IK ',num2str(f),' skipped (Invalid Trial)'])
        continue
    end
end

% Get a list of ID files
ID_dir = fullfile(fast_dir,'ID');
ID_files = dir(fullfile(ID_dir,'*.sto'));
ID_files = {ID_files.name}';

% Import the ID files into the gaits structure
i = 1;
for f = 1:length(ID_files)
    trial_num = sscanf(ID_files{f},strcat("ID_IK_","%f",".sto"));
    if ismember(trial_num,valid_trials)
        Data = ReadMotFile(fullfile(fast_dir,'ID',ID_files{f}));
        fast.ID.names = Data.names;
        fast.ID.ID_full.(strcat('gait_',num2str(i))) = Data.data;
        disp(['ID ',num2str(f),' imported'])
        i = i+1;
    else
        disp(['ID ',num2str(f),' skipped (Invalid Trial)'])
        continue
    end
end

