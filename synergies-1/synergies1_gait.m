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
valid_IDs = ["_06" "_07" "_08" "_11"];

%% Import EMG Data

% Select the metadata file
disp('Select the meta data file')
[metafile, metapath] = uigetfile('*.xlsx','Select the meta data file');
meta = readtable(fullfile(metapath,metafile));
cd (metapath)

disp('Select the EMG_max file')
[EMGmax_file, EMGmax_path] = uigetfile('*.xlsx','Select the meta data file');
EMG_max_values = readmatrix(fullfile(EMGmax_path,EMGmax_file));

% Select the folder containing the normal gait data
disp('Select the folder containing the normal gait data')
gait_dir = uigetdir('','Select the folder containing the normal gait data');
tic
cd(gait_dir)

% Get a list of EMG files
EMG_dir = fullfile(gait_dir,'EMG');
EMG_files = dir(fullfile(EMG_dir, '*.xlsx'));
EMG_files = {EMG_files.name}';

% Import the EMG files into the gaits structure
i = 1;
for f = 1:length(EMG_files)
    trial_num = sscanf(EMG_files{f},strcat("EMG_raw_","%f",".xlsx"));
    if ismember(trial_num,valid_trials)
        data = readtable(fullfile(gait_dir,'EMG',EMG_files{f}));
        gaits.EMG.names = data.Properties.VariableNames;
        gaits.EMG.EMG_raw.(strcat('gait_',num2str(i))) = table2array(data);
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
gaits.ICs = zeros(length(valid_trials),2);

i = 1;
for f = 1:height(meta)
    if contains(meta.ID{f},valid_IDs)
        gaits.ICs(i,1) = meta.IC1(f);
        gaits.ICs(i,2) = meta.IC2(f);
        i = i+1;
    else
        continue
    end
end

% Trim the EMG files
for f = 1:length(valid_trials)
    IC1 = gaits.ICs(f,1);
    i_first = find(gaits.EMG.EMG_raw.(strcat('gait_',num2str(f)))(:,1) == gaits.ICs(f,1));
    IC2 = gaits.ICs(f,2);
    i_last = find(gaits.EMG.EMG_raw.(strcat('gait_',num2str(f)))(:,1) == gaits.ICs(f,2));
    gaits.EMG.EMG_trim.(strcat('gait_',num2str(f))) = gaits.EMG.EMG_raw.(strcat('gait_',num2str(f)))(i_first:i_last,:);
    disp(['EMG gait_',num2str(f),' trimmed'])
end    

% Notch filter: Remove the 50Hz noise
fs = 1000;
fo = 50;  
q = 35; 
bw = (fo/(fs/2))/q;
[b,a] = iircomb(round(fs/fo),bw,'notch');


for f = 1:length(valid_trials)
    gaits.EMG.EMG_notch.(strcat('gait_',num2str(f))) = zeros(size(gaits.EMG.EMG_trim.(strcat('gait_',num2str(f)))));
    gaits.EMG.EMG_notch.(strcat('gait_',num2str(f)))(:,1) = gaits.EMG.EMG_trim.(strcat('gait_',num2str(f)))(:,1);
    for m = 2:length(gaits.EMG.names)
        EMG_trim = gaits.EMG.EMG_trim.(strcat('gait_',num2str(f)));
        gaits.EMG.EMG_notch.(strcat('gait_',num2str(f)))(:,m) = filter(b,a,EMG_trim(:,m));
    end
    disp(['EMG gait_',num2str(f),' notch filtered'])
end

% Bandpass filter (20-400 Hz)
order = 4; 
cutoff = [20 400];
[b,a] = butter(order/2, cutoff/(0.5*fs));

for f = 1:length(valid_trials)
    gaits.EMG.EMG_bandpass.(strcat('gait_',num2str(f))) = zeros(size(gaits.EMG.EMG_notch.(strcat('gait_',num2str(f)))));
    gaits.EMG.EMG_bandpass.(strcat('gait_',num2str(f)))(:,1) = gaits.EMG.EMG_notch.(strcat('gait_',num2str(f)))(:,1);
    for m = 2:length(gaits.EMG.names)
        EMG_notch = gaits.EMG.EMG_notch.(strcat('gait_',num2str(f)));
        gaits.EMG.EMG_bandpass.(strcat('gait_',num2str(f)))(:,m) = filtfilt(b,a,EMG_notch(:,m));
    end
    disp(['EMG gait_',num2str(f),' band-pass filtered'])
end

% Rectification
for f = 1:length(valid_trials)
    gaits.EMG.EMG_rect.(strcat('gait_',num2str(f))) = abs(gaits.EMG.EMG_bandpass.(strcat('gait_',num2str(f))));
    disp(['EMG gait_',num2str(f),' rectified'])
end

% 15 Hz low pass filter: Smoothing
order = 4;
cutoff = 15;
[b, a] = butter(order, cutoff/(0.5*fs));

for f = 1:length(valid_trials)
    gaits.EMG.EMG_smooth.(strcat('gait_',num2str(f))) = zeros(size(gaits.EMG.EMG_rect.(strcat('gait_',num2str(f)))));
    gaits.EMG.EMG_smooth.(strcat('gait_',num2str(f)))(:,1) = gaits.EMG.EMG_rect.(strcat('gait_',num2str(f)))(:,1);
    for m = 2:length(gaits.EMG.names)
        EMG_rect = gaits.EMG.EMG_rect.(strcat('gait_',num2str(f)));
        gaits.EMG.EMG_smooth.(strcat('gait_',num2str(f)))(:,m) = filtfilt(b,a,EMG_rect(:,m));
    end
    disp(['EMG gait_',num2str(f),' smoothed'])
end

% Normalize
gaits.EMG.EMG_max = EMG_max_values;

for f = 1:length(valid_trials)
    gaits.EMG.EMG_norm.(strcat('gait_',num2str(f))) = zeros(size(gaits.EMG.EMG_smooth.(strcat('gait_',num2str(f)))));
    gaits.EMG.EMG_norm.(strcat('gait_',num2str(f)))(:,1) = gaits.EMG.EMG_smooth.(strcat('gait_',num2str(f)))(:,1);
    for m = 2:length(gaits.EMG.names)
        EMG_smooth = gaits.EMG.EMG_smooth.(strcat('gait_',num2str(f)));
        gaits.EMG.EMG_norm.(strcat('gait_',num2str(f)))(:,m) = EMG_smooth(:,m)/max(gaits.EMG.EMG_max(:,m));
    end
    disp(['EMG gait_',num2str(f),' normalized'])
end

% Resample
for f = 1:length(valid_trials)
    gaits.EMG.EMG_resamp.(strcat('gait_',num2str(f))) = zeros(101,length(gaits.EMG.names));
    gaits.EMG.EMG_resamp.(strcat('gait_',num2str(f)))(:,1) = (0:100)';
    for m = 2:length(gaits.EMG.names)
        EMG_norm = gaits.EMG.EMG_norm.(strcat('gait_',num2str(f)));
        % sample values
        v = EMG_norm(:,m);
        % sample points
        x = (1:1:length(EMG_norm))';
        % query point
        xq = linspace(1,length(EMG_norm),101);
        % resample the kinetic data
        gaits.EMG.EMG_resamp.(strcat('gait_',num2str(f)))(:,m) = interp1(x,v,xq);
    end
    disp(['EMG gait',num2str(f),' resampled'])
end

%% Visualize EMG
X = (0:100)';
muscles = gaits.EMG.names(2:end);

for f = 1:length(valid_trials)
    t = tiledlayout('flow');
    for m = 1:(length(muscles))
        nexttile
        EMG_resamp = gaits.EMG.EMG_resamp.(strcat('gait_',num2str(f)));
        Y = EMG_resamp(:,m+1) * 100;
        plot(X,Y,'LineWidth',1)
        title(muscles{m})
        ylim([0 30])
        grid on
    end
    
    t.Title.String = 'Relative EMG Activation';
    t.Title.FontWeight = 'bold';

    subtitle = strcat('Normal Gait-',num2str(f));
    t.Subtitle.String = subtitle;

    t.XLabel.String = 'Movement Cycle [%]';
    t.XLabel.FontSize = 14;

    t.YLabel.String = 'Muscle Activation [%]';
    t.YLabel.FontSize = 14;

    savefig(fullfile(gait_dir,'Figures',strcat('normal_gait_',num2str(f),'_EMG.fig')))
end

%% NNMF

max_n_synergies = 6;

for f = 1:length(valid_trials)
    A = gaits.EMG.EMG_resamp.(strcat('gait_',num2str(f)))(:,2:9);
    for s = 1:max_n_synergies
        k = s;
        [W,H,D] = nnmf(A,k,"algorithm","mult","replicates",1000,'Options',statset('Display','final','MaxIter',100));
        gaits.EMG.nnmf.(strcat('gait_',num2str(f))).(strcat('k',num2str(s))).W = W;
        gaits.EMG.nnmf.(strcat('gait_',num2str(f))).(strcat('k',num2str(s))).H = H;
        gaits.EMG.nnmf.(strcat('gait_',num2str(f))).(strcat('k',num2str(s))).D = D;
        gaits.EMG.nnmf.(strcat('gait_',num2str(f))).(strcat('k',num2str(s))).recon = W*H;
    end
end

%% VAF
for f = 1:length(valid_trials)
    data = gaits.EMG.EMG_resamp.(strcat('gait_',num2str(f)))(:,2:9);
    dim_data = size(data);
    for k = 1:max_n_synergies
        data_rec = gaits.EMG.nnmf.(strcat('gait_',num2str(f))).(strcat('k',num2str(k))).recon;
        ursqr = zeros(1,dim_data(2));
        VAF = zeros(1,dim_data(2));
        for i = 1:dim_data(2)
            X = [data(:,i) data_rec(:,i)];
            ursqr(i) = sum(prod(X,2))^2 / (sum(data(:,i).^2)*sum(data_rec(:,i).^2));
            VAF(i) = 1 - sum((data(:,i)-data_rec(:,i)).^2)/sum(data(:,i).^2);
        end
        gaits.EMG.nnmf.(strcat('gait_',num2str(f))).ursqr(k,:) = ursqr;
        gaits.EMG.nnmf.(strcat('gait_',num2str(f))).VAF_muscles(k,:) = VAF;
        gaits.EMG.nnmf.(strcat('gait_',num2str(f))).VAF_trial(k,:) = mean(VAF);
    end
end

%% Minimum Synergies
for f = 1:length(valid_trials)
    gaits.EMG.nnmf.(strcat('gait_',num2str(f))).minSynergies = 999;
    for s = 1:max_n_synergies
        if ismember(1,gaits.EMG.nnmf.(strcat('gait_',num2str(f))).VAF_muscles(s,:) < .7) == 0 && gaits.EMG.nnmf.(strcat('gait_',num2str(f))).VAF_trial(k) > .9 && s < gaits.EMG.nnmf.(strcat('gait_',num2str(f))).minSynergies
            gaits.EMG.nnmf.(strcat('gait_',num2str(f))).minSynergies = s;
        else
            continue
        end
    end
    gaits.EMG.nnmf.min_synergies_summary(f) = gaits.EMG.nnmf.(strcat('gait_',num2str(f))).minSynergies;
end

%% Visualize Synergies for individual gaits
rows = 3;
cols = 2;
X1 = (0:100)';
X2 = categorical(gaits.EMG.names(2:9));
X2 = reordercats(X2,gaits.EMG.names(2:9));

for f = 1:length(valid_trials)
    t = tiledlayout(rows,cols);
    for r = 1:(rows)
        nexttile
        Y1 = gaits.EMG.nnmf.(strcat('gait_',num2str(f))).(strcat('k',num2str(rows))).W(:,r);
        hold on
        plot(X1,Y1)
        ylim([0 1])
        ylabel(strcat('Synergy',num2str(r)))

        if r == 1
            title('Activation Patterns')
        end

        grid on
        hold off
        
        nexttile
        hold on
        Y2 = gaits.EMG.nnmf.(strcat('gait_',num2str(f))).(strcat('k',num2str(rows))).H(r,:);
        bar(X2,Y2)
        ylim([0 1])
        
        if r == 1
            title('Activation Muscle Weightings')
        end    
        
        hold off
    end
    
    tTitle = strcat('Synergies: Normal Gait',num2str(f));
    t.Title.String = tTitle;
    t.Title.FontWeight = 'bold';

    savefig(fullfile(gait_dir,'Figures',strcat('normal_gait_',num2str(f),'_Synergies.fig')))
end

%% Visualize Synergies for all gaits

rows = 3;
cols = 4;
X1 = (0:100)';
X2 = categorical(gaits.EMG.names(2:9));
X2 = reordercats(X2,gaits.EMG.names(2:9));

t = tiledlayout(rows,cols);

for r = 1:(rows)    
    nexttile

    hold on
    Y1 = gaits.EMG.nnmf.gait_1.(strcat('k',num2str(rows))).W(:,r);
    Y2 = gaits.EMG.nnmf.gait_3.(strcat('k',num2str(rows))).W(:,r);
    plot(X1,Y1,X1,Y2)
    ylim([0 1])
    ylabel(strcat('Synergy',num2str(r)))

    if r == 1
        title('Activation Patterns (Left Gait)')
    end

    grid on
    hold off
    
    nexttile

    hold on
    Y1 = gaits.EMG.nnmf.gait_2.(strcat('k',num2str(rows))).W(:,r);
    Y2 = gaits.EMG.nnmf.gait_4.(strcat('k',num2str(rows))).W(:,r);
    plot(X1,Y1,X1,Y2)
    ylim([0 1])
    ylabel(strcat('Synergy',num2str(r)))

    if r == 1
        title('Activation Patterns (Right Gait)')
    end

    grid on
    hold off
    
    nexttile
    hold on
    Y3(:,1) = gaits.EMG.nnmf.gait_1.(strcat('k',num2str(rows))).H(r,:)';
    Y3(:,2) = gaits.EMG.nnmf.gait_3.(strcat('k',num2str(rows))).H(r,:)';
    bar(X2,Y3)
    ylim([0 1])

    if r == 1
        title('Activation Muscle Weightings (Left Gait)')
    end
    grid on
    hold off

    nexttile
    hold on
    Y3(:,1) = gaits.EMG.nnmf.gait_2.(strcat('k',num2str(rows))).H(r,:)';
    Y3(:,2) = gaits.EMG.nnmf.gait_4.(strcat('k',num2str(rows))).H(r,:)';
    bar(X2,Y3)
    ylim([0 1])

    if r == 1
        title('Activation Muscle Weightings (Right Gait)')
    end
    grid on
    hold off
end

t.Title.String = 'Synergies: All Gaits';
t.Title.FontWeight = 'bold';

savefig(fullfile(gait_dir,'Figures',strcat('All_Gaits_Synergies.fig')))

%% Import IK Data
% 
% % Get a list of IK files
% IK_dir = fullfile(gait_dir,'IK');
% IK_files = dir(fullfile(IK_dir,'*.mot'));
% IK_files = {IK_files.name}';
% 
% % Import the IK files into the gaits structure
% i = 1;
% for f = 1:length(IK_files)
%     trial_num = sscanf(IK_files{f},strcat("IK_","%f",".mot"));
%     if ismember(trial_num,valid_trials)
%         Data = ReadMotFile(fullfile(gait_dir,'IK',IK_files{f}));
%         gaits.IK.names = Data.names;
%         gaits.IK.IK_full.(strcat('gait_',num2str(i))) = Data.data;
%         disp(['IK ',num2str(f),' imported'])
%         i = i+1;
%     else
%         disp(['IK ',num2str(f),' skipped (Invalid Trial)'])
%         continue
%     end
% end
% 
% % Get a list of ID files
% ID_dir = fullfile(gait_dir,'ID');
% ID_files = dir(fullfile(ID_dir,'*.sto'));
% ID_files = {ID_files.name}';
% 
% % Import the ID files into the gaits structure
% i = 1;
% for f = 1:length(ID_files)
%     trial_num = sscanf(ID_files{f},strcat("ID_IK_","%f",".sto"));
%     if ismember(trial_num,valid_trials)
%         Data = ReadMotFile(fullfile(gait_dir,'ID',ID_files{f}));
%         gaits.ID.names = Data.names;
%         gaits.ID.ID_full.(strcat('gait_',num2str(i))) = Data.data;
%         disp(['ID ',num2str(f),' imported'])
%         i = i+1;
%     else
%         disp(['ID ',num2str(f),' skipped (Invalid Trial)'])
%         continue
%     end
% end
% 
