
%Author: Zach Roth <zachary.roth@student.kuleuven.be>
%Created: 31-October-2022

clear
close all
clc

% 'Select the folder containing ALL of the EMG trials'
disp('Select the folder containing all of the EMG trials')
EMG_dir = uigetdir('','Select the folder containing all of the EMG trials');
cd(EMG_dir)

% Get a list of all the EMG files
EMG_files = dir(fullfile(EMG_dir,'*.xlsx'));
EMG_files = {EMG_files.name}';

% Import the EMG files into a structure
for f = 1:length(EMG_files)
    EMG_max.(strcat('trial_',num2str(f))) = readmatrix(EMG_files{f});
    disp(['Trial ',num2str(f),' imported'])
end

% Get a list of column names
names = readtable(EMG_files{1});
names = names.Properties.VariableNames;

% Initialize the max array 
EMG_max.maxValues = zeros(1,width(EMG_max.trial_1));

% Find the max EMG values
for f = 1:length(EMG_files)
    for m = 2:width(EMG_max.maxValues)
        maxValue = max(EMG_max.(strcat('trial_',num2str(f)))(:,m));
        if maxValue > EMG_max.maxValues(1,m)
            EMG_max.maxValues(1,m) = maxValue;
        else
            continue
        end
    end
end

MaxValues = array2table(EMG_max.maxValues);
MaxValues.Properties.VariableNames = names;

disp('Select the folder where you want to save the max values')
savepath = uigetdir('','Select the folder where you want to save the max values');
filename = fullfile(savepath,'EMG_max.xlsx');
writetable(MaxValues,filename)
