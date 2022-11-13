
close all
clear
clc

%% Import Data

% Select the GaitData.mat structure from the Results Folder
disp('Select the GaitData.mat structure from the Results Folder')
[dataFile,dataPath] = uigetfile('*.mat');

% Import the GaitData structure 
load(fullfile(dataPath,dataFile));
cd(dataPath)

tic % Start Stopwatch Timer

% Get a list of the data types in the GaitData structure
% Find the indices of all data types that are not 'meta'
dataType_idx = find(~strcmp(fieldnames(GaitData),'meta'));
% Get a list of GaitData fieldnames
dataTypes = fieldnames(GaitData);
% Filter out the 'meta' fieldname
dataTypes = dataTypes(dataType_idx);

% Hard code the max number of synergies (= num EMG muscles - 1)
max_n_synergies = 7;

%% nnmf & VAF

for d = 1:length(dataTypes) % Loop over each data type
    for s = 1:2 % Loop over each side
        if s == 1
            side = 'L';
        else
            side = 'R';
        end
        
        A = GaitData.(dataTypes{d}).(side).concat; % A = the matrix to factorize

        for k = 1:max_n_synergies % k = Rank of factors (number of synergies)
            [W,H,D] = nnmf(A,k,"algorithm","mult","replicates",8192,'Options',statset('Display','final','MaxIter',256,'UseParallel',true));

            VAF = zeros(1,width(A)); % Initialize a VAF output array
            recon = W*H; % Reconstruct the data by multiplying the synergy activations patterns * muscle weightings
            for i = 1:width(A) % Calculate the VAF for each muscle
                X = [A(:,i) recon(:,i)];
                VAF(i) = 1 - sum((A(:,i)-recon(:,i)).^2)/sum(A(:,i).^2);
            end
            
            % Write the outputs of each synergy analysis to the NNMF struc
            NNMF.(dataTypes{d}).(side).(strcat('k',num2str(k))).W = W;
            NNMF.(dataTypes{d}).(side).(strcat('k',num2str(k))).H = H;
            NNMF.(dataTypes{d}).(side).(strcat('k',num2str(k))).D = D;
            NNMF.(dataTypes{d}).(side).(strcat('k',num2str(k))).recon = recon;
            NNMF.(dataTypes{d}).(side).VAF_muscles(k,:) = VAF;
            NNMF.(dataTypes{d}).(side).VAF_trial(k,:) = mean(VAF);
        end
    end
end

%% Minimum Synergies
for d = 1:length(dataTypes) % Loop over each data type
    for s = 1:2 % Loop over each side
        if s == 1
            side = 'L';
        else
            side = 'R';
        end
        % Set the min_synergies summary variable to an arbitrarily high
        % number
        NNMF.(dataTypes{d}).(side).minSynergies = 999;
        % Get the VAF data
        VAF_muscles = NNMF.(dataTypes{d}).(side).VAF_muscles;

        % Check if n synergies meets the 3 critieria listed below. If it
        % does, it become the new min_synergies
        for n = 1:max_n_synergies
            % Each muscle synergy accounts for > 70% VAF
            muscle_check = ismember(1,VAF_muscles(n,:) < .7) == 0;
            % The average VAF for all trials is > 90%
            trial_check = mean(VAF_muscles(n,:)) > .9;
            % The number of synergies is less than the previous trial that
            % satisfied the previous two conditions
            minSynergies_check = n < NNMF.(dataTypes{d}).(side).minSynergies;

            if muscle_check && trial_check && minSynergies_check
                NNMF.(dataTypes{d}).(side).minSynergies = n;
            else
                continue
            end
        end
        % Write the minSynergies to a Summary Table for convenience
        NNMF.minSynergies_summary(d,s) = NNMF.(dataTypes{d}).(side).minSynergies;
    end
end

%% Save the NNMF structure
% Create a save path to the results folder
saveName = fullfile(dataPath,'NNMF.mat');
% Save the GaitData Structure
save(saveName,"NNMF","-mat");

toc % Stop Stopwatch Timer