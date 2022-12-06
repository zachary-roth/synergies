%{
__/\\\\\_____/\\\__/\\\\____________/\\\\__/\\\\\\\\\\\\\\\____/\\\\\\\\\_____        
 _\/\\\\\\___\/\\\_\/\\\\\\________/\\\\\\_\/\\\///////////___/\\\///////\\\___       
  _\/\\\/\\\__\/\\\_\/\\\//\\\____/\\\//\\\_\/\\\_____________\///______\//\\\__      
   _\/\\\//\\\_\/\\\_\/\\\\///\\\/\\\/_\/\\\_\/\\\\\\\\\\\_______________/\\\/___     
    _\/\\\\//\\\\/\\\_\/\\\__\///\\\/___\/\\\_\/\\\///////_____________/\\\//_____    
     _\/\\\_\//\\\/\\\_\/\\\____\///_____\/\\\_\/\\\_________________/\\\//________   
      _\/\\\__\//\\\\\\_\/\\\_____________\/\\\_\/\\\_______________/\\\/___________  
       _\/\\\___\//\\\\\_\/\\\_____________\/\\\_\/\\\______________/\\\\\\\\\\\\\\\_ 
        _\///_____\/////__\///______________\///__\///______________\///////////////__

Author: Zach Roth <zachary.roth@student.kuleuven.be>

Synergies2 Workflow:
    digitize metadata (by hand)
    extract_EMG.m
    get_ICs.m
    MovementData2.m
  ->nmf2.m
    Visualization2.m

Description:
    This script is the fourth step in the synergies2 workflow. It is
    similar to nmf1, but only processes the EMG data.

Ensure that the data conforms to the following organization and naming
conventions:
    
    Synergies2 (Dir)
        Results (Dir created by MovementData2.m)
            MoveData.mat

Inputs:
    - MoveData.mat: a MATLAB structure containing the metadata and all
    intermediate processing steps for the muscle activations created by
    MovementData2.m (EMG only).

Outputs:
    - NMF.mat: a MATLAB structure containing the metadata (muscle names), 
    synergies (k=1-6), VAF for each synergy, and a k0 summary table (EMG 
    only).
%}

close all; clear; clc

%% IMPORT DATA
% Select the synergies 1,2,3 or 4 data folder
disp('Select the synergies 1,2,3 or 4 data folder')
synergiesPath = uigetdir('','Select the synergies 1,2,3 or 4 folder');
cd(synergiesPath)

tic % Start the stopwatch timer

% Import the MoveData.mat structure
load(fullfile(synergiesPath,"Results","MoveData.mat"))

% Get a list of subjects
subjects = fieldnames(MoveData);

%% Reformat Movement Data

for subj = 1:length(subjects)

    % Get a list of EMG muscles
    EMG_muscles = MoveData.(subjects{subj}).meta.EMG_muscles;
    % Get the names and indices of the LEFT EMG muscles
    EMG_muscles_left_idx = find(strncmpi('L',EMG_muscles,1));
    EMG_muscles_left = EMG_muscles(EMG_muscles_left_idx);
    % Get the names and indices of the RIGHT EMG muscles
    EMG_muscles_right_idx = find(strncmpi('R',EMG_muscles,1));
    EMG_muscles_right = EMG_muscles(EMG_muscles_right_idx);

    movements = fieldnames(MoveData.(subjects{subj}).movements);
    for move = 1:length(movements)
        % Skip Static, toestanding and unipod movements
        if movements{move} == "static" || movements{move} == "toestanding" || movements{move} == "unipod"
            continue
            % Concatenate the unilateral movements
        else
            % Get a list of all the trials
            trials = fieldnames(MoveData.(subjects{subj}).movements.(movements{move}).EMG.resample);
            concatEMG = [];
            %             concatCalc = [];
            for t = 1:length(trials)
                % Concatenate the trials
                concatEMG = vertcat(concatEMG, MoveData.(subjects{subj}).movements.(movements{move}).EMG.resample.(trials{t}));
                % Store the Left and Right trials
                if contains(movements{move},"_L")
                    NMF.(subjects{subj}).movements.(movements{move}).data.EMG = concatEMG;
                elseif contains(movements{move},"_R")
                    NMF.(subjects{subj}).movements.(movements{move}).data.EMG = concatEMG;
                else
                    % Create new L/R movement names
                    moveLeft = strcat(movements{move},"_L");
                    moveRight = strcat(movements{move},"_R");

                    % Store the split values in the NMF strucutre
                    NMF.(subjects{subj}).movements.(moveLeft).data.EMG = concatEMG(:,EMG_muscles_left_idx);
                    NMF.(subjects{subj}).movements.(moveRight).data.EMG = concatEMG(:,EMG_muscles_right_idx);
                end
                NMF.(subjects{subj}).meta.muscleNames.EMG_L = EMG_muscles_left;
                NMF.(subjects{subj}).meta.muscleNames.EMG_R = EMG_muscles_right;             
            end
        end
    end
end

%% NMF, VAF
for subj = 1:length(subjects)
    movements = fieldnames(NMF.(subjects{subj}).movements);

    for move = 1:length(movements)
        dataSources = fieldnames(NMF.(subjects{subj}).movements.(movements{move}).data);
        for source = 1:length(dataSources)
            A = NMF.(subjects{subj}).movements.(movements{move}).data.(dataSources{source});
            for k = 1:6
                % Non-negative matrix factorization
                [W,H,D] = nnmf(A,k,"algorithm","als","replicates",100,'Options',statset('Display','final','MaxIter',1000,'UseParallel',true));

                % VAF
                VAF_muscles = zeros(1,width(A)); % Initialize a VAF output array
                WH = W*H; % Reconstruct the data by multiplying the synergy activations patterns * muscle weightings
                for i = 1:width(A) % Calculate the VAF for each muscle
                    X = [A(:,i) WH(:,i)];
                    VAF_muscles(i) = 1 - sum((A(:,i)-WH(:,i)).^2)/sum(A(:,i).^2);
                end

                VAF_overall = mean(VAF_muscles);

                % Store Results
                % NMF
                NMF.(subjects{subj}).movements.(movements{move}).results.(dataSources{source}).(strcat('k',num2str(k))).W = W;
                NMF.(subjects{subj}).movements.(movements{move}).results.(dataSources{source}).(strcat('k',num2str(k))).H = H;
                NMF.(subjects{subj}).movements.(movements{move}).results.(dataSources{source}).(strcat('k',num2str(k))).D = D;
                NMF.(subjects{subj}).movements.(movements{move}).results.(dataSources{source}).(strcat('k',num2str(k))).WH = WH;
                % VAF
                NMF.(subjects{subj}).movements.(movements{move}).summary.(dataSources{source}).VAF_muscles(k,:) = VAF_muscles;
                NMF.(subjects{subj}).movements.(movements{move}).summary.(dataSources{source}).VAF_overall(k,:) = VAF_overall;
            end
        end
    end
end

%% Module Selection
for subj = 1:length(subjects)
    movements = fieldnames(NMF.(subjects{subj}).movements);
    for move = 1:length(movements)
        dataSources = fieldnames(NMF.(subjects{subj}).movements.(movements{move}).summary);
        ic = 1;
        for source = 1:length(dataSources)
            k0_VAF = 999;
            VAF = NMF.(subjects{subj}).movements.(movements{move}).summary.(dataSources{source}).VAF_overall;
            for k = 1:6
                if VAF(k) > 0.90 && k < k0_VAF
                    k0_VAF = k;
                end
            end
            % Store k0 results
            NMF.(subjects{subj}).k0{move,ic+1} = k0_VAF;
            ic = ic + 1;
        end
        % Store the movement
        NMF.(subjects{subj}).k0{move,1} = movements{move};
    end
end

%% Save the NMF structure
% Create a save path to the results folder
saveName = fullfile(synergiesPath,'Results','NMF.mat');
% Save the GaitData Structure
save(saveName,"NMF","-mat");

toc