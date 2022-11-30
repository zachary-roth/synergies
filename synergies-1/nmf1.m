%{

__/\\\\\_____/\\\__/\\\\____________/\\\\__/\\\\\\\\\\\\\\\_        
 _\/\\\\\\___\/\\\_\/\\\\\\________/\\\\\\_\/\\\///////////__       
  _\/\\\/\\\__\/\\\_\/\\\//\\\____/\\\//\\\_\/\\\_____________      
   _\/\\\//\\\_\/\\\_\/\\\\///\\\/\\\/_\/\\\_\/\\\\\\\\\\\_____     
    _\/\\\\//\\\\/\\\_\/\\\__\///\\\/___\/\\\_\/\\\///////______    
     _\/\\\_\//\\\/\\\_\/\\\____\///_____\/\\\_\/\\\_____________   
      _\/\\\__\//\\\\\\_\/\\\_____________\/\\\_\/\\\_____________  
       _\/\\\___\//\\\\\_\/\\\_____________\/\\\_\/\\\_____________ 
        _\///_____\/////__\///______________\///__\///______________

Author: Zach Roth (zachary.roth@student.kuleuven.be)
Created on: Nov 22 2022

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

    % PROVIDE (Hard-coded) a list of the calculated
    reducedMuscles_l = ["rect_fem_l","vas_lat_l","bifemlh_l","semiten_l","tib_ant_l","lat_gas_l","soleus_l","glut_med2_l"];
    reducedMuscles_r =["rect_fem_r","vas_lat_r","bifemlh_r","semiten_r","tib_ant_r","lat_gas_r","soleus_r","glut_med2_r"];

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
            concatCalc = [];
            for t = 1:length(trials)
                % Concatenate the trials
                concatEMG = vertcat(concatEMG, MoveData.(subjects{subj}).movements.(movements{move}).EMG.resample.(trials{t}));
                concatCalc = vertcat(concatCalc, MoveData.(subjects{subj}).movements.(movements{move}).calc.resample.(trials{t}));

                % Get a list of calculated muscle names
                calc_Muscles = MoveData.(subjects{subj}).meta.calc_Muscles.(movements{move});
                % Get the indices for the reduced set
                [~, reduced_Muscles_left_idx] = intersect(calc_Muscles,reducedMuscles_l);
                [~, reduced_Muscles_right_idx] = intersect(calc_Muscles,reducedMuscles_r);

                % Store the Left and Right trials
                if contains(movements{move},"_L")
                    NMF.(subjects{subj}).movements.(movements{move}).data.EMG = concatEMG;
                    NMF.(subjects{subj}).movements.(movements{move}).data.calc = concatCalc;
                    NMF.(subjects{subj}).movements.(movements{move}).data.calcReduced = concatCalc(:,reduced_Muscles_left_idx);
                elseif contains(movements{move},"_R")
                    NMF.(subjects{subj}).movements.(movements{move}).data.EMG = concatEMG;
                    NMF.(subjects{subj}).movements.(movements{move}).data.calc = concatCalc;
                    NMF.(subjects{subj}).movements.(movements{move}).data.calcReduced = concatCalc(:,reduced_Muscles_right_idx);
                else
                    % Get the indices of the Right and Left muscles
                    calc_Muscles_left_idx = find(endsWith(calc_Muscles,'l'));
                    calc_Muscles_left = calc_Muscles(calc_Muscles_left_idx);

                    calc_Muscles_right_idx = find(endsWith(calc_Muscles,'r'));
                    calc_Muscles_right = calc_Muscles(calc_Muscles_right_idx);

                    % Create new L/R movement names
                    moveLeft = strcat(movements{move},"_L");
                    moveRight = strcat(movements{move},"_R");

                    % Store the split values in the NMF strucutre
                    NMF.(subjects{subj}).movements.(moveLeft).data.EMG = concatEMG(:,EMG_muscles_left_idx);
                    NMF.(subjects{subj}).movements.(moveRight).data.EMG = concatEMG(:,EMG_muscles_right_idx);
                    NMF.(subjects{subj}).movements.(moveLeft).data.calc = concatCalc(:,calc_Muscles_left_idx);
                    NMF.(subjects{subj}).movements.(moveRight).data.calc = concatCalc(:,calc_Muscles_right_idx);
                    NMF.(subjects{subj}).movements.(moveLeft).data.calcReduced = concatCalc(:,reduced_Muscles_left_idx);
                    NMF.(subjects{subj}).movements.(moveRight).data.calcReduced = concatCalc(:,reduced_Muscles_right_idx);
                end
                NMF.(subjects{subj}).meta.muscleNames.EMG_L = EMG_muscles_left;
                NMF.(subjects{subj}).meta.muscleNames.EMG_R = EMG_muscles_right;
                NMF.(subjects{subj}).meta.muscleNames.calc_l = calc_Muscles_left;
                NMF.(subjects{subj}).meta.muscleNames.calc_r = calc_Muscles_right;
                NMF.(subjects{subj}).meta.muscleNames.calcReduced_l = reducedMuscles_l;
                NMF.(subjects{subj}).meta.muscleNames.calcReduced_r = reducedMuscles_r;
                
            end
        end
    end
end

%% NMF, VAF, AIC, BIC
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

                % BIC
                %                 m = length(A); % n time points
                %                 n = width(A); % n muscles
                %                 c = min(sqrt(m),sqrt(n));
                %
                %                 BIC1 = log((norm(WH-A,'fro')^2)) + k*((m+n)/m*n) * log(m*n/(m+n));
                %                 BIC2 = log((norm(WH-A,'fro')^2)) + k*((m+n)/m*n) * log(c^2);
                %                 BIC3 = log((norm(WH-A,'fro')^2)) + k*((m+n)/m*n) * (log(c^2))/(c^2);

                % Store Results
                % NMF
                NMF.(subjects{subj}).movements.(movements{move}).results.(dataSources{source}).(strcat('k',num2str(k))).W = W;
                NMF.(subjects{subj}).movements.(movements{move}).results.(dataSources{source}).(strcat('k',num2str(k))).H = H;
                NMF.(subjects{subj}).movements.(movements{move}).results.(dataSources{source}).(strcat('k',num2str(k))).D = D;
                NMF.(subjects{subj}).movements.(movements{move}).results.(dataSources{source}).(strcat('k',num2str(k))).WH = WH;
                % VAF
                NMF.(subjects{subj}).movements.(movements{move}).summary.(dataSources{source}).VAF_muscles(k,:) = VAF_muscles;
                NMF.(subjects{subj}).movements.(movements{move}).summary.(dataSources{source}).VAF_overall(k,:) = VAF_overall;

                % BIC
                %                 NMF.(subjects{subj}).movements.(movements{move}).summary.(dataSources{source}).BIC(k,1) = BIC1;
                %                 NMF.(subjects{subj}).movements.(movements{move}).summary.(dataSources{source}).BIC(k,2) = BIC2;
                %                 NMF.(subjects{subj}).movements.(movements{move}).summary.(dataSources{source}).BIC(k,3) = BIC3;
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