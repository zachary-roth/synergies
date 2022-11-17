
close all
clear
clc

disp('Select the NNMF.mat structure')
[file, path] = uigetfile('*.mat','Select the NNMF.mat structure');
load(fullfile(path,file))
cd(path)

load(fullfile(path,'GaitData.mat'))

% Get a list of the data types in the NNMF structure
% Find the indices of all data types that are not 'meta'
dataType_idx = find(~strcmp(fieldnames(NNMF),'minSynergies_summary'));
% Get a list of GaitData fieldnames
dataTypes = fieldnames(NNMF);
% Filter out the 'meta' fieldname
dataTypes = dataTypes(dataType_idx);

% Loop over the left and right side
for s = 1:2
    if s == 1
        side = 'L';
    else
        side = 'R';
    end
    % Loop over each data type
    for d = 1:length(dataTypes)
        min_synergies = NNMF.minSynergies_summary(d,s); % Get the min_synergies from the summary table
        % if no amount of synergies met the selection criteria, continue
        if min_synergies == 999
            continue
        end
        t = tiledlayout(min_synergies,2);
        % Loop over the number of synergies
        for k = 1:min_synergies
            % plot the activation patterns
            nexttile
            yW = NNMF.(dataTypes{d}).(side).(strcat('k',num2str(min_synergies))).W(:,k);
            hold on
            plot(yW)
            xlim([0 100])
            xticks([25 50 75 100])
            xlabel('Gait Cycle [%]')
            title(strcat('Activation Patterns Synergy:'," ",num2str(k)))
            grid on
            hold off

            nexttile
            yH = NNMF.(dataTypes{d}).(side).(strcat('k',num2str(min_synergies))).H(k,:);
            names = GaitData.(dataTypes{d}).(side).muscles;
            names = strrep(names, '_', ' ');
            X = categorical(names);
            X = reordercats(X,names);
            hold on
            bar(X,yH)
            title(strcat('Muscle Weighting Synergy:'," ",num2str(k)))
            grid on
            hold off
        end
        t.Title.String = strcat('Non-Negative Matrix Factorization Results:'," ", dataTypes{d});
        t.Title.FontWeight = 'bold';
        t.Title.FontSize = 16;

        t.Subtitle.String = strcat(side," ",'Gait,'," ", side," ",'Leg Muscles');
        t.Subtitle.FontSize = 14;

        saveName = strcat('nnmf_results_',dataTypes{d},'_',side);
        filename = fullfile(path,"Figures",saveName);
        savefig(filename)
    end
end