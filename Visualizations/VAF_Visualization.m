
close all
clear
clc

disp('Select the NNMF.mat structure')
[file, path] = uigetfile('*.mat','Select the NNMF.mat structure');
load(fullfile(path,file))
cd(path)

% Get a list of the data types in the NNMF structure
% Find the indices of all data types that are not 'meta'
dataType_idx = find(~strcmp(fieldnames(NNMF),'minSynergies_summary'));
% Get a list of GaitData fieldnames
dataTypes = fieldnames(NNMF);
% Filter out the 'meta' fieldname
dataTypes = dataTypes(dataType_idx);

left = tiledlayout(length(dataTypes),10);

for d = 1:length(dataTypes)
    % Total VAF
    nexttile([1 3])
    for K = 1:7
        hold on
        Y = NNMF.(dataTypes{d}).L.VAF_trial(K);
        b = bar(K,Y);
        if Y < .9
            b.FaceColor = [0.6350 0.0780 0.1840]; % Red
        else
            b.FaceColor = [0 0.4470 0.7410];
        end
    end

    yline(.9,'--','LineWidth',3)
    xticks(1:7)
    xlabel('n Synergies')
    ylabel('VAF [%]')
    title(strcat('Total VAF:',dataTypes{d}))
    grid on
    hold off

    % Muscle Synergies 
    
    for K = 1:7
        nexttile
        VAF_muscles = NNMF.(dataTypes{d}).L.VAF_muscles(K,:);
        for m = 1:width(VAF_muscles)
            Y = VAF_muscles(m);
            hold on
            b = bar(m,Y);
            if Y < .7
                b.FaceColor = [0.6350 0.0780 0.1840]; % Red
            else
                b.FaceColor = [0 0.4470 0.7410]; % Blue
            end
        end
        yline(.7,'--','LineWidth',3)
        title(strcat('k = ',num2str(K)))
        grid on
        hold off
    end
end

left.Title.String = 'Variability Accounted For in Trials (Left) and Muscles (Right)';
left.Title.FontWeight = 'bold';
left.Title.FontSize = 16;

left.Subtitle.String = strcat('Left Gait, Left Leg Muscles');
left.Subtitle.FontSize = 14;

filename = fullfile(path,"Figures","VAF_Left_Side");
savefig(filename)

right = tiledlayout(length(dataTypes),10);

for d = 1:length(dataTypes)
    % Total VAF
    nexttile([1 3])
    for K = 1:7
        hold on
        Y = NNMF.(dataTypes{d}).R.VAF_trial(K);
        b = bar(K,Y);
        if Y < .9
            b.FaceColor = [0.6350 0.0780 0.1840]; % Red
        else
            b.FaceColor = [0 0.4470 0.7410];
        end
    end

    yline(.9,'--','LineWidth',3)
    xticks(1:7)
    xlabel('n Synergies')
    ylabel('VAF [%]')
    title(strcat('Total VAF:',dataTypes{d}))
    grid on
    hold off

    % Muscle Synergies 
    
    for K = 1:7
        nexttile
        VAF_muscles = NNMF.(dataTypes{d}).R.VAF_muscles(K,:);
        for m = 1:width(VAF_muscles)
            Y = VAF_muscles(m);
            hold on
            b = bar(m,Y);
            if Y < .7
                b.FaceColor = [0.6350 0.0780 0.1840]; % Red
            else
                b.FaceColor = [0 0.4470 0.7410]; % Blue
            end
        end
        yline(.75,'--','LineWidth',3)
        ylim([0 1])
        title(strcat('k = ',num2str(K)))
        grid on
        hold off
    end
end

right.Title.String = 'Variability Accounted For in Trials (Left) and Muscles (Right)';
right.Title.FontWeight = 'bold';
right.Title.FontSize = 16;

right.Subtitle.String = strcat('Right Gait, Right Leg Muscles');
right.Subtitle.FontSize = 14;

filename = fullfile(path,"Figures","VAF_Right_Side");
savefig(filename)