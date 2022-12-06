%{
__/\\\________/\\\_________________________/\\\_        
 _\/\\\_______\/\\\_____________________/\\\\\\\_       
  _\//\\\______/\\\___/\\\______________\/////\\\_      
   __\//\\\____/\\\___\///___/\\\\\\\\\\_____\/\\\_     
    ___\//\\\__/\\\_____/\\\_\/\\\//////______\/\\\_    
     ____\//\\\/\\\_____\/\\\_\/\\\\\\\\\\_____\/\\\_   
      _____\//\\\\\______\/\\\_\////////\\\_____\/\\\_  
       ______\//\\\_______\/\\\__/\\\\\\\\\\_____\/\\\_ 
        _______\///________\///__\//////////______\///_ 

Author: Zach Roth <zachary.roth@student.kuleuven.be>

Synergies1/Synergies3 Workflow:
    MovementData1.m/MovementData3.m
    nmf1.m
  ->Visualization1.m

Description:
    This script is the final step in the synergies1 and synergies 3
workflows.

Ensure that the data conforms to the following organization and naming
conventions:
    
    Synergies1 or Synergies3(Dir)
            Results (Dir created by MovementData1.m or MovementData3.m)
                MoveData.mat
                NMF.mat

Inputs:
    - MoveData.mat: a MATLAB structure containing the metadata and all
    intermediate processing steps for the muscle activations
    - NMF.mat: a MATLAB structure containing the metadata (muscle names),
    synergies (k=1-6), VAF for each synergy, and a k0 summary table.

Outputs:
    Visualizations: Dir containing: 
        - EMG Validation Plots: Check for spikes in EMG which may indicate
        that the filtering was in sufficient or that outliers need to be
        removed. 
        - A_vs_WH Plots: Compare the observed activations with the
        reconstructed activations (using the appropraite k0). 
        - VAF Plot: Compare the overall VAF and VAF per muscle for each of
        the 6 computed synergies.
        - NMF Plots: Visualize the synergy activation patterns (W) and
        muscle weightings (H), using the appropriate k0.
%}

close all; clear; clc

%% IMPORT DATA
% Select the synergies 1 or 3 data folder
disp('Select the synergies 1 or 3 data folder')
synergiesPath = uigetdir('','Select the synergies 1 or 3 data folder');
cd(synergiesPath)

load(fullfile(synergiesPath,"Results","MoveData.mat"))
load(fullfile(synergiesPath,"Results","NMF.mat"))

% Delete the old Results Folder
visDir = fullfile(synergiesPath,"Results","Visualizations");
if exist(visDir,"dir") == 7
    rmdir(visDir,'s')
end

mkdir(visDir)

subjects = fieldnames(NMF);

for subj = 1:length(subjects)
    subj_visDir = fullfile(visDir,subjects{subj});
    mkdir(subj_visDir)

    %% Plot the low-pass filtered EMG data for validation
    EMG_validation_dir = fullfile(subj_visDir,"EMG_Validation");
    mkdir(EMG_validation_dir)
    
    movements = fieldnames(MoveData.(subjects{subj}).movements);
    for move = 1:length(movements)
        t = tiledlayout('flow');
        trials = fieldnames(MoveData.(subjects{subj}).movements.(movements{move}).EMG.lowPass);
        n_muscles = width(MoveData.(subjects{subj}).movements.(movements{move}).EMG.lowPass.(trials{1}));
        for musc = 1:n_muscles
            nexttile
            for f = 1:length(trials)
                lowPass = MoveData.(subjects{subj}).movements.(movements{move}).EMG.lowPass.(trials{f})(:,musc);
                hold on
                plot(lowPass);
            end % Trials Loop (Plot)
            
            % Get the plot titles and the Upper X Limits
            if contains(movements{move},"_R")
                muscles = MoveData.(subjects{subj}).meta.EMG_muscles(1:8);
                upperY = MoveData.(subjects{subj}).meta.maxEMG.maxEMG(1:8);
            elseif contains(movements{move},"_L")
                muscles = MoveData.(subjects{subj}).meta.EMG_muscles(9:16);
                upperY = MoveData.(subjects{subj}).meta.maxEMG.maxEMG(9:16);
            else
                muscles = MoveData.(subjects{subj}).meta.EMG_muscles;
                upperY = MoveData.(subjects{subj}).meta.maxEMG.maxEMG;
            end
                
            ylim([0 upperY(musc)])
            title(muscles{musc})
            grid on
            hold off
        end % Muscles Loop (Tile)
        plotTitle = strcat(subjects{subj}," ",movements{move});
        plotTitle = strrep(plotTitle,"_"," ");
        plotTitle = upper(plotTitle);
        t.Title.String = plotTitle;
        t.Title.FontSize = 16;
        t.Title.FontWeight = "Bold";
        
        t.XLabel.String = 'Time [ms]';
        t.XLabel.FontSize = 14;

        t.YLabel.String = 'Muscle Activation [mV]';
        t.YLabel.FontSize = 14;

        filename = strcat("EMG_Validation_",subjects{subj},"_",(movements{move}));
        filename = fullfile(EMG_validation_dir,filename);
        savefig(filename)
    end % Movments Loop (Tiled Layout)

    %% Plot the Observed Activations (A) against the Reconstructed Activations (WH)
    A_vs_WH_dir = fullfile(subj_visDir,"A_vs_WH");
    mkdir(A_vs_WH_dir)

    movements = fieldnames(NMF.(subjects{subj}).movements);
    for move = 1:length(movements)
        t = tiledlayout(4,2);
        A_EMG = NMF.(subjects{subj}).movements.(movements{move}).data.EMG;
        A_calcReduced = NMF.(subjects{subj}).movements.(movements{move}).data.calcReduced;
        
        k0_EMG = NMF.(subjects{subj}).k0{move,2};
        k0_calcReduced = NMF.(subjects{subj}).k0{move,4};
        
        WH_EMG = NMF.(subjects{subj}).movements.(movements{move}).results.EMG.(strcat("k",num2str(k0_EMG))).WH;
        WH_calcReduced = NMF.(subjects{subj}).movements.(movements{move}).results.calcReduced.(strcat("k",num2str(k0_calcReduced))).WH;
        
        muscles = NMF.(subjects{subj}).meta.muscleNames.calcReduced_r;
        muscles = strrep(muscles,"_r","");

        for musc = 1:8
            nexttile
            hold on
            plot(A_EMG(:,musc),"Color","#0072BD");
            plot(WH_EMG(:,musc),"Color","#0072BD","LineStyle",":");
            plot(A_calcReduced(:,musc),"Color","#A2142F");
            plot(WH_calcReduced(:,musc),"Color","#A2142F","LineStyle",":");
            xlim([0 100])
            tileTitle = muscles(musc);
            tileTitle = strrep(tileTitle,"_"," ");
            title(tileTitle)
            grid on
            hold off
        end

        plotTitle = strcat(subjects{subj}," ",movements{move});
        plotTitle = strrep(plotTitle,"_"," ");
        plotTitle = upper(plotTitle);
        t.Title.String = plotTitle;
        t.Title.FontSize = 16;
        t.Title.FontWeight = "Bold";

        t.XLabel.String = 'Movement Cycle [%]';
        t.XLabel.FontSize = 14;

        t.YLabel.String = 'Muscle Activation [%]';
        t.YLabel.FontSize = 14;

        leg = legend('A EMG','WH EMG','A CalcReduced','WH CalcReduced','Orientation','Horizontal');
        leg.Layout.Tile = 'north';

        filename = strcat("A_vs_WH_",subjects{subj},"_",(movements{move}));
        filename = fullfile(A_vs_WH_dir,filename);
        savefig(filename)
    end % A_vs_WH Loop

    %% VAF
    VAF_vis_dir = fullfile(subj_visDir,"VAF_Visualization");
    mkdir(VAF_vis_dir)

    movements = fieldnames(NMF.(subjects{subj}).movements);
    for move = 1:length(movements)
        dataTypes = fieldnames(NMF.(subjects{subj}).movements.(movements{move}).results);
        t = tiledlayout(length(dataTypes),7);
        t.TileSpacing = "tight";

        for d = 1:length(dataTypes)

            % Muscle Names
            if contains(movements{move},"_R")
                musclesEMG = NMF.(subjects{subj}).meta.muscleNames.EMG_R;
                musclesCalc = NMF.(subjects{subj}).meta.muscleNames.calc_r;
                musclesCalcReduced = NMF.(subjects{subj}).meta.muscleNames.calcReduced_r;
            elseif contains(movements{move},"_L")
                musclesEMG = NMF.(subjects{subj}).meta.muscleNames.EMG_L;
                musclesCalc = NMF.(subjects{subj}).meta.muscleNames.calc_l;
                musclesCalcReduced = NMF.(subjects{subj}).meta.muscleNames.calcReduced_l;
            end

            % Overall VAF
            nexttile
            for K = 1:6
                hold on
                Y = NMF.(subjects{subj}).movements.(movements{move}).summary.(dataTypes{d}).VAF_overall(K);
                b = bar(K,Y);
                if Y < .9
                    b.FaceColor = [0.6350 0.0780 0.1840]; % Red
                else
                    b.FaceColor = [0 0.4470 0.7410]; % Blue
                end
            end
            yline(.9,'--','LineWidth',3)
            xticks(1:6)
            xlabel('n Synergies')
            title(strcat('Overall VAF:',dataTypes{d}))
            grid on
            hold off
            
            % Muscle Synergies
            for K = 1:6
                nexttile
                VAF_muscles = NMF.(subjects{subj}).movements.(movements{move}).summary.(dataTypes{d}).VAF_muscles(K,:);
                for m = 1:width(VAF_muscles)
                    Y = VAF_muscles(m);
                    hold on
                    b = bar(m,Y);
                    if Y < .75
                        b.FaceColor = [0.6350 0.0780 0.1840]; % Red
                    else
                        b.FaceColor = [0 0.4470 0.7410]; % Blue
                    end
                end

                if dataTypes{d} == "EMG"
                    muscles = (musclesEMG);
                elseif dataTypes{d} == "calc"
                    muscles = (musclesCalc);
                elseif dataTypes{d} == "calcReduced"
                    muscles = (musclesCalcReduced);
                end

                muscles = strrep(muscles,"_"," ");
                xticks(1:length(VAF_muscles))
                xticklabels(muscles)
                
                yline(.75,'--','LineWidth',3)
                title(strcat('VAF k = ',num2str(K)))
                grid on
                hold off
            end
        end % Datatypes Loop
        
        movementTitle = strrep(movements{move},"_"," ");
        t.Title.String = strcat('Variability Accounted For in'," ",movementTitle," ","(",subjects{subj},")");
        t.Title.FontWeight = 'bold';
        t.Title.FontSize = 16;
        
        t.YLabel.String = 'Variability Accounted For [%]';
        t.YLabel.FontSize = 14;

        filename = strcat("VAF_",subjects{subj},"_",(movements{move}));
        filename = fullfile(VAF_vis_dir,filename);
        savefig(filename)
    end % VAF Loop

    %% NMF
    nmf_vis_dir = fullfile(subj_visDir,"NMF_Visualization");
    mkdir(nmf_vis_dir)

    movements = fieldnames(NMF.(subjects{subj}).movements);
    for move = 1:length(movements)
        dataTypes = fieldnames(NMF.(subjects{subj}).movements.(movements{move}).results);
        for d = 1:length(dataTypes)
            k0 = NMF.(subjects{subj}).k0{move,d+1};

            % Muscle Names
            if contains(movements{move},"_R")
                musclesEMG = NMF.(subjects{subj}).meta.muscleNames.EMG_R;
                musclesCalc = NMF.(subjects{subj}).meta.muscleNames.calc_r;
                musclesCalcReduced = NMF.(subjects{subj}).meta.muscleNames.calcReduced_r;
            elseif contains(movements{move},"_L")
                musclesEMG = NMF.(subjects{subj}).meta.muscleNames.EMG_L;
                musclesCalc = NMF.(subjects{subj}).meta.muscleNames.calc_l;
                musclesCalcReduced = NMF.(subjects{subj}).meta.muscleNames.calcReduced_l;
            end

            if k0 == 999
                continue
            else
                t = tiledlayout(k0,2);
                for k = 1:k0
                    nexttile
                    yW = NMF.(subjects{subj}).movements.(movements{move}).results.(dataTypes{d}).(strcat("k",num2str(k0))).W(:,k);
                    plot(yW)
                    xlim([0 100])
                    xticks([25 50 75 100])
                    xlabel('Gait Cycle [%]')
                    xline(60,"","Swing")
                    title(strcat('Activation Patterns Synergy:'," ",num2str(k)))
                    grid on
                    hold off

                    nexttile
                    yH = NMF.(subjects{subj}).movements.(movements{move}).results.(dataTypes{d}).(strcat("k",num2str(k0))).H(k,:);

                    if dataTypes{d} == "EMG"
                        muscles = (musclesEMG);
                    elseif dataTypes{d} == "calc"
                        muscles = (musclesCalc);
                    elseif dataTypes{d} == "calcReduced"
                        muscles = (musclesCalcReduced);
                    end

                    muscles = strrep(muscles,"_"," ");
                    X = categorical(muscles);
                    X = reordercats(X,muscles);

                    bar(X,yH)
                    title(strcat('Muscle Weighting Synergy:'," ",num2str(k)))
                    grid on
                    hold off
                end % k loop (Tile)
            end % k0 999 check loop

            movementTitle = strrep(movements{move},"_"," ");
            t.Title.String = strcat('Minimum Synergies for'," ",movementTitle," ","(",subjects{subj},","," ",dataTypes{d},")");
            t.Title.FontWeight = 'bold';
            t.Title.FontSize = 16;

            filename = strcat("nmf_",subjects{subj},"_",movements{move},"_",dataTypes{d});
            filename = fullfile(nmf_vis_dir,filename);
            savefig(filename)

        end % Datatypes Loop (Tiled Layout)
    end % NMF Loop
end % Subject Loop



