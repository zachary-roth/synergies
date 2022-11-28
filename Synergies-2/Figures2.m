
close all; clear; clc

%% Import Data

% Select the synergies-data folder
disp('Select the synergies2 data folder')
synergiesPath = uigetdir('','Select the synergies-data folder');
cd(synergiesPath)

load(fullfile(synergiesPath,"Results","MoveData.mat"))
load(fullfile(synergiesPath,"Results","NMF.mat"))


% Delete the old Figures Folder
figDir = fullfile(synergiesPath,"Results","Figures");
if exist(figDir,"dir") == 7
    rmdir(figDir,'s')
end
mkdir(figDir)

%% Visualize the Number 
% Trials per Subject
movements = ["gait_normal_L","gait_normal_R","gait_fast_L","gait_fast_R","s2s","squat","cmj"];
subjects = fieldnames(MoveData);

trials_arr = NaN(length(subjects),length(movements));

% Build the trials array
for subj = 1:length(subjects)
    for move = 1:length(movements)
        if ismember(movements{move},fieldnames(MoveData.(subjects{subj}).movements))
            n_trials = nnz(count(MoveData.(subjects{subj}).meta.trials.movement,movements{move}));
            trials_arr(subj,move) = n_trials;
        else
            continue
        end % Check Loop
    end % Movement Loop
end % Subject Loop

trial_tbl = array2table(trials_arr,'VariableNames',movements);
subj_tbl = array2table(subjects,"VariableNames","Subject");
trial_tbl = [subj_tbl trial_tbl];
writetable(trial_tbl,fullfile(figDir,"nTrials.xlsx"))

% Visualize Trials per Subject

movements_disp = strrep(movements,"_"," ");
movements_disp = upper(movements_disp);

x = categorical(subjects);
x = reordercats(x,{'CP3','CP4','CP5','CP8','CP15','CP16','CP18'});

t = tiledlayout('flow');
for move = 1:width(trials_arr)
    nexttile
    hold on
    bar(x,trials_arr(:,move))
    ylim([0 7])
    title(movements_disp{move})
    grid on
    hold off
end % Subject Loop

t.Title.String = "Trials Per Subject";
t.Title.FontSize = 16;
t.Title.FontWeight = 'bold';

t.XLabel.String = "Subjects";
t.XLabel.FontSize = 16;

t.YLabel.String = "Number of Trials";
t.YLabel.FontSize = 16;

filename = fullfile(figDir,"trials_per_subject");
savefig(filename)

% Visualize Trials overall

x = categorical(movements_disp);
x = reordercats(x,{'GAIT NORMAL L','GAIT NORMAL R','GAIT FAST L','GAIT FAST R','S2S','SQUAT','CMJ'});

figure
hold on
b = bar(x,trials_arr','stacked');
xlabel("Movments")
ylabel("Number of Trials")
title("Number of Trials Overall")
legend(subjects)
grid on
grid minor
hold off

filename = fullfile(figDir,"trials_overall");
savefig(filename)

%% Visualize k0
movements = ["gait_normal_L","gait_normal_R","gait_fast_L","gait_fast_R","s2s_L","s2s_R","squat_L","squat_R","cmj_L","cmj_R"];
subjects = fieldnames(NMF);

k0_arr = NaN(length(subjects),length(movements));

% Build the k0_array
for subj = 1:length(subjects)
    k0 = NMF.(subjects{subj}).k0;
    for move = 1:length(movements)
        if ismember(movements{move},k0(:,1))
        move_row = strcmp(k0(:,1),movements{move});
        k0_arr(subj,move) = k0{move_row,2};
        else
            continue
        end % Check Loop
    end % Movement Loop
end % Subjects Loop

k0_tbl = array2table(k0_arr,'VariableNames',movements);
subj_tbl = array2table(subjects,"VariableNames","Subject");
k0_tbl = [subj_tbl k0_tbl];
writetable(trial_tbl,fullfile(figDir,"k0.xlsx"))

movements_disp = strrep(movements,"_"," ");
movements_disp = upper(movements_disp);

t = tiledlayout('flow');
for move = 1:width(k0_arr)
    nexttile
    hold on
    histogram(k0_arr(:,move))
    xlim([0 6])
    xticks(0:6)
    ylim([0 4])
    yticks(0:4)
    title(movements_disp{move})
    subtitle(strcat("n=",num2str(nnz(~isnan(k0_arr(:,move))))))
    grid on
    hold off
end % Tile Layout 

t.Title.String = "Minimum Synergies per Movement";
t.Title.FontSize = 20;
t.Title.FontWeight = 'bold';

t.XLabel.String = "Minimum Synergies Required to Explain 90% VAF";
t.XLabel.FontSize = 16;
t.XLabel.FontWeight = 'bold';

t.YLabel.String = "Number of Subjects";
t.YLabel.FontSize = 16;
t.YLabel.FontWeight = 'bold';

filename = fullfile(figDir,"k0_Histograms");
savefig(filename)

%% Symmetry Plot

% k0 symetries per subject
patterns = ["gait_normal","gait_fast","s2s","squat","cmj"];

patterns_disp = strrep(patterns,"_"," ");
patterns_disp = upper(patterns_disp);

t = tiledlayout('flow');

for subj = 1:length(subjects)
    nexttile
    for pat = 1:length(patterns)
        Li = find(contains(movements,patterns{pat}),1,"first");
        Ri = find(contains(movements,patterns{pat}),1,"last");
        yL = k0_arr(subj,Li);
        yR = k0_arr(subj,Ri);
        x = [1; 2];
        y = [yL;yR]+.05*pat;
        hold on
        plot(x,y)
    end % Pattern Loop
    xlim([.75 2.25])
    xticks(1:2)
    xticklabels(["k0 Left", "k0 Right"])
    ylim([0 6])
    yticks([])
    title(subjects{subj})
    hold off
end % Subject Loop

leg = legend(patterns_disp);
leg.Layout.Tile = 'north';
leg.Orientation = 'horizontal';

t.Title.String = "Left/Right Synergy Symmetries per Subject";
t.Title.FontSize = 20;
t.Title.FontWeight = 'bold';

t.XLabel.String = "Movements (Left and Right Legs)";
t.XLabel.FontSize = 16;
t.XLabel.FontWeight = 'bold';

t.YLabel.String = "Min Synergies";
t.YLabel.FontSize = 16;
t.YLabel.FontWeight = 'bold';

filename = fullfile(figDir,"symmetry_subj");
savefig(filename)

% k0 symetries per movement
t = tiledlayout('flow');

for pat = 1:length(patterns)
    nexttile
    for subj = 1:length(subjects)
        Li = find(contains(movements,patterns{pat}),1,"first");
        Ri = find(contains(movements,patterns{pat}),1,"last");
        yL = k0_arr(subj,Li);
        yR = k0_arr(subj,Ri);
        x = [1; 2];
        y = [yL;yR];
        hold on
        plot(x,y)
    end % Pattern Loop
    xlim([.75 2.25])
    xticks(1:2)
    xticklabels(["k0 Left", "k0 Right"])
    ylim([0 6])
    yticks([])
    title(patterns_disp{pat})
    hold off
end % Subject Loop

leg = legend(subjects);
leg.Layout.Tile = 'north';
leg.Orientation = 'horizontal';

t.Title.String = "Left/Right Synergy Symmetries per Movement";
t.Title.FontSize = 20;
t.Title.FontWeight = 'bold';

t.XLabel.String = "Movements (Left and Right Legs)";
t.XLabel.FontSize = 16;
t.XLabel.FontWeight = 'bold';

t.YLabel.String = "Min Synergies";
t.YLabel.FontSize = 16;
t.YLabel.FontWeight = 'bold';

filename = fullfile(figDir,"symmetry_move");
savefig(filename)