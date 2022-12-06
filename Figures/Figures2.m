
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

%% Visualize the Number of Trials
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
writetable(k0_tbl,fullfile(figDir,"k0.xlsx"))

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

%% Visualize Symmetry

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
        y = [yL;yR]+.1*pat;
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

%% Visualize k0 (Affected vs Unaffected)
movements = k0_tbl.Properties.VariableNames(2:end);
movements_aff = strrep(movements,"_L","_Affected");
movements_aff = strrep(movements_aff,"_R","_Unaffected");


% Hard code a list of movement patterns
patterns = ["gait_normal","gait_fast","s2s","squat","cmj"];

% Intialize output arrays
k0_aff_arr = NaN(size(k0_arr));
k0_sym_arr = NaN(size(k0_arr));

% Create the k0_aff and k0_sym arrays (convert _L, _R -> _Affected,_Unaffected)
for subj = 1:length(subjects)
    for pat = 1:length(patterns)
        % Get L& R Movement Patterns
        mvmt_LR = movements(contains(movements,patterns(pat)));
        % Check that there is a k0 fir both sides
        if length(mvmt_LR) == 1
            continue
        end % Both side check
        
        % Get the k0 values for each movement
        k0_mvmt_L = k0_tbl.(mvmt_LR{1})(subj);
        k0_mvmt_R = k0_tbl.(mvmt_LR{2})(subj);
        
        % Determine the affected side
        if k0_mvmt_L == k0_mvmt_R
           k0_sym_arr(subj,pat*2-1) = k0_mvmt_L;
           k0_sym_arr(subj,pat*2) = k0_mvmt_R;
        elseif k0_mvmt_L < k0_mvmt_R
            k0_aff_arr(subj,pat*2-1) = k0_mvmt_L;
            k0_aff_arr(subj,pat*2) = k0_mvmt_R;
        elseif k0_mvmt_L > k0_mvmt_R
            k0_aff_arr(subj,pat*2-1) = k0_mvmt_R;
            k0_aff_arr(subj,pat*2) = k0_mvmt_L;
        end % Affected Side Loop 
    end % Patterns Loop
end % Subj Loop

% Convert the arrays to tables
subj_tbl = array2table(subjects,"VariableNames","Subject");

k0_aff_tbl = array2table(k0_aff_arr,"VariableNames",movements_aff);
k0_aff_tbl = [subj_tbl k0_aff_tbl];
writetable(k0_aff_tbl,fullfile(figDir,"k0_aff.xlsx"))

k0_sym_tbl = array2table(k0_sym_arr,"VariableNames",movements);
k0_sym_tbl = [subj_tbl k0_sym_tbl];
writetable(k0_sym_tbl,fullfile(figDir,"k0_sym.xlsx"))

% Visualize K0 (affected/unaffected)
movements_disp = strrep(movements_aff,"_"," ");
movements_disp = upper(movements_disp);

t = tiledlayout('flow');
for move = 1:width(k0_aff_arr)
    nexttile
    hold on
    histogram(k0_aff_arr(:,move))
    xlim([0 6])
    xticks(0:6)
    ylim([0 4])
    yticks(0:4)
    title(movements_disp{move})
    subtitle(strcat("n=",num2str(nnz(~isnan(k0_aff_arr(:,move))))))
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

filename = fullfile(figDir,"k0_aff_Histograms");
savefig(filename)

%% Normal Gait, Affected Side, 2 Synergies

W_data = [{NMF.CP15.movements.gait_normal_L.results.EMG.k2.W(:,1)},{NMF.CP15.movements.gait_normal_L.results.EMG.k2.W(:,2)};
          {NMF.CP16.movements.gait_normal_L.results.EMG.k2.W(:,1)},{NMF.CP16.movements.gait_normal_L.results.EMG.k2.W(:,2)};
           {NMF.CP3.movements.gait_normal_L.results.EMG.k2.W(:,1)}, {NMF.CP3.movements.gait_normal_L.results.EMG.k2.W(:,2)}];

H_data = [{NMF.CP15.movements.gait_normal_L.results.EMG.k2.H(1,:)},{NMF.CP15.movements.gait_normal_L.results.EMG.k2.H(2,:)};
          {NMF.CP16.movements.gait_normal_L.results.EMG.k2.H(1,:)},{NMF.CP16.movements.gait_normal_L.results.EMG.k2.H(2,:)};
           {NMF.CP3.movements.gait_normal_L.results.EMG.k2.H(1,:)}, {NMF.CP3.movements.gait_normal_L.results.EMG.k2.H(2,:)}];
        
% Muscle Names
muscles = ["REF","VAL","BIF","MEH","TIA","GAS","SOL","GLU"];
X = categorical(muscles);
X = reordercats(X,muscles);

t = tiledlayout(2,2); %Adjust for n Synergies

for k = 1:2 % n Synergies
    nexttile
    for subj = 1:height(W_data)
        yW = W_data{subj,k};
        hold on
        plot(yW)
        xlim([0 100])
        xticks([25 50 75 100])
        xlabel('Gait Cycle [%]')
        xline(60,"","Swing")
        ylim([0 1.5])
        title(strcat('Activation Patterns Synergy:'," ",num2str(k)))
        grid on
        hold off
    end

    nexttile
        hold on
        yH = vertcat(H_data{:,k})';
        bar(X,yH)
        ylim([0 1])
        title(strcat('Muscle Weighting Synergy:'," ",num2str(k)))
        grid on
        hold off
end

t.Title.String = "Minimum Synergies for Normal Gait, Afffected Side (EMG)";
t.Title.FontWeight = 'bold';
t.Title.FontSize = 16;

leg = legend("CP15","CP16","CP3");
leg.Layout.Tile = 'north';
leg.Orientation = 'horizontal';


filename = "min_synergies_normal_gait_affected_side_EMG";
filename = fullfile(figDir,filename);
savefig(filename)

%% Normal Gait, Unaffected Side, 4 Synergies

W_data = [{NMF.CP16.movements.gait_normal_L.results.EMG.k4.W(:,1)},{NMF.CP16.movements.gait_normal_L.results.EMG.k4.W(:,2)},{NMF.CP16.movements.gait_normal_L.results.EMG.k4.W(:,3)},{NMF.CP16.movements.gait_normal_L.results.EMG.k4.W(:,4)};
           {NMF.CP4.movements.gait_normal_L.results.EMG.k4.W(:,1)}, {NMF.CP4.movements.gait_normal_L.results.EMG.k4.W(:,2)}, {NMF.CP4.movements.gait_normal_L.results.EMG.k4.W(:,3)},{NMF.CP4.movements.gait_normal_L.results.EMG.k4.W(:,4)};
           {NMF.CP8.movements.gait_normal_L.results.EMG.k4.W(:,3)}, {NMF.CP8.movements.gait_normal_L.results.EMG.k4.W(:,2)}, {NMF.CP8.movements.gait_normal_L.results.EMG.k4.W(:,1)},{NMF.CP8.movements.gait_normal_L.results.EMG.k4.W(:,4)}];

H_data = [{NMF.CP16.movements.gait_normal_L.results.EMG.k4.H(1,:)},{NMF.CP16.movements.gait_normal_L.results.EMG.k4.H(2,:)},{NMF.CP16.movements.gait_normal_L.results.EMG.k4.H(3,:)},{NMF.CP16.movements.gait_normal_L.results.EMG.k4.H(4,:)};
           {NMF.CP4.movements.gait_normal_L.results.EMG.k4.H(1,:)}, {NMF.CP4.movements.gait_normal_L.results.EMG.k4.H(2,:)}, {NMF.CP4.movements.gait_normal_L.results.EMG.k4.H(3,:)},{NMF.CP4.movements.gait_normal_L.results.EMG.k4.H(4,:)};
           {NMF.CP8.movements.gait_normal_L.results.EMG.k4.H(3,:)}, {NMF.CP8.movements.gait_normal_L.results.EMG.k4.H(2,:)}, {NMF.CP8.movements.gait_normal_L.results.EMG.k4.H(1,:)},{NMF.CP8.movements.gait_normal_L.results.EMG.k4.H(4,:)}];

t = tiledlayout(4,2); %Adjust for n Synergies
for k = 1:4 % n Synergies
    nexttile
    for subj = 1:height(W_data)
        yW = W_data{subj,k};
        hold on
        plot(yW)
        xlim([0 100])
        xticks([25 50 75 100])
        xlabel('Gait Cycle [%]')
        xline(60,"","Swing")
        ylim([0 1])
        title(strcat('Activation Patterns Synergy:'," ",num2str(k)))
        grid on
        hold off
    end

    nexttile
        hold on
        yH = vertcat(H_data{:,k})';
        bar(X,yH)
        ylim([0 1])
        title(strcat('Muscle Weighting Synergy:'," ",num2str(k)))
        grid on
        hold off
end

movementTitle = strrep(movements{move},"_"," ");
t.Title.String = "Minimum Synergies for Normal Gait, Unafffected Side (EMG)";
t.Title.FontWeight = 'bold';
t.Title.FontSize = 16;

leg = legend("CP16","CP4","CP8");
leg.Layout.Tile = 'north';
leg.Orientation = 'horizontal';

filename = "min_synergies_normal_gait_unaffected_side_EMG";
filename = fullfile(figDir,filename);
savefig(filename)

%% Fast Gait, Affected Side, 2 Synergies
W_data = [{NMF.CP15.movements.gait_fast_L.results.EMG.k2.W(:,1)},{NMF.CP15.movements.gait_fast_L.results.EMG.k2.W(:,2)};
           {NMF.CP16.movements.gait_fast_L.results.EMG.k2.W(:,1)}, {NMF.CP16.movements.gait_fast_L.results.EMG.k2.W(:,2)};
           {NMF.CP3.movements.gait_fast_L.results.EMG.k2.W(:,1)}, {NMF.CP3.movements.gait_fast_L.results.EMG.k2.W(:,2)}];

H_data = [{NMF.CP15.movements.gait_fast_L.results.EMG.k2.H(1,:)},{NMF.CP15.movements.gait_fast_L.results.EMG.k2.H(2,:)};
           {NMF.CP16.movements.gait_fast_L.results.EMG.k2.H(1,:)}, {NMF.CP16.movements.gait_fast_L.results.EMG.k2.H(2,:)};
           {NMF.CP3.movements.gait_fast_L.results.EMG.k2.H(1,:)}, {NMF.CP3.movements.gait_fast_L.results.EMG.k2.H(2,:)}];

t = tiledlayout(2,2); %Adjust for n Synergies
for k = 1:2 % n Synergies
    nexttile
    for subj = 1:height(W_data)
        yW = W_data{subj,k};
        hold on
        plot(yW)
        xlim([0 100])
        xticks([25 50 75 100])
        xlabel('Gait Cycle [%]')
        xline(60,"","Swing")
        ylim([0 2])
        title(strcat('Activation Patterns Synergy:'," ",num2str(k)))
        grid on
        hold off
    end

    nexttile
        hold on
        yH = vertcat(H_data{:,k})';
        ylim([0 1])
        bar(X,yH)
        title(strcat('Muscle Weighting Synergy:'," ",num2str(k)))
        grid on
        hold off
end

movementTitle = strrep(movements{move},"_"," ");
t.Title.String = "Minimum Synergies for Fast Gait, Afffected Side (EMG)";
t.Title.FontWeight = 'bold';
t.Title.FontSize = 16;

leg = legend("CP15","CP16","CP3");
leg.Layout.Tile = 'north';
leg.Orientation = 'horizontal';

filename = "min_synergies_fast_gait_affected_side_EMG";
filename = fullfile(figDir,filename);
savefig(filename)

%% Normal Gait, Unaffected Side, 4 Synergies
W_data = [{NMF.CP16.movements.s2s_L.results.EMG.k3.W(:,2)},{NMF.CP16.movements.s2s_L.results.EMG.k3.W(:,3)},{NMF.CP16.movements.s2s_L.results.EMG.k3.W(:,1)},;
           {NMF.CP16.movements.s2s_R.results.EMG.k3.W(:,1)}, {NMF.CP16.movements.s2s_R.results.EMG.k3.W(:,2)}, {NMF.CP16.movements.s2s_R.results.EMG.k3.W(:,3)},;
           {NMF.CP8.movements.s2s_L.results.EMG.k3.W(:,1)}, {NMF.CP8.movements.s2s_L.results.EMG.k3.W(:,3)}, {NMF.CP8.movements.s2s_L.results.EMG.k3.W(:,2)};
           {NMF.CP8.movements.s2s_R.results.EMG.k3.W(:,2)}, {NMF.CP8.movements.s2s_R.results.EMG.k3.W(:,1)}, {NMF.CP8.movements.s2s_R.results.EMG.k3.W(:,3)},];

H_data = [{NMF.CP16.movements.s2s_L.results.EMG.k3.H(2,:)},{NMF.CP16.movements.s2s_L.results.EMG.k3.H(3,:)},{NMF.CP16.movements.s2s_L.results.EMG.k3.H(1,:)};
           {NMF.CP16.movements.s2s_R.results.EMG.k3.H(1,:)}, {NMF.CP16.movements.s2s_R.results.EMG.k3.H(2,:)}, {NMF.CP16.movements.s2s_R.results.EMG.k3.H(3,:)};
           {NMF.CP8.movements.s2s_L.results.EMG.k3.H(1,:)}, {NMF.CP8.movements.s2s_L.results.EMG.k3.H(3,:)}, {NMF.CP8.movements.s2s_L.results.EMG.k3.H(2,:)};
           {NMF.CP8.movements.s2s_R.results.EMG.k3.H(2,:)}, {NMF.CP8.movements.s2s_R.results.EMG.k3.H(1,:)}, {NMF.CP8.movements.s2s_R.results.EMG.k3.H(3,:)}];

t = tiledlayout(3,2); %Adjust for n Synergies
for k = 1:3 % n Synergies
    nexttile
    for subj = 1:height(W_data)
        yW = W_data{subj,k};
        hold on
        plot(yW)
        xlim([0 100])
        xticks([25 50 75 100])
        xlabel('Movement Cycle [%]')
        ylim([0 1])
        title(strcat('Activation Patterns Synergy:'," ",num2str(k)))
        grid on
        hold off
    end

    nexttile
        hold on
        yH = vertcat(H_data{:,k})';
        bar(X,yH)
        ylim([0 1])
        title(strcat('Muscle Weighting Synergy:'," ",num2str(k)))
        grid on
        hold off
end

movementTitle = strrep(movements{move},"_"," ");
t.Title.String = "Minimum Synergies for Sit-to-Stand, No Affected Side (EMG)";
t.Title.FontWeight = 'bold';
t.Title.FontSize = 16;

leg = legend("CP16 (L)","CP16 (R)","CP8 (L)","CP8 (R)");
leg.Layout.Tile = 'north';
leg.Orientation = 'horizontal';

filename = "min_synergies_s2s_no_affected_side_EMG";
filename = fullfile(figDir,filename);
savefig(filename)
