
close all
clear
clc

disp('Select the SYNERGIES 3 GaitData.mat structure')
[file, path] = uigetfile('*.mat','Select the GaitData.mat structure');
load(fullfile(path,file));
cd(path)

MRI.calcReduced = GaitData.calcReduced;

disp('Select the SYNERGIES 1 GaitData.mat structure')
[synergies1file, synergies1path] = uigetfile('*.mat','Select the GaitData.mat structure');
load(fullfile(synergies1path,synergies1file));
cd(path)

GaitData.MRI.calcReduced = MRI.calcReduced;

left = tiledlayout(4,2);

for m = 1:8
    nexttile
    EMG = GaitData.EMG.L.concat(:,m);
    generic_model = GaitData.calcReduced.L.concat(:,m);
    MRI_model = GaitData.MRI.calcReduced.L.concat(:,m);
    hold on
    plot(EMG)
    plot(generic_model)
    plot(MRI_model)
    xlim([0 101])
    xline(60,'--','Swing')
    xticks([25 50 75 100])
    xlabel('Gait Cycle [%]')
    ylabel('Relative Activation [%]')
    title(GaitData.EMG.L.muscles(m))
    grid on
    hold off
end

left.Title.String = 'Observed vs Calculated Muscle Activations';
left.Title.FontWeight = 'bold';
left.Title.FontSize = 16;

left.Subtitle.String = strcat('Left Gait, Left Leg Muscles');
left.Subtitle.FontSize = 14;

leg = legend('EMG','Generic Model','MRI-Based Model','Orientation', 'Horizontal');
leg.Layout.Tile = 'north';

filename = fullfile(path,"Figures","Validation_EMG_Calc_Left_Side");
savefig(filename)

right = tiledlayout(4,2);

for m = 1:8
    nexttile
    EMG = GaitData.EMG.R.concat(:,m);
    generic_model = GaitData.calcReduced.R.concat(:,m);
    MRI_model = GaitData.MRI.calcReduced.R.concat(:,m);
    hold on
    plot(EMG)
    plot(generic_model)
    plot(MRI_model)
    xlim([0 101])
    xline(60,'--','Swing')
    xticks([25 50 75 100])
    xlabel('Gait Cycle [%]')
    ylabel('Relative Activation [%]')
    title(GaitData.EMG.L.muscles(m))
    grid on
    hold off
end

right.Title.String = 'Observed vs Calculated Muscle Activations';
right.Title.FontWeight = 'bold';
right.Title.FontSize = 16;

right.Subtitle.String = strcat('Right Gait, Right Leg Muscles');
right.Subtitle.FontSize = 14;

leg = legend('EMG','Generic Model','MRI-Based Model','Orientation', 'Horizontal');
leg.Layout.Tile = 'north';

filename = fullfile(path,"Figures","Validation_EMG_Calc_Right_Side");
savefig(filename)