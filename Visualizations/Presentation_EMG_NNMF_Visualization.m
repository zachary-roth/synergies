
close all
clear
clc

disp('Select the NNMF.mat structure')
[file, path] = uigetfile('*.mat','Select the NNMF.mat structure');
load(fullfile(path,file))
cd(path)

muscles = [{'REF'}, {'VAL'}, {'BIF'}, {'MEH'}, {'TIA'}, {'GAS'}, {'SOL'}, {'GLU'}];
X = categorical(muscles);
X = reordercats(X,muscles);

t = tiledlayout(3,2);

nexttile
Y1 = NNMF.EMG.L.k3.W(:,1);
Y2 = NNMF.EMG.R.k3.W(:,3);

hold on
plot(Y1)
plot(Y2)
xlim([0 100])
title('Activation Patterns')
ylabel('Synergy 1')
grid on
hold off

nexttile
for h = 1:8
    Y(h,1) = NNMF.EMG.L.k3.H(1,h);
    Y(h,2) = NNMF.EMG.R.k3.H(3,h);
end

hold on
bar(X,Y)
ylim([0 1])
title('Muscle Weightings')
grid on
hold off

nexttile
Y1 = NNMF.EMG.L.k3.W(:,2);
Y2 = NNMF.EMG.R.k3.W(:,2);

hold on
plot(Y1)
plot(Y2)
xlim([0 100])
ylabel('Synergy 2')
grid on
hold off

nexttile
for h = 1:8
    Y(h,1) = NNMF.EMG.L.k3.H(2,h);
    Y(h,2) = NNMF.EMG.R.k3.H(2,h);
end

hold on
bar(X,Y)
ylim([0 1])
grid on
hold off

nexttile
Y1 = NNMF.EMG.L.k3.W(:,3);
Y2 = NNMF.EMG.R.k3.W(:,1);

hold on
plot(Y1)
plot(Y2)
xlim([0 100])
xlabel('Gait Cycle [%]')
ylabel('Synergy 3')
grid on
hold off

nexttile
for h = 1:8
    Y(h,1) = NNMF.EMG.L.k3.H(3,h);
    Y(h,2) = NNMF.EMG.R.k3.H(1,h);
end

hold on
bar(X,Y)
xlabel('Muscles')
ylim([0 1])
grid on
hold off

leg = legend('Left Gaits/Muscles','Right Gaits/Muscles','Orientation','Horizontal');
leg.Layout.Tile = 'north';

t.Title.String = "Muscle Synergies: EMG";
t.Title.FontWeight = 'bold';
t.Title.FontSize = 16;

filename = fullfile(path,"Figures","EMG_NNMF_Presentation");
savefig(filename)
