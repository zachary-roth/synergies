
close all
clear
clc

disp('Select the NNMF.mat structure')
[file, path] = uigetfile('*.mat','Select the NNMF.mat structure');
load(fullfile(path,file))
cd(path)

load(fullfile(path,'GaitData.mat'))

muscles = GaitData.calc.L.muscles;
muscles = strrep(muscles,"_l","");
muscles = strrep(muscles,"_"," ");
X = categorical(muscles);
X = reordercats(X,muscles);

%%
t = tiledlayout(7,2);

%% 1,1
nexttile
Y1 = NNMF.calc.L.k7.W(:,3);
Y2 = NNMF.calc.R.k5.W(:,1);

hold on
plot(Y1)
plot(Y2)
xlim([0 100])
title('Activation Patterns')
ylabel('Synergy 1')
grid on
hold off

% 1,2
nexttile
for h = 1:43
    Y(h,1) = NNMF.calc.L.k7.H(3,h);
    Y(h,2) = NNMF.calc.R.k5.H(1,h);
end

hold on
bar(X,Y)
ylim([0 1])
title('Muscle Weightings')
grid on
hold off

%% 2,1
nexttile
Y1 = NNMF.calc.L.k7.W(:,1);
Y2 = NNMF.calc.R.k5.W(:,2);

hold on
plot(Y1)
plot(Y2)
xlim([0 100])
ylabel('Synergy 2')
grid on
hold off

% 2,2
nexttile
for h = 1:43
    Y(h,1) = NNMF.calc.L.k7.H(2,h);
    Y(h,2) = NNMF.calc.R.k5.H(2,h);
end

hold on
bar(X,Y)
ylim([0 1])
grid on
hold off

%% 3,1
nexttile
Y1 = NNMF.calc.L.k7.W(:,2);
Y2 = NNMF.calc.R.k5.W(:,3);

hold on
plot(Y1)
plot(Y2)
xlim([0 100])
ylabel('Synergy 3')
grid on
hold off

% 3,2
nexttile
for h = 1:43
    Y(h,1) = NNMF.calc.L.k7.H(2,h);
    Y(h,2) = NNMF.calc.R.k5.H(3,h);
end

hold on
bar(X,Y)
ylim([0 1])
grid on
hold off

%% 4,1
nexttile
Y1 = NNMF.calc.L.k7.W(:,4);
Y2 = NNMF.calc.R.k5.W(:,4);

hold on
plot(Y1)
plot(Y2)
xlim([0 100])
ylabel('Synergy 4')
grid on
hold off

% 4,2
nexttile
for h = 1:43
    Y(h,1) = NNMF.calc.L.k7.H(4,h);
    Y(h,2) = NNMF.calc.R.k5.H(4,h);
end

hold on
bar(X,Y)
ylim([0 1])
grid on
hold off

%% 5,1
nexttile
Y1 = NNMF.calc.L.k7.W(:,6);
Y2 = NNMF.calc.R.k5.W(:,5);

hold on
plot(Y1)
plot(Y2)
xlim([0 100])
ylabel('Synergy 5')
grid on
hold off

% 5,2
nexttile
for h = 1:43
    Y(h,1) = NNMF.calc.L.k7.H(6,h);
    Y(h,2) = NNMF.calc.R.k5.H(5,h);
end

hold on
bar(X,Y)
ylim([0 1])
grid on
hold off

%% 6,1
nexttile
Y1 = NNMF.calc.L.k7.W(:,5);


hold on
plot(Y1)
xlim([0 100])
ylabel('Synergy 6')
grid on
hold off

% 6,2
nexttile
for h = 1:43
    Y(h,1) = NNMF.calc.L.k7.H(5,h);
    Y(h,2) = 0;
end

hold on
bar(X,Y)
ylim([0 1])
grid on
hold off

%% 7,1
nexttile
Y1 = NNMF.calc.L.k7.W(:,7);


hold on
plot(Y1)
xlim([0 100])
xlabel('Gait Cycle [%]')
ylabel('Synergy 7')
grid on
hold off

% 6,2
nexttile
for h = 1:43
    Y(h,1) = NNMF.calc.L.k7.H(7,h);
    Y(h,2) = 0;
end

hold on
bar(X,Y)
xlabel('Muscles')
ylim([0 1])
grid on
hold off

%%
leg = legend('Left Gaits/Muscles','Right Gaits/Muscles','Orientation','Horizontal');
leg.Layout.Tile = 'north';

t.Title.String = "Muscle Synergies: Calculated Activations";
t.Title.FontWeight = 'bold';
t.Title.FontSize = 16;

filename = fullfile(path,"Figures","calc_NNMF_Presentation");
savefig(filename)
