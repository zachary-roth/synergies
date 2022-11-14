
close all
clear
clc

disp('Select the GaitData.mat structure')
[file, path] = uigetfile('*.mat','Select the GaitData.mat structure');
load(fullfile(path,file))
cd(path)

raw = GaitData.EMG.L.raw.raw1(:,2);
trim = GaitData.EMG.L.trim.trim1(:,1);
notch = GaitData.EMG.L.notch.notch1(:,1);
bandbass = GaitData.EMG.L.bandPass.bandPass1(:,1);
rect = GaitData.EMG.L.rect.rect1(:,1);
lowPass = GaitData.EMG.L.lowPass.lowPass1(:,1);
norm = GaitData.EMG.L.norm.norm1(:,1);
resample = GaitData.EMG.L.resample.resample1(:,1);
concat = GaitData.EMG.L.concat(:,1);

t = tiledlayout(2,8);

nexttile
hold on
plot(raw)
xlim([0 length(raw)])
ylim([-0.5 0.5])
xlabel('Time [ms]')
ylabel('Muscle Activation [mV]')
title('Raw EMG (Entire Gait Trial)')
grid on
hold off

nexttile
hold on
plot(trim)
xlim([0 length(trim)])
ylim([-0.3 0.3])
xlabel('Time [ms]')
ylabel('Muscle Activation [mV]')
title('Trimmed EMG (One Gait Cycle)')
grid on
hold off

nexttile
hold on
plot(notch)
xlim([0 length(notch)])
ylim([-0.3 0.3])
xlabel('Time [ms]')
ylabel('Muscle Activation [mV]')
title('Notch Filter (50 Hz)')
grid on
hold off

nexttile
hold on
plot(bandbass)
xlim([0 length(bandbass)])
ylim([-0.3 0.3])
xlabel('Time [ms]')
ylabel('Muscle Activation [mV]')
title('Bandpass Filter (20-400 Hz)')
grid on
hold off

nexttile
hold on
plot(rect)
xlim([0 length(rect)])
ylim([0 0.3])
xlabel('Time [ms]')
ylabel('Muscle Activation [mV]')
title('Rectification')
grid on
hold off

nexttile
hold on
plot(lowPass)
xlim([0 length(lowPass)])
ylim([0 0.12])
xlabel('Time [ms]')
ylabel('Muscle Activation [mV]')
title('Lowpass Filter (15 Hz)')
grid on
hold off

nexttile
hold on
plot(norm)
xlim([0 length(norm)])
ylim([0 0.12])
xlabel('Time [ms]')
ylabel('Relative Activation [%]')
title('Normalize')
grid on
hold off

nexttile
hold on
plot(resample)
xlim([0 length(resample)])
ylim([0 0.12])
xlabel('Gait Cycle [%]')
ylabel('Relative Activation [%]')
title('Resample')
grid on
hold off

nexttile([1 8])
hold on
plot(concat)
xlim([0 length(concat)])
xticks([25 50 75 101 126 151 176 202 226 252 277])
ylim([0 0.15])
xlabel('Gait Cycle [%]')
ylabel('Relative Activation [%]')
title('Concatenation')
xline(101,':','Gait 2')
xline(202,':','Gait 3')
grid on
hold off

t.Title.String = 'EMG Processing';
t.Title.FontWeight = 'bold';
t.Title.FontSize = 16;

t.Subtitle.String = 'Left Gait, Left Rectus Femoris';
t.Subtitle.FontSize = 14;

filename = fullfile(path,"Figures","EMG_Processing");
savefig(filename)