close all; clear; clc

%%
disp('Select the synergies-data folder')
dataDir = uigetdir("",'Select the synergies-data folder');

plot_titles = [
    "m. Rectus Femoris";
    "m. Vastus Lateralis"; 
    "m. Biceps Femoris";
    "m. Semitendinosus";
    "m. Tibialis Anterior";
    "m. Gastrocnemius";
    "m. Soleus";
    "m. Gluteus Medius"
    ];

%%
synergies1_NMF_path = fullfile(dataDir,"Synergies1","Results","NMF.mat");
load(synergies1_NMF_path)

gen_EMG_A = NMF.CP3.movements.gait_normal_L.data.EMG;
gen_EMG_WH = NMF.CP3.movements.gait_normal_L.results.EMG.k2.WH;

gen_calcReduced_A = NMF.CP3.movements.gait_normal_L.data.calcReduced;
gen_calcReduced_WH = NMF.CP3.movements.gait_normal_L.results.calcReduced.k5.WH;

t = tiledlayout(2,4);

for m = 1:8
    nexttile
    hold on

    % EMG
    plot(gen_EMG_A(:,m),"Color","#0072BD")
    plot(gen_EMG_WH(:,m),"LineStyle","--","Color","#0072BD")

    % gen
    plot(gen_calcReduced_A(:,m),"Color","#D95319")
    plot(gen_calcReduced_WH(:,m),"LineStyle","--","Color","#D95319")

    xlim([0 100])

    xline(60,"","Swing","LabelOrientation","horizontal")

    ylim([0 1])

    title(plot_titles(m))
    grid on
    hold off

end % Tile Loop (Muscle)

leg = legend("A Observed","WH Observed","A Estimated (Gen. Model)","WH Estimated (Gen. Model)");
leg.Layout.Tile = "north";
leg.Orientation = "horizontal";

t.Title.String = 'Observed vs Estimated Activations, with Reconstructions';
t.Title.FontWeight = 'bold';
t.Title.FontSize = 24;

t.Subtitle.String = strcat('Left Gait, Left Leg Muscles');
t.Subtitle.FontSize = 18;

t.XLabel.String = "Gait Cycle [%]";
t.XLabel.FontSize = 22;

t.YLabel.String = "Relative Muscle Activation [%]";
t.YLabel.FontSize = 22;

t.TileSpacing = "compact";
t.Padding = "tight";

filename = fullfile(dataDir,"Synergies1","Results","Visualizations","Presentation_A_WH_gen");
savefig(filename)


%%

synergies3_NMF_path = fullfile(dataDir,"Synergies3","Results","NMF.mat");
load(synergies3_NMF_path)

MRI_calcReduced_A = NMF.CP3.movements.gait_normal_L.data.calcReduced;
MRI_calcReduced_WH = NMF.CP3.movements.gait_normal_L.results.calcReduced.k5.WH;



t = tiledlayout(2,4);

for m = 1:8
    nexttile
    hold on

    % EMG
    plot(gen_EMG_A(:,m),"Color","#0072BD")
    plot(gen_EMG_WH(:,m),"LineStyle","--","Color","#0072BD")

    % gen
    plot(gen_calcReduced_A(:,m),"Color","#D95319")
    plot(gen_calcReduced_WH(:,m),"LineStyle","--","Color","#D95319")

    % MRI
    plot(MRI_calcReduced_A(:,m),"Color","#77AC30")
    plot(MRI_calcReduced_WH(:,m),"LineStyle","--","Color","#77AC30")

    xlim([0 100])

    xline(60,"","Swing","LabelOrientation","horizontal")

    ylim([0 1])

    title(plot_titles(m))
    grid on
    hold off

end % Tile Loop (Muscle)

leg = legend("A Observed","WH Observed","A Estimated (Gen. Model)","WH Estimated (Gen. Model)","A Estimated (MRI Model)","WH Estimated (MRI Model)");
leg.Layout.Tile = "north";
leg.Orientation = "horizontal";

t.Title.String = 'Observed vs Estimated Activations, with Reconstructions';
t.Title.FontWeight = 'bold';
t.Title.FontSize = 24;

t.Subtitle.String = strcat('Left Gait, Left Leg Muscles');
t.Subtitle.FontSize = 18;

t.XLabel.String = "Gait Cycle [%]";
t.XLabel.FontSize = 22;

t.YLabel.String = "Relative Muscle Activation [%]";
t.YLabel.FontSize = 22;

t.TileSpacing = "compact";
t.Padding = "tight";

filename = fullfile(dataDir,"Synergies3","Results","Visualizations","Presentation_A_WH_MRI");
savefig(filename)

% EMG k0 = 2, gen k0 = 5, MRI k0 = 5   "#77AC30"