close all; clear; clc;

% Select the synergies Data folder
disp('Select the synergies Data folder')
synergiesPath = uigetdir('','Select the synergies Data folder');
cd(synergiesPath)

dataDir = fullfile(synergiesPath,"Data");
synergies2Dir = fullfile(synergiesPath,"Synergies2");

subjectsAll = dir(fullfile(dataDir,"CP*"));
subjectsAll = {subjectsAll.name};

subjects = dir(fullfile(synergies2Dir,"CP*"));
subjects = {subjects.name};

for subj = 1:length(subjects)
    % Get the index of the Synergies2 subject from the list of all subjects
    iAll = find(strcmp(subjects{subj},subjectsAll));

    % Import the meta data
    meta = readtable(fullfile(synergies2Dir,subjects{subj},strcat(subjects{subj},"_trials")));

    for f = 1:height(meta)
        % Import the c3d file
        c3dName = strcat(subjects{subj},"_T0",meta.ID{f},".c3d");
        c3dfile = fullfile(dataDir,subjects{subj},"Vicon",c3dName);
        [Markers,MLabels,VideoFrameRate,AnalogSignals,ALabels,AUnits,AnalogFrameRate,Event,ParameterGroup,CameraInfo]= readC3D(c3dfile, []);

        % Fill the marker gaps
        MarkersNaN = Markers;
        MarkersNaN(Markers==0)= NaN;
        gaplim = 10; % can be adjusted
        MarkersFilled = FillGapsMarkers(MarkersNaN,VideoFrameRate,gaplim);
        Markers = MarkersFilled;

        if meta.side(f) == "both"
            meta.IC1(f) = 0.015;
            meta.IC2(f) = length(AnalogSignals)/AnalogFrameRate;
        else
            if meta.side(f) == "R"
                idx_H = find(ismember(MLabels, 'RHEE'));
                Heel = Markers(:,idx_H*3);
            elseif  meta.side(f) == "L"
                idx_H = find(ismember(MLabels, 'LHEE'));
                Heel = Markers(:,idx_H*3);
            end % Heel Loop

            % Find all the possible ICs based on the heel marker data
            possibleIC_TF = islocalmin(Heel,"MinProminence",2);
            possibleICs = find(possibleIC_TF);
            possibleIC_heelHeights = Heel(possibleIC_TF);

            % Find the possible IC1s based on the force plate data

            % Set a minimum Force Threshold
            threshold = 1;

            fpNames = ["Force.Fy1","Force.Fy2"];
            fp_Check = zeros(2,4); % IC1 | IC1_heelHeight | diff_fp | fp_closestIndex

            fp_data = [];

            for fp = 1:2
                % Get the GRFy (left or right) col
                fp_idx = ismember(ALabels, fpNames(fp));

                % Resample the GRFy data
                x = (1:1:length(AnalogSignals))'; % sample points
                v = AnalogSignals(:,fp_idx);
                xq = linspace(1,length(AnalogSignals),length(Markers)); % query point
                vq = interp1(x,v,xq);

                % Write the resampled data to a temp storage array
                fp_data(:,fp) = vq;

                IC1 = find(fp_data(:,fp) > threshold, 1, 'first');
                IC1_heelHeight = Heel(IC1); %
                [diff_fp,fp_closestIndex] = min(abs(Heel(possibleIC_TF)-IC1_heelHeight));
                fp_Check(fp,1) = IC1;
                fp_Check(fp,2) = IC1_heelHeight;
                fp_Check(fp,3) = diff_fp;
                fp_Check(fp,4) = fp_closestIndex;
            end

            if fp_Check(1,3) < fp_Check(2,3)
                fp_data = fp_data(:,1);
                IC1 = fp_Check(1,1);
                IC1_heelHeight = fp_Check(1,2);
            else
                fp_data = fp_data(:,2);
                IC1 = fp_Check(2,1);
                IC1_heelHeight = fp_Check(2,2);
            end
            
            % Plot the Kinematic and Force Plate Data
            x = (1:length(Heel));
            
            fig = figure;
            subplot(2,1,1)
            hold on
            plot(x, Heel,'LineWidth',1)
            plot(x(possibleIC_TF),Heel(possibleIC_TF),"*r")
            plot(IC1,IC1_heelHeight,"^g",'MarkerFaceColor','g','MarkerSize',10)
            xlim([0 length(Heel)])
            xlabel('Time')
            ylabel('Heel Marker Height')
            titleString = (strcat(subjects{subj}," ",meta.movement{f}," ","Trial_",meta.ID{f}));
            titleString = strrep(titleString,"_"," ");
            titleString = upper(titleString);
            title(titleString)
            subtitle('Heel Marker')
            hold off

            subplot(2,1,2)
            hold on
            plot(fp_data)
            xlim([0 length(Heel)])
            xlabel('Time')
            ylabel('Force [N]')
            subtitle('GRFy')
            hold off
            
            fig.WindowState = 'maximized';

            [IC2,IC2_heelHeight] = ginput;

            %meta.IC1(f) = round(IC1/100,2);
            %meta.IC2(f) = round(IC2/100,2);

            close all

        end % Side Loop
    end % File Loop
    %filename = fullfile(synergies2Dir,subjects{subj},strcat(subjects{subj},"_trials.xlsx"));
    %writetable(meta,filename)
end % Subject Loop

