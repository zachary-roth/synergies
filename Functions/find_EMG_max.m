
function [maxEMG_array,maxEMG_struc] = find_EMG_max(subjDir)
% Find the maximum Activation values across all trials
% Inputs
%
%

% Outputs
% maxEMG_array = a 1 x muscle array of the max EMG values 
% maxEMG_struct = a structure containing the rectified EMG signals for each
%                 trial

% Get a list of EMG files
EMG_dir = fullfile(subjDir,'EMG');
EMG_files = dir(fullfile(EMG_dir, '*.xlsx'));
EMG_files = {EMG_files.name};

% Loop over the EMG trials
for f = 1:length(EMG_files)
    
    % Import the raw EMG files
    raw = readmatrix(fullfile(EMG_dir,EMG_files{f}));
    if isnan(raw)
        continue
    end
        
    % Demean
    demean = raw - mean(raw,1);
    
    % Notch filter: Remove the 50Hz noise
    % Create the notch filter
    fs = 1000;
    fo = 50;
    q = 35;
    bw = (fo/(fs/2))/q;
    [b,a] = iircomb(round(fs/fo),bw,'notch');
    
    % Initialize an output array
    notch = zeros(size(demean));
    
    % Filter each column
    for musc = 1:width(demean)
        notch(:,musc) = filter(b,a,demean(:,musc));
    end

    % Bandpass filter (20-400 Hz)
    % Create the bandpass filter
    order = 4;
    cutoff = [20 400];
    [b,a] = butter(order/2, cutoff/(0.5*fs));
    
    % Intialize an output array
    bandPass = zeros(size(notch));
    
    % Filter each column
    for musc = 1:width(notch)
        bandPass(:,musc) = filtfilt(b,a,notch(:,musc));
    end

    % Rectification
    rect = abs(bandPass);
    
    trialName = strrep(EMG_files{f},".xlsx","");
    maxEMG_struc.(trialName) = rect;

    % Vertically concatenate all the EMG files into a single array
    if f == 1
        rect_concat = rect;
    else
        rect_concat = [rect_concat; rect];
    end
end

maxEMG_array = max(rect_concat);


