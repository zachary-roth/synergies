
function [EMG_max_all] = find_EMG_max(EMG_dir,EMG_files)
% Find the maximum Activation values across all trials

% Loop over the EMG trials
for f = 1:length(EMG_files)
    
    % Import the raw EMG files
    raw = readmatrix(fullfile(EMG_dir,EMG_files{f}));
    if isnan(raw)
        continue
    end
        
    % Notch filter: Remove the 50Hz noise
    % Create the notch filter
    fs = 1000;
    fo = 50;
    q = 35;
    bw = (fo/(fs/2))/q;
    [b,a] = iircomb(round(fs/fo),bw,'notch');
    
    % Initialize an output array
    notch = zeros(size(raw));
    
    % Filter each column
    for m = 1:width(raw)
        notch(:,m) = filter(b,a,raw(:,m));
    end

    % Bandpass filter (20-400 Hz)
    % Create the bandpass filter
    order = 4;
    cutoff = [20 400];
    [b,a] = butter(order/2, cutoff/(0.5*fs));
    
    % Intialize an output array
    bandPass = zeros(size(notch));
    
    % Filter each column
    for m = 1:width(notch)
        bandPass(:,m) = filtfilt(b,a,notch(:,m));
    end

    % Rectification
    rect = abs(bandPass);
    
    % Vertically concatenate all the EMG files into a single array
    if f == 1
        rect_concat = rect;
    else
        rect_concat = [rect_concat; rect];
    end
end

EMG_max_all = max(rect_concat);