
close all
clear
clc

% Import Data

[calcFile,calcPath] = uigetfile('*.mat');
calcStruc = load(fullfile(calcPath,calcFile));

% Tranpose the calculated activations to get a time x muscle matrix
calc = calcStruc.Results.MActivation.genericMRS';

% Resample the calculated activations to 0-100% of the gait cycle

% Initialize an output array
calc_resample = zeros(101,width(calc));

for m = 1:width(calc) 
    v = calc(:,m); % sample values
    x = (1:1:length(calc))'; % sample points
    xq = linspace(1,length(calc),101); % query point
    calc_resample(:,m) = interp1(x,v,xq); % resample the data
end

results = struct;

%% NNMF

max_n_synergies = 6;
A = calc_resample;

for k = 1:max_n_synergies
    [W,H,D] = nnmf(A,k,"algorithm","mult","replicates",1000,'Options',statset('Display','final','MaxIter',100));
    results.nnmf.(strcat('k',num2str(k))).W = W;
    results.nnmf.(strcat('k',num2str(k))).H = H;
    results.nnmf.(strcat('k',num2str(k))).D = D;
    results.nnmf.(strcat('k',num2str(k))).recon = W*H;
end

%% VAF
dim_data = size(calc_resample);
for k = 1:max_n_synergies
    data_rec = results.nnmf.(strcat('k',num2str(k))).recon;
    ursqr = zeros(1,dim_data(2));
    VAF = zeros(1,dim_data(2));
    for i = 1:dim_data(2)
        X = [calc_resample(:,i) data_rec(:,i)];
        ursqr(i) = sum(prod(X,2))^2 / (sum(calc_resample(:,i).^2)*sum(data_rec(:,i).^2));
        VAF(i) = 1 - sum((calc_resample(:,i)-data_rec(:,i)).^2)/sum(calc_resample(:,i).^2);
    end
    results.nnmf.ursqr(k,:) = ursqr;
    results.nnmf.VAF_muscles(k,:) = VAF;
    results.nnmf.VAF_trial(k,:) = mean(VAF);
end


%% Minimum Synergies
    results.nnmf.minSynergies = 999;
    for s = 1:max_n_synergies
        if ismember(1,results.nnmf.VAF_muscles(s,:) < .7) == 0 && results.nnmf.VAF_trial(k) > .9 && s < results.nnmf.minSynergies
            results.nnmf.minSynergies = s;
        else
            continue
        end
    end
    results.nnmf.min_synergies_summary = results.nnmf.minSynergies;