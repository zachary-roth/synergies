function [ursqr, VAF] = rsqr_uncentered_VAF(data,data_rec)
% RSQRUncentered: This function calculates the uncetered correlation coefficient between data 
% and data_rec (both row vectors)
%
% Syntax:   r_sqr = rsqr_uncentered(data,data_rec)
%
% Input:
% data      Array   matrix of observed data  (e.g., data = [mus pert_dir])
% data_rec  Array   matrix of reconstructed/predicted data (e.g., data_rec = [mus pert_dir])
%
% Output:
% ursqr     Array   matrix with uncentered correlation coefficients
%
% Calls:
% std_mean0.m 
%       
% Created: May 24, 2006 (Gelsy Torres-Oviedo)
% Last Modified: July 10, 2006 (Torrence Welch)
% Last Modification: fix ursqr calculation
%*********************************************************************************************************

% Shift dimensions for regression purposes
%     each data channel in separate columns
warning off
data = data';
data_rec = data_rec';

% % Cluster 3.0 method
% dim_data = size(data);
% for i = 1:dim_data(2)
%     % Calculate standard deviation assuming mean of 0
%     datastd0 = std_mean0(data(:,i));
%     datarecstd0 = std_mean0(data_rec(:,i));
% 
%     X = [data(:,i)/datastd0 data_rec(:,i)/datarecstd0];
%     n = length(X);
%     ur(i) = sum(prod(X,2))/n;
%     ursqr(i) = ur(i)^2;
% end

% Zar book method
dim_data = size(data);
for i = 1:dim_data(2)
    X = [data(:,i) data_rec(:,i)];
    % prod(X,2) is product of all elements in row
    ursqr(i) = sum(prod(X,2))^2 / (sum(data(:,i).^2)*sum(data_rec(:,i).^2));
    VAF(i) = 1 - sum((data(:,i)-data_rec(:,i)).^2)/sum(data(:,i).^2);
end

ursqr = ursqr';
VAF = VAF';
return
%========================
% end rsqr_uncentered.m