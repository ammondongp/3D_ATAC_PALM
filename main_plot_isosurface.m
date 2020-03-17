%
% Copyright (c) 2015, Yarpiz (www.yarpiz.com)
% All rights reserved. Please read the "license.txt" for license terms.
%
% Project Code: YPML110
% Project Title: Implementation of DBSCAN Clustering in MATLAB
% Publisher: Yarpiz (www.yarpiz.com)
% 
% Developer: S. Mostapha Kalami Heris (Member of Yarpiz Team)
% 
% Contact Info: sm.kalami@gmail.com, info@yarpiz.com
%

clc;
clear;
close all;

[file,folder]=uigetfile('*');
FileQuant=fullfile(folder, file);

data = load(FileQuant);
dataX=zeros(length(data(:,1)),3);

scale_xyz=100;

x_min=min(data(:,1));
y_min=min(data(:,2));
z_min=min(data(:,3));

ref=[scale_xyz, x_min, y_min, z_min];

X0 = ceil((max(data(:,1)) - x_min)/scale_xyz);
Y0 = ceil((max(data(:,2)) - y_min)/scale_xyz);
Z0 = ceil((max(data(:,3)) - z_min)/scale_xyz);

dataX(:,1) = (data(:,1) - x_min)./scale_xyz;
dataX(:,2) = (data(:,2) - y_min)./scale_xyz;
dataX(:,3) = (data(:,3) - z_min)./scale_xyz;

% for normalization
[kX vol]= convhulln(dataX);


%% Run DBSCAN Clustering Algorithm

epsilon=1.5; %default 1.5
MinPts=10; % default 10
IDX=DBSCAN(dataX,epsilon,MinPts);

IDX_nonzero_length = length(nonzeros(IDX));

%% Plot 2D projection of clusters

PlotClusterinResult(dataX, IDX);
title(['DBSCAN Clustering (\epsilon = ' num2str(epsilon) ', MinPts = ' num2str(MinPts) ')']);

filename = strcat(FileQuant(1: end-4), '_IDX.txt');
dlmwrite(filename, IDX, 'Delimiter','\t');

%% Plot isosurface of each cluster and Record the estimated radius of each cluster
[cluster_identified, vol, vol_nucleus, dens, dens_nucleus] = Plot_clustering_isosurface(dataX, IDX);

radi=nthroot(vol*3./(4*pi),3)*100;
% vol=vol(vol>0);
% 
% radius=zeros(length(vol),1); % unit: micrometer
% for i=1:length(vol)
%     radius(i,1)=nthroot(vol(i)*3/(4*pi), 3)*0.1;
% end
% filename = strcat(FileQuant, '_cluster_radius_record.txt');
% dlmwrite(filename, radius, 'Delimiter','\t');
% 
% % filename1 = strcat(FileQuant, '_cluster_density_record.txt');
% % dlmwrite(filename1, length(vol)/(vol_nucleus/1000), 'Delimiter','\t');

%% Record Sampled dots within each cluster for ViSP visualization

% Fv_sampling = cluster_dots_record(dataX, IDX, ref);
% filename3 = strcat(FileQuant, '_cluster_dots_record.txt');
% dlmwrite(filename3, Fv_sampling, 'Delimiter','\t');