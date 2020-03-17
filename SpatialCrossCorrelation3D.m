clear;
clc;
disp('Load ATAC data:');
[file,folder]=uigetfile('*');
FileQuant=fullfile(folder, file);
Spots= load(FileQuant);

disp('Load CTCF data:');
[file2,folder2]=uigetfile('*');
FileQuant2=fullfile(folder2, file2);
Spots2= load(FileQuant2);

name = FileQuant(1:end-3);

number_X = 5000; % number of random points to generate
Nsim = 4; % number of simulations for calculating 3D shape autocorrelation
dr = 50; %nm

%% Select Spots within the nucleus

X=Spots(:,1:3);
Y=Spots2(:,1:3);
% %Select xy data
% xy=figure;
% scatter(X(:,1), X(:,2), '.');
% 
% hxy = impoly();
% nodesxy = getPosition(hxy);
% close(xy);
% idx1 = inpoly([X(:,1), X(:,2)], nodesxy);
% 
% X=X(idx1, :);
% 
% 
% %Select xz data
% xz=figure;
% scatter(X(:,1), X(:,3), '.');
% 
% hxz = impoly();
% nodesxz = getPosition(hxz);
% close(xz);
% idx2 = inpoly([X(:,1), X(:,3)], nodesxz);
% X=X(idx2, :);
% 
% 
% %Select yz data
% yz=figure;
% scatter(X(:,2), X(:,3), '.');
% 
% 
% hyz = impoly();
% nodesyz = getPosition(hyz);
% close(yz);
% idx3 = inpoly([X(:,2), X(:,3)], nodesyz);
% X=X(idx3, :);
% 
% %saving data 
% Spots = Spots(idx1,:);
% Spots = Spots(idx2,:);
% Spots = Spots(idx3,:);
% % filename = strcat(name, '_trimed.3d');
% % dlmwrite(filename, Spots, 'Delimiter','\t');
xx1=min(min(X(:,1)),min(Y(:,1)));
xx2=min(min(X(:,2)),min(Y(:,2)));
xx3=min(min(X(:,3)),min(Y(:,3)));

X(:,1)=X(:,1)-xx1;
X(:,2)=X(:,2)-xx2;
X(:,3)=X(:,3)-xx3;

Y(:,1)=Y(:,1)-xx1;
Y(:,2)=Y(:,2)-xx2;
Y(:,3)=Y(:,3)-xx3;

% p = randperm(length(X(:,1)),number);

len_X=length(X(:,1));
Xs=X;
%% Generate shape autocorrelation curve for the hull
width = max([max(X(:,1)),max(X(:,2)),max(X(:,3)),max(Y(:,1)),max(Y(:,2)),max(Y(:,3))]);
[~, vol_s]= boundary(Xs);

% determine the number of localizations to be analyzed after normalizing to a specific density: number_X/(3E11)
PointsNum = round(number_X*vol_s/(3E10));

% X=X(1:PointsNum, :);
[k vol]= boundary(X); %conv hull of real data

% calculate the number for uniform sampling 
number = round(PointsNum*width^3/vol);

% RandomAll = [];
CorrR = {};
rR = {};
len=zeros(1, Nsim);

% sampling uniform distribution for Nsim rounds
parfor i = 1:Nsim
    random1 = rand(3, number)'*width;
    random2 = rand(3, number)'*width;
    InR1 = inhull(random1, X(k,:));
    InR2 = inhull(random2, X(k,:));
    random1 = random1(InR1,:);
    random2 = random2(InR2,:);
    [CorrR{i} rR{i}] = spatialxcorr_3D_without_edge(random1, random2, dr, vol);
    len(i) = length(CorrR{i});
%     RandomAll = vertcat(RandomAll, random);
end

mlen=min(len);
Corr_R = zeros(1, mlen);
r_R = zeros(1, mlen);
for i = 1:Nsim
    Corr_R = Corr_R + CorrR{i}(1:mlen);
    r_R = r_R + rR{i}(1:mlen);
end
Corr_R = Corr_R/Nsim;
r_R = r_R/Nsim;

%% Generate shape autocorrelation curve for the real data
[Corr r] = spatialxcorr_3D_without_edge(X, Y, dr, vol);
lenC = length(Corr);
if lenC>=mlen
    Corr = Corr(1:mlen);
    FCorr = Corr./Corr_R;
    Fr = r_R;
else
    Corr_R = Corr_R(1:lenC);
    FCorr = Corr./Corr_R;
    Fr = r;
end

%% Plotting data
figure, plot(Fr, FCorr, '-xg');
hold on;
plot(Fr, Corr, '-xb');
plot(Fr, Corr_R, '-xr');

xlim([0 1000]);

data_plot=[Fr; FCorr; Corr_R]';
filename = strcat(name, '_CrossCorrelation.txt');
dlmwrite(filename, data_plot, 'Delimiter','\t');
%% random control
% RCorr = CorrR{2}(1:mlen)./Corr_R;
% figure, plot(Fr, RCorr);
% xlim([0 2000]);