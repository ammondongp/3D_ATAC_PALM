%       According to Peebles and Hauser (1974), we define the pair correlation function g(r) measures 
% the prob?ability dP of finding an enhancer site in a volume element dV at a separation r from 
% another enhancer site:
%       ¦¤P=ng(r) ¦¤V
% where n is the mean number density of the enhancers in the nucleus. 
% 
%       In practice, the pair correlation function can be estimated from a sample of objects counting 
% the pairs of objects with different separations r [Peebles & Hauser [4] estimator]: 
%       g(r)=(NR*DD(r))/(N*RR(r))
% where DD(r) and RR(r) are counts of pairs of enhancers (in bins of separation) in the data 
% catalog and in the random catalog, respectively. 
%
%       The random catalog consists of uniformly distributed positions in the same volume defined
% by data catalog 3D convex hull. To reduce the noise, we computationally generate the random 
% catalog that has a size 10 times greater than that of the data catalog. The normalizing coefficients
% containing the num?bers of points in the initial (N) and random (Tarantino et al.) catalogs are 
% included in the estimator.
%
%       Here, non-redundant pair wise Euclidean distance set within each catalog can be constructed
% by: 
%       dst(i, j)(i¡Ùj) = ||ri-rj||¡¬
% We define: 
%       C(dst(i, j),r) = 0 (when dst(i, j) > r) or C(dst(i, j),r) = 1 (when dst(i, j) <= r)
% 
%       The bin size of the g(r) distribution function is ¦¤r.
% Then,
%       g(r)=(NR*DD(r))/(N*RR(r))=(NR/N)*(QDD/QRR)
% in which
%       QDD = ¦²(CDD*(dst(i, j), r+¦¤r/2))- ¦²(CDD*(dst(i, j), r-¦¤r/2))
% and
%       QRR = ¦²(CRR*(dst(i, j), r+¦¤r/2))- ¦²(CRR*(dst(i, j), r-¦¤r/2))
%
% Please refer to:
%   Liu et al., 3D imaging of Sox2 enhancer clusters in embryonic stem cells, Elife (2014)
%   --> Materials and methods --> 3D pair correlation function and calculation


%%prefacts3
clear;
clc;

[file,folder]=uigetfile('*');
FileQuant=fullfile(folder, file);
Spots = load(FileQuant);%'cell_05.3d'); % 3d files to calculate autocorrelation function
name = FileQuant(1:end-3);
number_X = 5000; % number of random points to generate
Nsim = 4; % number of simulations for calculating 3D shape autocorrelation
dr = 50; %nm

%% Select Spots within the nucleus
X=Spots(:,1:3);

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

X(:,1)=X(:,1)-min(X(:,1));
X(:,2)=X(:,2)-min(X(:,2));
X(:,3)=X(:,3)-min(X(:,3));

% p = randperm(length(X(:,1)),number);

len_X=length(X(:,1));
Xs=X;
%% Generate shape autocorrelation curve for the hull
width = max([max(X(:,1)),max(X(:,2)),max(X(:,3))]);
[~, vol_s]= boundary(Xs);

% determine the number of localizations to be analyzed after normalizing to a specific density: number_X/(3E11)
PointsNum = round(number_X*vol_s/(3E11));
PointsNum =min(PointsNum, size(X,1));

X=X(1:PointsNum, :);
[k vol]= boundary(X); %conv hull of real data

% calculate the number for uniform sampling 
number = round(PointsNum*width^3/vol);

RandomAll = [];
CorrR = {};
rR = {};
len=zeros(1, Nsim);

% sampling uniform distribution for Nsim rounds
parfor i = 1:Nsim
    random = rand(3, number)'*width;
    InR = inhull(random, X(k,:));
    random = random(InR,:);
    [CorrR{i} rR{i}] = spatialxcorr_3D_without_edge(random, random, dr, vol);
    len(i) = length(CorrR{i});
    RandomAll = vertcat(RandomAll, random);
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
[Corr r] = spatialxcorr_3D_without_edge(X, X, dr, vol);
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
filename = strcat(name, '_autocorrelation.txt');
dlmwrite(filename, data_plot, 'Delimiter','\t');
%% random control
% RCorr = CorrR{2}(1:mlen)./Corr_R;
% figure, plot(Fr, RCorr);
% xlim([0 2000]);