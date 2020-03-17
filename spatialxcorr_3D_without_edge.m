%% 3D spatial cross correlation function
%x, y: 3D array
%dr: bin

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


function [Corr r Density_y width] = spatialxcorr_3D_without_edge(x, y, dr, vol)
% pair-wise distance
c = pdist2(x, y);

% max width
width = max([max(y(:,1)) - min(y(:,1)),max(y(:,2)) - min(y(:,2)),max(y(:,3)) - min(y(:,3))]);

Density_y = size(y,1)/vol;

[m n]= size(x);

steps = floor(2000/dr);

for i = 1:steps,
    r(i) = dr*i - dr/2;
    logic = ((c <= r(i)+ dr/2) & (c > r(i) - dr/2));
    Corr(i) = nnz(logic)/(size(y,1)-1)/(Density_y*pi*(4*(r(i)^2)*dr+1/3*dr^3));
end

end