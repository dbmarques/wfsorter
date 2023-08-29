% WFSORTERSTATSFCN() - Calcualtes estimates of sorting quality statistics.
%   Calculations are mostly as performed and describe as in Plexon's 
%   Offline Sorting. Calculates F(MANOVA), J3, Pseudo-F, Davies-Bouldin 
%   index, and Dunn index.
% 
%   Usage:
%       [Fmanova,p,J3,PseudoF,DB,Dunn] = wfsorterstatsfcn(X,group,onlysorted)
% 
%   Inputs:
%       X = multivariate data/spikes (observations, variables)
%       group = index vector of groups/units
%       onlysorted = uses (0) all or (1) ignores group==0. [default: 0]
% 
%   Outputs:
%       Fmanova = F-value from MANOVA
%       p = p-value from MANOVA
%       J3 = J3 index
%       PseudoF = Pseudo-F value
%       DB = Davies-Bouldin index
%       Dunn = Dunn index
%  
%   References:
%       https://plexon.com/wp-content/uploads/2017/06/Offline-Sorter-v4.3-User-Guide.pdf
%           pg. 325-328.
%       Wheeler, Bruce C., Automatic Discrimination of Single Units in
%           Methods for Neural Ensemble Recordings, ed. by Nicolelis, M., 
%           CRC Press, Boca Raton, 1999.
% 
%   Author: Danilo Benette Marques, 2023.

function [Fmanova,pvalue,J3,PseudoF,DB,Dunn] = sortstatsfcn(X,group,onlysorted)

if nargin<3
    onlysorted = 0;
end

if onlysorted %Calculates stats ignoring group==0 ("unsorted")
X = X(group>0,:);
group = group(group>0,:);
end

%Get all estimates
[Fmanova,pvalue] = fmanova(X,group);
[J3,PseudoF] = j3(X,group);
[DB] = daviesbouldin(X,group);
[Dunn] = dunn(X,group);

% %% SORTING QUALITY STATISTICS FUNCTIONS

%--------------------------------------------------------------------------
% FMANOVA() - Calculates F-statistic from multivariate one-way analysis of
%   variance (MANOVA). Designed to evaluate spike sorting quality.
% 
%   Usage:
%       [F,pvalue,df] = fmanova(X,group)
% 
%   Inputs:
%       X = multivariate data (observations, variables)
%       group = index vector of groups
% 
%   Outputs:
%       F = F-statistic
%       pvalue = p-value
%       df = degrees of freedom. df(1) = between groups, df(2) =
%           observations. Report as F(df(1),df(2)) = F.
% 
%   Author: Danilo Benette Marques, 2023.

function [F,pvalue,df] = fmanova(X,group)

% Multivariate one-way analysis of variance
[d,pvalue,stats] = manova1(X,group);
pvalue = pvalue(d); %needs review!

%Wilks' lambda
lambda = stats.lambda(1); %needs confirmation!

df(1) = stats.dfB; %g - 1
df(2) = stats.dfW; %N - g

g = numel(unique(group)); %number of groups
p = size(X,2); %number of dependent variables
N = size(X,1); %total number of observations

%Converts Wilks' lambda to F-statistic
F = WilksLambda2F(lambda,g, p, N);

%Yitian Shao's UCSB script (Rencher, 2003)
function F = WilksLambda2F(lambda, g, p, N)
% Converts Wilks' Lambda to F-value
% 
%   Usage:
%       F = WilksLambda2F(lambda, g, p, N)
% 
%   Inputs:
%       lambda = Wilks' Lambda 
%       g = number of groups (classes) (k)
%       p = number of dependent variables (features) (p)              
%       N =  total number of obseravations (of all classes) (k*n)
%
%   Outputs:
%       F = F-value
%
%   Reference: 
%       Rencher, Alvin C. Methods of multivariate analysis. Vol. 492. 
%           John Wiley & Sons, 2003., p156-163.
%   Original: https://www.mathworks.com/matlabcentral/fileexchange/67208-maxwellre-convert-wilk-lambda-to-f-value
%   Author: Yitian Shao (yitianshao@ucsb.edu). Created on 05/03/2018

vH = g-1; %(g-1)
vE = N-g; %g*(N-1)

% F statistic is exact if min(p,vH) = 1 or 2, otherwise it is approximate.
if (vH == 1) % exact F-value
    S = 1;
    DF_num = p;
    DF_denom = vE -p +1;
elseif (vH == 2) % exact F-value
    S = 2;
    DF_num = 2*p;
    DF_denom = 2*(vE -p +1);
elseif (p == 1) % exact F-value
    S = 1;
    DF_num = vH;
    DF_denom = vE;
elseif (p == 2) % exact F-value
    S = 2;
    DF_num = 2*vH;
    DF_denom = 2*(vE-1);
else % approximate F-value
    S = sqrt( (p^2*vH^2 -4)/(p^2+vH^2-5) );
    DF_num = p*vH;
    DF_denom = S*(vE+vH -0.5*(p+vH+1)) -0.5*(p*vH-2);
end

Y = lambda^(1/S);
F = ((1-Y)/Y)*(DF_denom/DF_num);
end

end

%--------------------------------------------------------------------------
% J3() - Calculates J3 index as performed and described in Plexon's Offline
%   Sorter. Also calculates Pseudo-F. Designed to evaluate spike sorting
%   quality.
% 
%   Usage:
%       [J3,PseudoF] = j3(X,group)
% 
%   Inputs:
%       X = multivariate data (observations, variables)
%       group = index vector of groups/units
% 
%   Outputs:
%       J3 = J3 index
%       PseudoF = Pseudo-F
%  
%   References:
%       https://plexon.com/wp-content/uploads/2017/06/Offline-Sorter-v4.3-User-Guide.pdf
%           pg. 327.
%       Wheeler, Bruce C., Automatic Discrimination of Single Units in
%           Methods for Neural Ensemble Recordings, ed. by Nicolelis, M., 
%           CRC Press, Boca Raton, 1999.
% 
%   Author: Danilo Benette Marques, 2023.

function [J3,PseudoF] = j3(X,group)

units = unique(group);
for u = 1:numel(units)
   f = X(group==units(u),:); %points of feature space of unit u
   mu = mean(X(group==units(u),:),1); %center of unit u
   E1 = sum((f-mu).^2,2); %Euclidean distance squared
   
   %J1  is a measure of the average distance in feature space between points in a cluster (f) from their center (m).
   j1(u) = sum(E1,1);
   
   N = sum(group==units(u)); %number of points in unit u
   m = mean(X,1); %grand center of all points in all units
   E2 = sum((mu-m).^2,2); %Euclidean distance squared
   
   %J2 is a measure of the average distance between unit clusters.
   j2(u) = N*E2;
    
end
J1 = sum(j1);
J2 = sum(j2);

J3 = J2/J1;

% Calculates Pseudo-F
% It is essentially J3 that has been adjusted for the number of waveforms 
% and the number of units
N = size(X,1);
g = numel(unique(group));

PseudoF = (N-g)/(g-1) * J3;

end

%--------------------------------------------------------------------------
% DAVIES-BOULDIN - Calculates Davies-Bouldin index using MATLAB
%   evalclusters built-in functions
function [DB] = daviesbouldin(X,group)
    
   %Substitutes zeros and puts k in order
   clust = zeros(size(group));
   units = unique(group);
   for k = 1:numel(units)
       clust(group==units(k),:) = k;
   end
   
   evalK = evalclusters(X,clust,'daviesbouldin');
   DB = evalK.CriterionValues;
end

%--------------------------------------------------------------------------
% DUNN() - Calculates Dunn index as performed and described in Plexon's
%   Offline Sorter. Designed to evaluate spike sorting
% 
%   Usage:
%       [Dunn] = dunn(X,group)
% 
%   Inputs:
%       X = multivariate data (observations, variables)
%       group = index vector of groups/units
% 
%   Outputs:
%       Dunn = Dunn index
%  
%   References:
%       https://plexon.com/wp-content/uploads/2017/06/Offline-Sorter-v4.3-User-Guide.pdf
%           pg. 328.
% 
%   Author: Danilo Benette Marques, 2023.

function [Dunn] = dunn(X,group)

units = unique(group);

%d(i,j) is the Euclidean distance between the centroids of unit i and unit j
for i = 1:numel(units)
    centroidi = mean(X(group==units(i),:),1);
    
    for j = 1:numel(units)
    centroidj = mean(X(group==units(j),:),1);
    
    dij(i,j) = pdist2(centroidi, centroidj); 
    end
end
dij(logical(eye(size(dij,1))))=NaN; %i==j --> NaN

%d'(k) is the average distance of each point in unit k to the centroid of unit k
for k = 1:numel(units)
    centroidk = mean(X(group==units(k),:),1);
    dk(k) = mean(pdist2(X(group==units(k),:),centroidk),1); %average distance of each point in unit k to the centroid of unit k
end
maxkdk = max(dk);

Dunn = min(min(dij/maxkdk,[],2));
end

end