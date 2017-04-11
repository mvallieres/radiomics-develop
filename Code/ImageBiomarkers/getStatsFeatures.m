function stats = getStatsFeatures(vol)
% -------------------------------------------------------------------------
% AUTHOR(S): 
% - Martin Vallieres <mart.vallieres@gmail.com>
% -------------------------------------------------------------------------
% HISTORY:
% - Creation: April 2017
% -------------------------------------------------------------------------
% DISCLAIMER:
% "I'm not a programmer, I'm just a scientist doing stuff!"
% -------------------------------------------------------------------------
% STATEMENT:
% This file is part of <https://github.com/mvallieres/radiomics-develop/>, 
% a private repository dedicated to the development of programming code for
% new radiomics applications.
% --> Copyright (C) 2017  Martin Vallieres
%     All rights reserved.
%
% This file is written on the basis of a scientific collaboration for the 
% "radiomics-develop" team.
%
% By using this file, all members of the team acknowledge that it is to be 
% kept private until public release. Other scientists willing to join the 
% "radiomics-develop" team is however highly encouraged. Please contact 
% Martin Vallieres for this matter.
% -------------------------------------------------------------------------

% - vol: 3D volume, NON-QUANTIZED, with NaNs outside the region of interest
% --> vol: continous imaging intensity distribution


% INITIALIZATION
X = vol(~isnan(vol(:)));
Nv = numel(X); 


% STARTING COMPUTATION

% Mean
u = sum(X)/Nv;
stats.Fstat_mean = u;

% Variance
var = sum((X - u).^2)/Nv;
stats.Fstat_var = var;

% Skewness
skew = 0;
if var ~= 0
    skew = (sum((X - u).^3)/Nv) / var^(3/2);
end
stats.Fstat_skew = skew;

% Kurtosis
kurt = 0;
if var ~= 0
    kurt = (sum((X - u).^4)/Nv) / var^2 - 3;
end
stats.Fstat_kurt = kurt;

% Median
med = median(X);
stats.Fstat_median = med;

% Minimum grey level
stats.Fstat_min = min(X);

% 10th percentile
p10 = calcPercentile(X,10);
stats.Fstat_P10 = p10;

% 90th percentile
p90 = calcPercentile(X,90);
stats.Fstat_P90 = p90;

% Maximum grey level
stats.Fstat_max = max(X);

% Interquantile range
p75 = calcPercentile(X,75); p25 = calcPercentile(X,25);
stats.Fstat_iqr = p75 - p25;

% Range
stats.Fstat_range = max(X) - min(X);

% Mean absolute deviation
stats.Fstat_mad = sum(abs(X - u))/Nv;

% Robust mean absolute deviation
X_10_90 = X(X >= p10 & X <= p90);
N_10_90 = numel(X_10_90);
stats.Fstat_rmad = sum(abs(X_10_90 - mean(X_10_90)))/N_10_90;

% Median absolute deviation
stats.Fstat_medmad = sum(abs(X - med))/Nv;

% Coefficient of variation
stats.Fstat_cov = sqrt(var)/u;

% Quartile coefficient of dispersion
stats.Fstat_qcod = (p75 - p25)/(p75 + p25);

% Energy
energy = sum(X.^2);
stats.Fstat_energy = energy;

% Root mean square
stats.Fstat_rms = sqrt(energy/Nv);

end