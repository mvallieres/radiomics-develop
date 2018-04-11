function stats = getStatsFeatures(vol,intensity)
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
% - intensity (optional): If 'arbitrary', some feature will not be computed. If
%   'definite', all feature will be computed. If not present as an argument,
%   all features will be computed. Here, 'filter' is the same as
%   'arbitrary'.


% PRELIMINARY
if nargin < 2
    definite = true;
else
    if strcmp(intensity,'arbitrary')
        definite = false;
    elseif strcmp(intensity,'definite')
        definite = true;
    elseif strcmp(intensity,'filter')
        definite = false;
    else
        error('Second argument must either be "arbitrary" or "definite" or "filter"')
    end 
end

% INITIALIZATION
X = vol(~isnan(vol(:)));
Nv = numel(X); 


% STARTING COMPUTATION

% Mean
if definite
    u = sum(X)/Nv;
    stats.Fstat_mean = u;
else
    stats.Fstat_mean = [];
end

% Variance
if definite
    var = sum((X - u).^2)/Nv;
    stats.Fstat_var = var;
else
    stats.Fstat_var = [];
end

% Skewness
if definite
    skew = 0;
    if var ~= 0
        skew = (sum((X - u).^3)/Nv) / var^(3/2);
    end
    stats.Fstat_skew = skew;
else
    stats.Fstat_skew = [];
end

% Kurtosis
if definite
    kurt = 0;
    if var ~= 0
        kurt = (sum((X - u).^4)/Nv) / var^2 - 3;
    end
    stats.Fstat_kurt = kurt;
else
    stats.Fstat_kurt = [];
end

% Median
if definite
    med = median(X);
    stats.Fstat_median = med;
else
    stats.Fstat_median = [];
end

% Minimum grey level
if definite
    stats.Fstat_min = min(X);
else
    stats.Fstat_min = [];
end

% 10th percentile
if definite
    p10 = calcPercentile(X,10);
    stats.Fstat_P10 = p10;
else
    stats.Fstat_P10 = [];
end

% 90th percentile
if definite
    p90 = calcPercentile(X,90);
    stats.Fstat_P90 = p90;
else
    stats.Fstat_P90 = [];
end

% Maximum grey level
if definite
    stats.Fstat_max = max(X);
else
    stats.Fstat_max = [];
end

% Interquantile range
if definite
    p75 = calcPercentile(X,75); p25 = calcPercentile(X,25);
    stats.Fstat_iqr = p75 - p25;
else
    stats.Fstat_iqr = [];
end

% Range
if definite
    stats.Fstat_range = max(X) - min(X);
else
    stats.Fstat_range = [];
end

% Mean absolute deviation
if definite
    stats.Fstat_mad = sum(abs(X - u))/Nv;
else
    stats.Fstat_mad = [];
end

% Robust mean absolute deviation
if definite
    X_10_90 = X(X >= p10 & X <= p90);
    N_10_90 = numel(X_10_90);
    stats.Fstat_rmad = sum(abs(X_10_90 - mean(X_10_90)))/N_10_90;
else
    stats.Fstat_rmad = [];
end

% Median absolute deviation
if definite
    stats.Fstat_medad = sum(abs(X - med))/Nv;
else
    stats.Fstat_medad = [];
end

% Coefficient of variation
if definite
    stats.Fstat_cov = sqrt(var)/u;
else
    stats.Fstat_cov = [];
end

% Quartile coefficient of dispersion
if definite
    stats.Fstat_qcod = (p75 - p25)/(p75 + p25);
else
    stats.Fstat_qcod = [];
end

% Energy
if definite
    energy = sum(X.^2);
    stats.Fstat_energy = energy;
else
    stats.Fstat_energy = [];
end

% Root mean square
if definite
    stats.Fstat_rms = sqrt(energy/Nv);
else
    stats.Fstat_rms = [];
end

end