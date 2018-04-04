function intHist = getIntHistFeatures(vol)
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

% - vol: 3D volume, QUANTIZED (e.g. nBins = 100, levels = [1, ..., max]), with NaNs outside the region of interest


% INITIALIZATION
X = vol(~isnan(vol(:)));
Nv = numel(X);

% CONSTRUCTION OF HISTOGRAM AND ASSOCIATED NUMBER OF GRAY-LEVELS
levels = 1:max(X); % Always defined from 1 to the maximum value of the volume to remove any ambiguity
Ng = numel(levels); % Number of gray-levels
H = zeros(1,Ng); % The histogram of X
for i = 1:Ng
    H(i) = sum(X == i); % == i or == levels(i) is equivalent since levels = 1:max(X), and Ng = numel(levels)
end
p = H./Nv; % Occurence probability for each grey level bin i


% STARTING COMPUTATION

% Intensity histogram mean
u = levels * p';
intHist.Fih_mean = u;

% Intensity histogram variance
var = ((levels -u).^2) * p';
intHist.Fih_var = var;

% Intensity histogram skewness
skew = 0;
if var ~= 0
    skew = (((levels -u).^3) * p') / var^(3/2);
end
intHist.Fih_skew = skew;

% Intensity histogram kurtosis
kurt = 0;
if var ~= 0
    kurt = (((levels -u).^4) * p') / var^2 - 3;
end
intHist.Fih_kurt = kurt;

% Intensity histogram median
med = median(X);
intHist.Fih_median = med;

% Intensity histogram minimum grey level
intHist.Fih_min = min(X);

% Intensity histogram 10th percentile
p10 = calcPercentile(X,10);
intHist.Fih_P10 = p10;

% Intensity histogram 90th percentile
p90 = calcPercentile(X,90);
intHist.Fih_P90 = p90;

% Intensity histogram maximum grey level
intHist.Fih_max = max(X);

% Intensity histogram mode
mode = find(H == max(H)); % levels = 1:max(X), so the index of the ith bin of H is the same as i
if numel(mode) > 1
    dist = abs(mode - u);
    [~,indMin] = min(dist);
    mode = mode(indMin);
end
intHist.Fih_mode = mode;

% Intensity histogram interquantile range
p75 = calcPercentile(X,75); p25 = calcPercentile(X,25);
intHist.Fih_iqr = p75 - p25; % Since X goes from 1:max(X), all with integer values, the result is an integer

% Intensity histogram range
intHist.Fih_range = max(X) - min(X);

% Intensity histogram mean absolute deviation
intHist.Fih_mad = sum(abs(X - u))/Nv;

% Intensity histogram robust mean absolute deviation
X_10_90 = X(X >= p10 & X <= p90);
N_10_90 = numel(X_10_90);
intHist.Fih_rmad = sum(abs(X_10_90 - mean(X_10_90)))/N_10_90;

% Intensity histogram median absolute deviation
intHist.Fih_medmad = sum(abs(X - med))/Nv;

% Intensity histogram coefficient of variation
intHist.Fih_cov = sqrt(var)/u;

% Intensity histogram quartile coefficient of dispersion
intHist.Fih_qcod = (p75 - p25)/(p75 + p25);

% Intensity histogram entropy
p = p(p > 0);
intHist.Fih_entropy = -sum(p.*log2(p));

% Intensity histogram uniformity
intHist.Fih_uniformity = sum(p.^2);

% Calculation of histogram gradient
histGrad = zeros(1,Ng); histGrad(1) = H(2) - H(1); histGrad(end) = H(end) - H(end-1);
for i = 2:(Ng-1)
    histGrad(i) = (H(i+1) - H(i-1))/2;
end

% Maximum histogram gradient
[intHist.Fih_max_grad,indMax] = max(histGrad);

% Maximum histogram gradient grey level
intHist.Fih_max_grad_gl = levels(indMax);

% Minimum histogram gradient
[intHist.Fih_min_grad,indMin] = min(histGrad);

% Minimum histogram gradient grey level
intHist.Fih_min_grad_gl = levels(indMin);

end