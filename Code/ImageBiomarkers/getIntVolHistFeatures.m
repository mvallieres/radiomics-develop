function intVolHist = getIntVolHistFeatures(vol,wb,userSetRange)
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

% - vol: 3D volume, QUANTIZED, with NaNs outside the region of interest
% 1) Naturally discretised volume can be kept as is (e.g. HU values of CT scans)
% 2) All other volumes with continuous intensity distribution should be
% quantized (e.g., nBins = 100), with levels = [min, ..., max]

% -> Third argument is optional. It needs to be used when a priori
% discretising with "FBS" or "FBSequal".


% INITIALIZATION
X = vol(~isnan(vol(:)));
if nargin == 3
    if ~isempty(userSetRange)
        minVal = userSetRange(1);
        maxVal = userSetRange(2);
    else
        minVal = min(X);
        maxVal = max(X);
    end
else
    minVal = min(X);
    maxVal = max(X);
end
if maxVal == Inf
    maxVal = max(X);
end
if minVal == -Inf
    minVal = min(X);
end
levels = (minVal:wb:maxVal)'; % Vector of grey-levels
Ng = numel(levels);
Nv = numel(X);

% Calculating fractional volume
fractVol = zeros(Ng,1);
for i = 1:Ng
    fractVol(i) = 1 - sum(X < levels(i))/Nv;
end

% Calculating intensity fraction
fractInt = (levels - min(levels))/(max(levels) - min(levels));

% Volume at intensity fraction 10 
V10 = findVX(fractInt,fractVol,10);
intVolHist.Fivh_V10 = V10;

% Volume at intensity fraction 90
V90 = findVX(fractInt,fractVol,90);
intVolHist.Fivh_V90 = V90;

% Intensity at volume fraction 10 --> For initial arbitrary intensities, we will always be discretising (1000 bins). So intensities are definite here.
I10 = findIX(levels,fractVol,10);
intVolHist.Fivh_I10 = I10;

% Intensity at volume fraction 90 --> For initial arbitrary intensities, we will always be discretising (1000 bins). So intensities are definite here.
I90 = findIX(levels,fractVol,90);
intVolHist.Fivh_I90 = I90; 

% Volume at intensity fraction difference V10-V90
intVolHist.Fivh_V10minusV90 = V10 - V90;

% Intensity at volume fraction difference I10-I90 --> For initial arbitrary intensities, we will always be discretising (1000 bins). So intensities are definite here.
intVolHist.Fivh_I10minusI90 = I10 - I90;

% Area under IVH curve
intVolHist.Fivh_auc = trapz(fractVol)/(Ng - 1);

end