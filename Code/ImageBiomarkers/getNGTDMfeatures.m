function ngtdm = getNGTDMfeatures(vol,distCorrection)
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

% - vol: 3D volume, isotropically resampled, quantized (e.g. Ng = 32, levels = [1, ..., Ng]), with NaNs outside the region of interest
% - distCorrection: % Set this variable to true in order to use discretization length difference corrections as used here: https://doi.org/10.1088/0031-9155/60/14/5471. Set this variable to false to replicate IBSI results.


% GET THE NGTDM MATRIX
levels = 1:max(vol(~isnan(vol(:)))); % Correct definition, without any assumption
if nargin == 2
    [NGTDM,countValid] = getNGTDMmatrix(vol,levels,distCorrection);
else
    [NGTDM,countValid] = getNGTDMmatrix(vol,levels);
end
nTot = sum(countValid);
countValid = countValid./nTot; % Now representing the probability of gray-level occurences
NL = length(NGTDM);
Ng = sum(countValid~=0);
pValid = find(countValid>0);
nValid = length(pValid);



% COMPUTING TEXTURES

% Coarseness
coarseness = ((countValid')*NGTDM)^(-1);
if coarseness > 10^6
    coarseness = 10^6;
end
ngtdm.Fngt_coarseness = coarseness;

% Contrast
if Ng == 1
    ngtdm.Fngt_contrast = 0;
else
    val = 0;
    for i = 1:NL
        for j = 1:NL
            val = val + countValid(i)*countValid(j)*(i-j)^2;
        end
    end
    ngtdm.Fngt_contrast = val*sum(NGTDM)/(Ng*(Ng-1)*nTot);
end

% Busyness
if Ng == 1
    ngtdm.Fngt_busyness = 0;
else
    denom = 0;
    for i = 1:nValid
        for j = 1:nValid
            denom = denom + abs(pValid(i)*countValid(pValid(i))-pValid(j)*countValid(pValid(j)));
        end
    end
    ngtdm.Fngt_busyness = ((countValid')*NGTDM)/denom;
end

% Complexity
val = 0;
for i = 1:nValid
    for j = 1:nValid
        val = val + (abs(pValid(i)-pValid(j))/(nTot*(countValid(pValid(i)) + countValid(pValid(j)))))*(countValid(pValid(i))*NGTDM(pValid(i)) + countValid(pValid(j))*NGTDM(pValid(j)));
    end
end
ngtdm.Fngt_complexity = val;

% Strength
if sum(NGTDM) == 0
    ngtdm.Fngt_strength = 0;
else
    val = 0;
    for i = 1:nValid
        for j = 1:nValid
            val = val + (countValid(pValid(i))+countValid(pValid(j)))*(pValid(i)-pValid(j))^2;
        end
    end
    ngtdm.Fngt_strength = val/sum(NGTDM);
end

end