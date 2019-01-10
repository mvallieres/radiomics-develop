function [GLDZM] = getGLDZMmatrix(ROIOnlyInt,mask,levels)
% -------------------------------------------------------------------------
% AUTHOR(S): 
% - Martin Vallieres <mart.vallieres@gmail.com>
% -------------------------------------------------------------------------
% HISTORY:
% - Creation: April 2017
% - Revision I: August 2017
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

% - vol: 3D volume, isotropically resampled, quantized, with NaNs outside the region of interest
% - levels: should be removed at some point, no longer needed if we always
%   quantize our volume such that levels = 1,2,3,4,...,max(quantized Volume). 
%   So simply calculate levels = 1:max(ROIOnly(~isnan(ROIOnly(:))))
%   directly in this function


% COMPUTATION OF DISTANCE MAP
mask = padarray(mask,[1,1,1],0);
perimeter = bwperim(mask,6); % Computing the smallest ROI edge possible. Source of difference?
perimeter = perimeter(2:end-1,2:end-1,2:end-1); % Removing the padding.
mask = mask(2:end-1,2:end-1,2:end-1); % Removing the padding
distMap = bwdist(perimeter,'cityblock') + 1; % +1 according to the definition of the IBSI


% INITIALIZATION
Ng = length(levels); % Since levels is always defined as 1,2,3,4,...,max(quantized Volume)
levelTemp = max(levels) + 1;
ROIOnlyInt(isnan(ROIOnlyInt)) = levelTemp;
distInit = max(distMap(mask == 1)); % Since the ROI morph always encompasses ROI int, using the mask as defined from ROI morph does not matter since we want to find the maximal possible distance.
GLDZM = zeros(Ng,distInit);


% COMPUTATION OF GLDZM
temp = ROIOnlyInt;
for i = 1:Ng
    temp(ROIOnlyInt~=levels(i)) = 0;
    temp(ROIOnlyInt==levels(i)) = 1;
    connObjects = bwconncomp(temp,26);
    nZone = length(connObjects.PixelIdxList);
    for j = 1:nZone
        col = min(distMap(connObjects.PixelIdxList{j}));
        GLDZM(i,col) = GLDZM(i,col) + 1;
    end
end


% REMOVE UNECESSARY COLUMNS
stop = find(sum(GLDZM),1,'last');
GLDZM(:,(stop+1):end) = [];

end