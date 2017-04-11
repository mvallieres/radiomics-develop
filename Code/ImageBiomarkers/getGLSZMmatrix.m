function [GLSZM] = getGLSZMmatrix(ROIOnly,levels)
% -------------------------------------------------------------------------
% [GLSZM] = getGLSZM(ROIOnly,levels)
% -------------------------------------------------------------------------
% DESCRIPTION:
% This function computes the Gray-Level Size Zone Matrix (GLSZM) of the 
% region of interest (ROI) of an input volume. The input volume is assumed 
% to be isotropically resampled. The zones of different sizes are computed 
% using 26-voxel connectivity.
%
% --> This function is compatible with 2D analysis (language not adapted in the text)
% -------------------------------------------------------------------------
% REFERENCE:
% [1] Thibault, G., Fertil, B., Navarro, C., Pereira, S., Cau, P., Levy, 
%     N., Mari, J.-L. (2009). Texture Indexes and Gray Level Size Zone 
%     Matrix. Application to Cell Nuclei Classification. In Pattern 
%     Recognition and Information Processing (PRIP) (pp. 140–145).
% -------------------------------------------------------------------------
% INPUTS:
% - ROIonly: Smallest box containing the ROI, with the imaging data ready 
%            for texture analysis computations. Voxels outside the ROI are 
%            set to NaNs.
% - levels: Vector containing the quantized gray-levels in the tumor region
%           (or reconstruction levels of quantization).
%
% ** 'ROIonly' and 'levels' should be outputs from 'prepareVolume.m' **
% -------------------------------------------------------------------------
% OUTPUTS:
% - GLSZM: Gray-Level Size Zone Matrix of 'ROIOnly'.
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


% % PRELIMINARY
% nLevel = length(levels);
% if nLevel > 100
%     adjust = 10000;
% else
%     adjust = 1000;
% end
levelTemp = max(levels) + 1;
ROIOnly(isnan(ROIOnly)) = levelTemp;
levels = [levels,levelTemp];


% % QUANTIZATION EFFECTS CORRECTION
% % In case (for example) we initially wanted to have 64 levels, but due to
% % quantization, only 60 resulted.
% uniqueVect = round(levels*adjust)/adjust;
% ROIOnly = round(ROIOnly*adjust)/adjust;
uniqueVect = levels;
NL = length(levels) - 1;


% INITIALIZATION
nInit = numel(ROIOnly(:)); % THIS NEEDS TO BE CHANGED. THE ARRAY INITIALIZED COULD BE TOO BIG!
GLSZM = zeros(NL,nInit);


% COMPUTATION OF GLSZM
temp = ROIOnly;
for i = 1:NL
    temp(ROIOnly~=uniqueVect(i)) = 0;
    temp(ROIOnly==uniqueVect(i)) = 1;
    connObjects = bwconncomp(temp,26);
    nZone = length(connObjects.PixelIdxList);
    for j = 1:nZone
        col = length(connObjects.PixelIdxList{j});
        GLSZM(i,col) = GLSZM(i,col) + 1;
    end
end


% REMOVE UNECESSARY COLUMNS
stop = find(sum(GLSZM),1,'last');
GLSZM(:,(stop+1):end) = [];

end