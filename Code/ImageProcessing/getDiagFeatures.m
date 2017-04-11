function diag = getDiagFeatures(volObj,roiObj,type)
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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FOR THE IMAGE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~strcmp(type,'reSeg')

    % Image dimension x
    diag.(['image_',type,'_dimX']) = volObj.spatialRef.ImageSize(2);

    % Image dimension y
    diag.(['image_',type,'_dimY']) = volObj.spatialRef.ImageSize(1);

    % Image dimension z
    diag.(['image_',type,'_dimZ']) = volObj.spatialRef.ImageSize(3);

    % Voxel dimension x
    diag.(['image_',type,'_voxDimX']) = volObj.spatialRef.PixelExtentInWorldX;

    % Voxel dimension y
    diag.(['image_',type,'_voxDimY']) = volObj.spatialRef.PixelExtentInWorldY;

    % Voxel dimension z
    diag.(['image_',type,'_voxDimZ']) = volObj.spatialRef.PixelExtentInWorldZ;

    % Mean intensity
    diag.(['image_',type,'_meanInt']) = mean(volObj.data(:));

    % Minimum intensity
    diag.(['image_',type,'_minInt']) = min(volObj.data(:));

    % Maximum intensity
    diag.(['image_',type,'_maxInt']) = max(volObj.data(:));

end
% -------------------------------------------------------------------------



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FOR THE ROI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[boxBound] = computeBoundingBox(roiObj.data); 
Xgl = volObj.data(:); Xgl = Xgl(roiObj.data(:) == 1);

% Map dimension x
diag.(['roi_',type,'_dimX']) = roiObj.spatialRef.ImageSize(2);

% Map dimension y
diag.(['roi_',type,'_dimY']) = roiObj.spatialRef.ImageSize(1);

% Map dimension z
diag.(['roi_',type,'_dimZ'])= roiObj.spatialRef.ImageSize(3);

% Bounding box dimension x
diag.(['roi_',type,'_boxBoundDimX']) = boxBound(2,2) - boxBound(2,1) + 1;

% Bounding box dimension y
diag.(['roi_',type,'_boxBoundDimY']) = boxBound(1,2) - boxBound(1,1) + 1;

% Bounding box dimension z
diag.(['roi_',type,'_boxBoundDimZ']) = boxBound(3,2) - boxBound(3,1) + 1;

% Voxel dimension x
diag.(['roi_',type,'_voxDimX']) = roiObj.spatialRef.PixelExtentInWorldX;

% Voxel dimension y
diag.(['roi_',type,'_voxDimY']) = roiObj.spatialRef.PixelExtentInWorldY;

% Voxel dimension z
diag.(['roi_',type,'_voxDimZ']) = roiObj.spatialRef.PixelExtentInWorldZ;

% Voxel number
diag.(['roi_',type,'_voxNumb']) = numel(Xgl);

% Mean intensity
diag.(['roi_',type,'_meanInt']) = mean(Xgl);

% Minimum intensity
diag.(['roi_',type,'_minInt']) = min(Xgl);

% Maximum intensity
diag.(['roi_',type,'_maxInt']) = max(Xgl);
% -------------------------------------------------------------------------

end