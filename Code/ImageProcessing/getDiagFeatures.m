function diag = getDiagFeatures(volObj,roiObj_Int,roiObj_Morph,type)
% -------------------------------------------------------------------------
% AUTHOR(S): Martin Vallieres <mart.vallieres@gmail.com>                   
% -------------------------------------------------------------------------
% HISTORY:                                                                 
% - Creation: February 2017                                                
% -------------------------------------------------------------------------
% STATEMENT:                                                                
% --> Copyright (C) 2017  Martin Vallieres                                 
%                                                                          
% This program is written on the basis of a scientific collaboration       
% between the following parties:                                           
% - Martin Vallieres <mart.vallieres@gmail.com>                            
% - Alex Zwanenburg <alexander.zwanenburg@nct-dresden.de>                  
%                                                                          
% By using this file, all parties acknowledge that it is to be kept private
% until publication of a potential scientific study.                       
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
[boxBound_Int] = computeBoundingBox(roiObj_Int.data);
[boxBound_Morph] = computeBoundingBox(roiObj_Morph.data); 
Xgl_Int = volObj.data(:); Xgl_Int = Xgl_Int(roiObj_Int.data(:) == 1);
Xgl_Morph = volObj.data(:); Xgl_Morph = Xgl_Morph(roiObj_Morph.data(:) == 1);

% Map dimension x
diag.(['roi_',type,'_Int_dimX']) = roiObj_Int.spatialRef.ImageSize(2);

% Map dimension y
diag.(['roi_',type,'_Int_dimY']) = roiObj_Int.spatialRef.ImageSize(1);

% Map dimension z
diag.(['roi_',type,'_Int_dimZ'])= roiObj_Int.spatialRef.ImageSize(3);

% Bounding box dimension x
diag.(['roi_',type,'_Int_boxBoundDimX']) = boxBound_Int(2,2) - boxBound_Int(2,1) + 1;

% Bounding box dimension y
diag.(['roi_',type,'_Int_boxBoundDimY']) = boxBound_Int(1,2) - boxBound_Int(1,1) + 1;

% Bounding box dimension z
diag.(['roi_',type,'_Int_boxBoundDimZ']) = boxBound_Int(3,2) - boxBound_Int(3,1) + 1;

% Bounding box dimension x
diag.(['roi_',type,'_Morph_boxBoundDimX']) = boxBound_Morph(2,2) - boxBound_Morph(2,1) + 1;

% Bounding box dimension y
diag.(['roi_',type,'_Morph_boxBoundDimY']) = boxBound_Morph(1,2) - boxBound_Morph(1,1) + 1;

% Bounding box dimension z
diag.(['roi_',type,'_Morph_boxBoundDimZ']) = boxBound_Morph(3,2) - boxBound_Morph(3,1) + 1;

% Voxel number
diag.(['roi_',type,'_Int_voxNumb']) = numel(Xgl_Int);

% Voxel number
diag.(['roi_',type,'_Morph_voxNumb']) = numel(Xgl_Morph);

% Mean intensity
diag.(['roi_',type,'_meanInt']) = mean(Xgl_Int);

% Minimum intensity
diag.(['roi_',type,'_minInt']) = min(Xgl_Int);

% Maximum intensity
diag.(['roi_',type,'_maxInt']) = max(Xgl_Int);
% -------------------------------------------------------------------------

end