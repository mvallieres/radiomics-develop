function [ROImask] = computeROI(ROI_XYZ,spatialRef,orientation,scanType,interp)
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


% HERE, ONLY THE DIMENSION OF SLICES IS ACTAULLY INTERPOLATED --> THIS IS
% THE ONLY RESOLUTION INFO WE CAN GET FROM THE RTstruct XYZ POINTS. 
% WE ASSUME THAT THE FUNCTION "poly2mask.m" WILL CORRECTLY CLOSE ANY 
% POLYGON IN THE IN-PLANE DIMENSION, EVEN IF WE GO FROM LOWER TO HIGHER 
% RESOLUTION (e.g. RTstruct created on PET and applied to CT)
% --> ALLOWS TO INTERPOLATE A RTstruct CREATED ON ANOTHER IMAGING VOLUME
%     WITH DIFFERENT RESOLUTIONS, BUT FROM THE SAME FRAM OF REFERENCE 
%     (e.g. T1w and T2w in MR scans, PET/CT, etc.)
% --> IN THE IDEAL AND RECOMMENDED CASE, A SPECIFIC RTstruct WAS CREATED AND
%     SAVED FOR EACH IMAGING VOLUME (SAFE PRACTICE)


% USING INTERPOLATION
while strcmp(interp,'interp')

    % Initialization
    if strcmp(orientation,'Axial')
        dimIJK = 3; dimXYZ = 3; direction = 'Z'; resXYZ = [spatialRef.PixelExtentInWorldX,spatialRef.PixelExtentInWorldY,0]; % Only the resolution in 'Z' will be changed
        dir1 = 'X'; dir2 = 'Y';
    elseif strcmp(orientation,'Sagittal')
        dimIJK = 1; dimXYZ = 2; direction = 'Y'; resXYZ = [spatialRef.PixelExtentInWorldX,0,spatialRef.PixelExtentInWorldZ]; % Only the resolution in 'Y' will be changed
        dir1 = 'X'; dir2 = 'Z';
    elseif strcmp(orientation,'Coronal')
        dimIJK = 2; dimXYZ = 1; direction = 'X'; resXYZ = [0,spatialRef.PixelExtentInWorldX,spatialRef.PixelExtentInWorldZ]; % Only the resolution in 'X' will be changed
        dir1 = 'Y'; dir2 = 'Z';
    end

    % Creating new imref3d object for sample points (with slice dimension similar to original volume where RTstruct was created)
    sliceSpacing = double(findSpacing(ROI_XYZ(:,dimIJK),scanType)); % Slice spacing in mm
    newSize = ceil(spatialRef.(['ImageExtentInWorld',direction])/sliceSpacing); %  Using "round" would yield the closest new size we can get. But using "ceil" is safer.
    resXYZ(dimXYZ) = sliceSpacing; 
    sz = spatialRef.ImageSize; sz(dimIJK) = newSize;
    newSpatialRef = imref3d(sz,resXYZ(1),resXYZ(2),resXYZ(3));
    newSpatialRef.([dir1,'WorldLimits']) = spatialRef.([dir1,'WorldLimits']); newSpatialRef.([dir2,'WorldLimits']) = spatialRef.([dir2,'WorldLimits']);
    diff = newSpatialRef.(['ImageExtentInWorld',direction]) - spatialRef.(['ImageExtentInWorld',direction]);
    if abs(diff) >= 0.01 % Sampled and queried volume are considered "different".
        newLimit = spatialRef.([direction,'WorldLimits'])(1) - (diff)/2;
        newSpatialRef.([direction,'WorldLimits']) = newSpatialRef.([direction,'WorldLimits']) - (newSpatialRef.([direction,'WorldLimits'])(1)-newLimit); % Sampled volume is now centered on queried volume.    
    else % Less than a 0.001 mm, sampled and queried volume are considered to be the same. At this point, spatialRef and newSpatialRef may have differed due to data manipulation, so we simply compute the ROI mask with spatialRef (i.e. simply using "poly2mask.m"), without performing interpolation.
        interp = 'noInterp';
        break % Getting out of the "while" statement
    end
    
    % Getting sampled points
    sz = newSpatialRef.ImageSize;
    [Xi,Yi,Zi] = meshgrid(1:sz(2),1:sz(1),1:sz(3));
    [X,Y,Z] = intrinsicToWorld(newSpatialRef,Xi,Yi,Zi);

    % Getting sampled mask --> newSpatialRef is actually for the "old original" volume (the sampled volume) on which the RTstruct was created.
    %V = getPolyMask(ROI_XYZ,newSpatialRef,orientation); % Using the poly2mask.m function. 
    V = getPolygonMask(ROI_XYZ,newSpatialRef,orientation); % Using the inpolygon.m function. To be further tested.

    % Getting query points (Xq,Yq,Zq) of output ROImask
    szQ = spatialRef.ImageSize;
    [Xqi,Yqi,Zqi] = meshgrid(1:szQ(2),1:szQ(1),1:szQ(3));
    [Xq,Yq,Zq] = intrinsicToWorld(spatialRef,Xqi,Yqi,Zqi);
    
    % Getting queried mask
    Vq = interp3(X,Y,Z,V,Xq,Yq,Zq,'cubic',0); 
    ROImask = Vq; ROImask(Vq < 0.5) = 0; ROImask(Vq >= 0.5) = 1;
    
    % Getting out of the "while" statement
    interp = 'NoMoreInterp';
end


% SIMPLY USING "poly2mask.m" or "inpolygon.m". "inpolygon.m" is slower, but
% apparently more accurate.
if strcmp(interp,'noInterp')
    %ROImask = getPolyMask(ROI_XYZ,spatialRef,orientation); % Using the poly2mask.m function.
    ROImask = getPolygonMask(ROI_XYZ,spatialRef,orientation); % Using the inpolygon.m function. To be further tested.
end

end
