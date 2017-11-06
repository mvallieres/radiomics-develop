function [ROImask] = getPolygonMaskWithEdges(ROI_XYZ,spatialRef,orientation)
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


% April 5, 2017.
% IN CONTRAST TO getPolyMask THIS FUNCTION HAS NOT BEEN TESTED FOR SAGITTAL AND CORONAL SET OF IMAGES.
% NEEDS TO BE DONE TO BE TOTALLY SURE THIS IS WORKING.


sz = spatialRef.ImageSize;
ROImask = zeros(sz(1),sz(2),sz(3)); edges = zeros(sz(1),sz(2),sz(3));
[X,Y,Z] = worldToIntrinsic(spatialRef,ROI_XYZ(:,1),ROI_XYZ(:,2),ROI_XYZ(:,3)); % X,Y,Z in intrinsic image coordinates
points = [X,Y,Z];
if strcmp(orientation,'Axial')
    a = 1; b = 2; c = 3;
elseif strcmp(orientation,'Sagittal')
    a = 2; b = 3; c = 1;
elseif strcmp(orientation,'Coronal')
    a = 1; b = 3; c = 2;
end
K = round(points(:,c)); % Must assign the points to one slice
slices = unique(K);
[xq,yq] = meshgrid(1:sz(2),1:sz(1));
for k = 1:numel(slices)
    ind = find(K == slices(k));
    ROImask(:,:,slices(k)) = or(ROImask(:,:,slices(k)),inpolygon(xq,yq,points(ind,a),points(ind,b)));
    edgesIND = sub2ind(sz,int16(points(ind,b)),int16(points(ind,a)),int16(repmat(slices(k),[numel(ind),1])));
    edges(edgesIND) = 1;
end

ROImask = ROImask + edges;  % This allows to make sure that the actual XYZ points of the RTstruct are included in the ROI. It appears that this is not always the case with the function inpolygon.m.
ROImask(ROImask >= 1) = 1;
end
