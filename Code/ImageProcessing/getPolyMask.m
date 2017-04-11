function [ROImask] = getPolyMask(ROI_XYZ,spatialRef,orientation)
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


sz = spatialRef.ImageSize;
ROImask = zeros(sz(1),sz(2),sz(3));
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
for k = 1:numel(slices)
    ind = find(K == slices(k));
    ROImask(:,:,slices(k)) = or(ROImask(:,:,slices(k)),poly2mask(points(ind,a),points(ind,b),sz(1),sz(2)));
end

end