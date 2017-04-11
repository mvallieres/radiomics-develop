function [XYZ,faces,vertices] = getMesh(mask,res)
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

% - mask: Contains only 0's and 1's.
% - res: [a,b,c] vector specfying the resolution of the volume in mm.  % XYZ resolution (world), or JIK resolution (intrinsic matlab).
%
% - XYZ: [nPoints X 3] matrix of three column vectors, defining the [X,Y,Z]
%        positions of the points in the ROI (1's) of the mask volume.
%        --> In mm.
% - faces: [nPoints X 3] matrix of three column vectors, defining the [X,Y,Z]
%          positions of the faces of the isosurface of the mask (output
%          from "isosurface.m" function of MATLAB).
%          --> These are more precisely indexes to "vertices".
% - vertices: [nPoints X 3] matrix of three column vectors, defining the [X,Y,Z]
%             positions of the vertices of the isosurface of the mask (output
%             from "isosurface.m" function of MATLAB).
%             --> In mm.
%
% --> IMPORTANT: Make sure the "mask" is padded with a layer of 0's in all
%                dimensions to reduce potential isosurface computation errors.


% Getting the grid of X,Y,Z positions, where the coordinate reference
% system (0,0,0) is located at the upper left corner of the first voxel
% (-0.5: half a voxel distance). For the whole volume defining the mask, no
% matter if it is a 1 or a 0.
[X,Y,Z] = meshgrid(res(1).*((1:size(mask,2))-0.5),res(2).*((1:size(mask,1))-0.5),res(3).*((1:size(mask,3))-0.5));

% Getting the isosurface of the mask
[faces,vertices] = isosurface(X,Y,Z,mask,0.5);

% Getting the X,Y,Z positions of the ROI (i.e. 1's) of the mask
XYZ = [X(:),Y(:),Z(:)];
XYZ = XYZ(mask(:)==1,:);

end