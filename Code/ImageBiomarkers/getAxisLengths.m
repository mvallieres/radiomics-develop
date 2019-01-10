function [major,minor,least] = getAxisLengths(XYZ)
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

% - XYZ: [nPoints X 3] matrix of three column vectors, defining the [X,Y,Z]
%        positions of the points in the ROI (1's) of the mask volume.
%        --> In mm.


% Getting the geometric centre of mass
com_geom = sum(XYZ)/size(XYZ,1); % [1 X 3] vector

% Subtracting the centre of mass
XYZ(:,1) = XYZ(:,1) - com_geom(1);
XYZ(:,2) = XYZ(:,2) - com_geom(2);
XYZ(:,3) = XYZ(:,3) - com_geom(3);

% Getting the covariance matrix
covMat = cov(XYZ);

% Getting the eigenvalues
eigVal = eig(covMat);
major = eigVal(3);
minor = eigVal(2);
least = eigVal(1);

end