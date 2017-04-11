function area = getMeshArea(faces,vertices)
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

% - faces: [nPoints X 3] matrix of three column vectors, defining the [X,Y,Z]
%          positions of the faces of the isosurface or convex hull of the mask 
%          (output from "isosurface.m" or "convhull.m" functions of MATLAB).
%          --> These are more precisely indexes to "vertices"
% - vertices: [nPoints X 3] matrix of three column vectors, defining the [X,Y,Z]
%             positions of the vertices of the isosurface of the mask (output
%             from "isosurface.m" function of MATLAB).
%             --> In mm.


% Getting two vectors of edges for each face
a = vertices(faces(:,2),:) - vertices(faces(:,1),:);
b = vertices(faces(:,3),:) - vertices(faces(:,1),:);

% Calculating the surface area of each face and summing it up all at once.
c = cross(a,b,2);
area = 1/2 * sum(sqrt(sum(c.^2,2)));

end