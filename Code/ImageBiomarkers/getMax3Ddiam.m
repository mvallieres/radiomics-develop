function  sizeROI = getMax3Ddiam(faces,vertices)
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


% Finding the max distance between all pair or points of the convex hull
max = 0;
nPoints = size(faces,1);
for i = 1:nPoints
    for j = (i + 1):nPoints
        dist = (vertices(faces(i,1),1) - vertices(faces(j,1),1))^2 + (vertices(faces(i,2),2) - vertices(faces(j,2),2))^2 + (vertices(faces(i,3),3) - vertices(faces(j,3),3))^2;
        if dist > max
            max = dist;
        end
    end
end

sizeROI = sqrt(max);

end