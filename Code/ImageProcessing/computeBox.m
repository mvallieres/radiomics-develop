function [vol,roi,newSpatialRef] = computeBox(vol,roi,spatialRef,boxString)
% -------------------------------------------------------------------------
% AUTHOR(S): 
% - Martin Vallieres <mart.vallieres@gmail.com>
% -------------------------------------------------------------------------
% HISTORY:
% - Creation: March 2018
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

if contains(boxString,'box')
    comp = strcmp(boxString,'box');
    [boxBound] = computeBoundingBox(roi);
    if ~comp
        indBox = strfind(boxString,'box');
        if indBox(1) == 1 % Addition of a certain number of voxels in all dimensions
            nV = str2double(boxString((indBox(1)+3):end));
            nV = [nV,nV,nV];
        else % Multiplication of the size of the box
            factor = str2double(boxString(1:(indBox(1)-1)));
            sizeBox = [(boxBound(1,2) - boxBound(1,1) + 1),(boxBound(2,2) - boxBound(2,1) + 1),(boxBound(3,2) - boxBound(3,1) + 1)];
            newBox = sizeBox * factor;
            nV = round((newBox - sizeBox)/2);
        end

        ok = 0;
        while ~ok
            border = zeros(3,2);
            border(1,1) = boxBound(1,1) - nV(1); border(1,2) = boxBound(1,2) + nV(1);
            border(2,1) = boxBound(2,1) - nV(2); border(2,2) = boxBound(2,2) + nV(2);
            border(3,1) = boxBound(3,1) - nV(3); border(3,2) = boxBound(3,2) + nV(3);
            check1 = border(:,1) > 0; check1 = sum(check1(:));
            check2 = border(1,2) <= size(vol,1);
            check3 = border(2,2) <= size(vol,2);
            check4 = border(3,2) <= size(vol,3);
            check = check1 + check2 + check3 + check4;
            if check == 6
                ok = 1;
            else
                nV = floor(nV./2);
                if sum(nV(:)) == 0
                    ok = 1;
                    nV = [0,0,0];
                end
            end
        end
    else
        nV = [0,0,0]; % Will compute the smallest bounding box possible
    end
    boxBound(1,1) = boxBound(1,1) - nV(1); boxBound(1,2) = boxBound(1,2) + nV(1);
    boxBound(2,1) = boxBound(2,1) - nV(2); boxBound(2,2) = boxBound(2,2) + nV(2);
    boxBound(3,1) = boxBound(3,1) - nV(3); boxBound(3,2) = boxBound(3,2) + nV(3);
    vol = vol(boxBound(1,1):boxBound(1,2),boxBound(2,1):boxBound(2,2),boxBound(3,1):boxBound(3,2));
    roi = roi(boxBound(1,1):boxBound(1,2),boxBound(2,1):boxBound(2,2),boxBound(3,1):boxBound(3,2));

    res = [spatialRef.PixelExtentInWorldX,spatialRef.PixelExtentInWorldY,spatialRef.PixelExtentInWorldZ]; % Resolution in mm, nothing has changed here in terms of resolution. % XYZ format here.
    sizeBox = [(boxBound(1,2) - boxBound(1,1) + 1),(boxBound(2,2) - boxBound(2,1) + 1),(boxBound(3,2) - boxBound(3,1) + 1)]; % IJK, as required by imref3d
    [Xlimit,Ylimit,Zlimit] = intrinsicToWorld(spatialRef,boxBound(2,1),boxBound(1,1),boxBound(3,1));
    newSpatialRef = imref3d(sizeBox,res(1),res(2),res(3));
    newSpatialRef.XWorldLimits = newSpatialRef.XWorldLimits - (newSpatialRef.XWorldLimits(1)-(Xlimit-res(1)/2)); % The limit is defined as the border of the first pixel
    newSpatialRef.YWorldLimits = newSpatialRef.YWorldLimits - (newSpatialRef.YWorldLimits(1)-(Ylimit-res(2)/2));
    newSpatialRef.ZWorldLimits = newSpatialRef.ZWorldLimits - (newSpatialRef.ZWorldLimits(1)-(Zlimit-res(3)/2));
elseif strcmp(boxString,'full')
    newSpatialRef = spatialRef;
end

end
