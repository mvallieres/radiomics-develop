function [volObj,roiObj] = getROI(sData,contourNumber,box,interp)
% -------------------------------------------------------------------------
% function [volObj,roiObj] = getROI(sData,contourNumber,box,interp)
% -------------------------------------------------------------------------
% DESCRIPTION:
% Computes the ROI box (+ smallest box containing the region of interest) 
% and associated mask from a 'sData' file.
% -------------------------------------------------------------------------
% INPUTS:
% - sData: Cell of structures organizing the data.
% - contourNumber: Which contour to use in the computation of the tumor
%                  box. See sData{2}.scan.contour(contourNumber). If the
%                  argument is a single number, the mask associated to that
%                  contour will be created. If it is a vector, the masks
%                  associated to each number will be appended.
% - box:  String specifying the size if the box containing the region of interest.
%         --> 'full': Full imaging data as output.
%         --> 'box' computes the smallest bounding box.
%         --> Ex: 'box10': 10 voxels in all three dimensions are added to
%             the smallest bounding box. The number after 'box' defines the
%             number of voxels to add.
%         --> Ex: '2box': Computes the smallest box and outputs double its 
%             size. The number before 'box' defines the multiplication in
%             size.
% - interp: (optional). String specifying if we are to compute the ROI from
%           XYZ points solely using the function "inpolygon.m" of MATLAB.
%           This function can be considered safe when the RTstruct has been
%           saved specifically for the volume of interest. Otherwise, an
%            interpolation process (default; no argument) prior to 
%           "inpolygon.m" in the slice axis direction is recommended.
%         --> Ex: 'noInterp' (this is the only option, don't put a fourth argument).
%                 As a consequence: No Interpolation is performed in the
%                 slice axis dimensions (not recommended).
% -------------------------------------------------------------------------
% OUTPUTS:
% - volObj: 3D array of imaging data defining the smallest box
%           containing the region of interest.
%           --> In 'vol' format: vol.data is the 3D array, vol.spatialRef
%               is its associated imref3d object.
% - roiObj: 3D array of 1's and 0's defining the ROI in ROIbox.
%           --> In 'vol' format: vol.data is the 3D array, vol.spatialRef
%               is its associated imref3d object.
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




% PARSING OF ARGUMENTS
if ~strcmp(box,'full') && isempty(strfind(box,'box')) % FOR MATLAB 2016b and above, the last check can be replaced by "~contains(box,'box')"
    error('The third argument must either be "full" or contain the word "box"')
end
if nargin == 4
    if ~strcmp(interp,'noInterp') % The only option for the fourth argument is 'noInterp'
        interp = 'interp';
    end
else
    interp = 'interp';
end


% INTIALIZATIONS
nContour = numel(contourNumber);
ROImask_cell = cell(1,nContour);
if strcmp(sData{2}.type,'PTscan') || strcmp(sData{2}.type,'CTscan') || strcmp(sData{2}.type,'MRscan')
    scan = 'scan'; volume = 'volume';
elseif strcmp(sData{2}.type,'PTsim')
    scan = 'model'; volume = 'activity';
elseif strcmp(sData{2}.type,'MRsim')
    % DO SOMETHING!
end
spatialRef = sData{2}.(scan).(volume).spatialRef;
vol = sData{2}.(scan).(volume).data;


% COMPUTING ALL MASKS
for c = 1:nContour
    contour = contourNumber(c);
    
    % GETTING THE XYZ POINTS FROM THE sData STRUCTURE
    ROI_XYZ = sData{2}.scan.contour(contour).points_XYZ;

    % APPLYING ROTATION TO XYZ POINTS (if necessary --> MRscan)
    if isfield(sData{2}.scan.volume,'scanRot')
        ROI_XYZ = (sData{2}.scan.volume.scanRot * ROI_XYZ')';
    end

    % APPLYING TRANSLATION IF SIMULATION STRUCTURE AS INPUT (software STAMP utility)
    if isfield(sData{2}.(scan).(volume),'transScanToModel')
        translation = sData{2}.(scan).(volume).transScanToModel;
        ROI_XYZ(:,1) = ROI_XYZ(:,1) + translation(1);
        ROI_XYZ(:,2) = ROI_XYZ(:,2) + translation(2);
        ROI_XYZ(:,3) = ROI_XYZ(:,3) + translation(3);
    end

    % COMPUTING THE ROI MASK 
    % Problem here in computeROI.m: If the volume is a full-body CT and the slice
    % interpolation process occurs, a lot of RAM will be used in interp3.m.
    % One solution could be to a priori compute the bounding box before
    % computing the ROI (using XYZ points). But we still want the user to
    % be able to fully use the "box" argument, so we are fourré...TO SOLVE!
    ROImask_cell{c} = computeROI(ROI_XYZ,spatialRef,sData{2}.scan.orientation,sData{2}.type,interp);
end


% COMBINING ALL MASKS
roi = zeros(size(ROImask_cell{1}));
for c = 1:nContour
    roi = roi + ROImask_cell{c};
end
roi(roi >= 1) = 1;
roi(roi < 1) = 0;


% COMPUTING THE BOUNDING BOX
if ~isempty(strfind(box,'box')) %  FOR MATLAB 2016b and above, this can be modifed with the function "contains"
    comp = strcmp(box,'box');
    [boxBound] = computeBoundingBox(roi);
    if ~comp
        indBox = strfind(box,'box');
        if indBox(1) == 1 % Addition of a certain number of voxels in all dimensions
            nV = str2double(box((indBox(1)+3):end));
            nV = [nV,nV,nV];
        else % Multiplication of the size of the box
            factor = str2double(box(1:(indBox(1)-1)));
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
elseif strcmp(box,'full')
    newSpatialRef = spatialRef;
end


% ARRANGE OUTPUT
volObj = struct; roiObj = struct;
volObj.data = vol ; volObj.spatialRef = newSpatialRef;
roiObj.data = roi; roiObj.spatialRef = newSpatialRef;

end