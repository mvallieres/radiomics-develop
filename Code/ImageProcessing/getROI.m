function [volObj,roiObj] = getROI(sData,contourString,boxString,interp)
% -------------------------------------------------------------------------
% function [volObj,roiObj] = getROI(sData,contourNumber,box,interp)
% -------------------------------------------------------------------------
% DESCRIPTION:
% Computes the ROI box (+ smallest box containing the region of interest) 
% and associated mask from a 'sData' file.
% -------------------------------------------------------------------------
% INPUTS:
% - sData: Cell of structures organizing the data.
% - contourString: In the form '2' or '3-5+3. To be detailed shortly.
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
if ~strcmp(boxString,'full') && ~contains(boxString,'box')
    error('The third argument must either be "full" or contain the word "box"')
end
if nargin == 4
    if ~strcmp(interp,'noInterp') % The only option for the fourth argument is 'noInterp'
        interp = 'interp';
    end
else
    interp = 'interp';
end
[contourNumber,operations] = parseContourString(contourString); 


% INTIALIZATIONS
nContour = numel(contourNumber);
ROImask_cell = cell(1,nContour);
if strcmp(sData{2}.type,'PTscan') || strcmp(sData{2}.type,'CTscan') || strcmp(sData{2}.type,'MRscan') || strcmp(sData{2}.type,'ADCscan')
    scan = 'scan'; volume = 'volume';
elseif strcmp(sData{2}.type,'PTsim')
    scan = 'model'; volume = 'activity';
elseif strcmp(sData{2}.type,'MRsim')
    % DO SOMETHING!
end
spatialRef = sData{2}.(scan).(volume).spatialRef;
vol = single(sData{2}.(scan).(volume).data);


% COMPUTING ALL MASKS
for c = 1:nContour
    contour = contourNumber(c);
    
    % GETTING THE XYZ POINTS FROM THE sData STRUCTURE
    ROI_XYZ = sData{2}.scan.contour(contour).points_XYZ;

    % APPLYING ROTATION TO XYZ POINTS (if necessary --> MRscan)
    if isfield(sData{2}.scan.volume,'scanRot')
        ROI_XYZ = (sData{2}.scan.volume.scanRot * ROI_XYZ(:,1:3)')';
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
    ROImask_cell{c} = single(computeROI(ROI_XYZ,spatialRef,sData{2}.scan.orientation,sData{2}.type,interp));
end


% APPLYING OPERATIONS ON ALL MASKS
roi = ROImask_cell{1};
for c = 2:nContour
    if strcmp(operations(c-1),'+')
        roi = roi + ROImask_cell{c};
    elseif strcmp(operations(c-1),'-')
        roi = roi - ROImask_cell{c};
    end
    roi(roi >= 1) = 1;
    roi(roi < 1) = 0;
end


% COMPUTING THE BOUNDING BOX
[vol,roi,newSpatialRef] = computeBox(vol,roi,spatialRef,boxString);


% ARRANGE OUTPUT
volObj = struct; roiObj = struct;
volObj.data = vol ; volObj.spatialRef = newSpatialRef;
roiObj.data = roi; roiObj.spatialRef = newSpatialRef;

end