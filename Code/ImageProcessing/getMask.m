function [volObj,roiObj] = getMask(sData,contourString,formatData,boxString)
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

% PARSING OF ARGUMENTS
if ~strcmp(boxString,'full') && ~contains(boxString,'box')
    error('The fourth argument must either be "full" or contain the word "box"')
end
if ~strcmp(formatData,'nii') && ~strcmp(formatData,'nrrd') && ~strcmp(formatData,'img')
    error('The third argument must either be ''nii'' or ''nrrd'' or ''img''')
end
[contourNumber,operations] = parseContourString(contourString);


% INTIALIZATIONS
nContour = numel(contourNumber);
ROImask_cell = cell(1,nContour);
spatialRef = sData{2}.scan.volume.spatialRef; % GETTING THE ONE FROM THE DICOM DATA. DICOM DATA THUS MUST BE PRESENT!
vol = single(sData{2}.(formatData).volume.data);


% COMPUTING ALL MASKS
for c = 1:nContour
    contour = contourNumber(c);
    ROImask_cell{c} = single(sData{2}.(formatData).mask(contour).data);
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
