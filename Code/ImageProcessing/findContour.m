function [contourString] = findContour(sData,nameROI,nameStructureSet)
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


% THIRD ARGUMENT IS OPTIONAL. HOWEVER, IF TWO ROIs WITH THE SAME NAME BUT
% FROM DIFFERENT STRUCTURE SETS ARE PRESENT, THE ALGORITHM WILL JUST OUTPUT
% THE FIRST CONTOUR IT FINDS, SO BEWARE. BETTER TO ALWAYS SUPPLY A
% STRUCTURE SET NAME.

delimiters = {'+','-'};
if isfield(sData{2},'nrrd') % Used by default if present
    nContourData = numel(sData{2}.nrrd.mask);
elseif isfield(sData{2},'img') % Used as second default if present
    nContourData = numel(sData{2}.img.mask);
else % Otherwise we use DICOM data
    nContourData = numel(sData{2}.scan.contour);
end

[nameROI,vectPlusMinus] = getSepROInames(nameROI,delimiters);
contourNumber = zeros(1,numel(nameROI));
if ~isempty(nameStructureSet)
    [nameStructureSet,~] = getSepROInames(nameStructureSet,delimiters);
    if numel(nameROI) ~= numel(nameStructureSet)
        error('The numbers of defined ROI names and Structure Set names are not the same')
    end
end

for i = 1:numel(nameROI)
    for j = 1:nContourData
        if isfield(sData{2},'nrrd') % Used by default if present
            nameTemp = sData{2}.nrrd.mask(j).name;
        elseif isfield(sData{2},'img') % Used as second default if present
            nameTemp = sData{2}.img.mask(j).name;
        else % Otherwise we use DICOM data
            nameTemp = sData{2}.scan.contour(j).name;
        end
        if strcmp(nameTemp,nameROI{i})
            if ~isempty(nameStructureSet)
                nameSetTemp = sData{2}.scan.contour(j).nameSet; % FOR DICOM + RTSTRUCT
                if strcmp(nameSetTemp,nameStructureSet{i})
                    contourNumber(i) = j;
                    break
                end
            else
                contourNumber(i) = j;
                break
            end
        end
    end
end
nROI = numel(contourNumber);
contourString = [num2str(contourNumber(1))];
for i = 2:nROI
    if vectPlusMinus(i-1) == 1
        sign = '+';
    elseif vectPlusMinus(i-1) == -1
        sign = '-';
    end
    contourString = [contourString,sign,num2str(contourNumber(i))];
end

end
