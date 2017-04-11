function contourNumber = findContour(sData,nameROI,nameStructureSet)
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

delimiter = ',';
nContourData = numel(sData{2}.scan.contour);

nameROI = getSepROInames(nameROI,delimiter);
contourNumber = zeros(1,numel(nameROI));
if nargin == 3
    nameStructureSet = getSepROInames(nameStructureSet,delimiter);
    if numel(nameROI) ~= numel(nameStructureSet)
        error('The numbers of defined ROI names and Structure Set names are not the same')
    end
end

for i = 1:numel(nameROI);
    for j = 1:nContourData
        nameTemp = sData{2}.scan.contour(j).name;
        if strcmp(nameTemp,nameROI{i})
            if nargin == 3
                nameSetTemp = sData{2}.scan.contour(j).nameSet;
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
ind = ~contourNumber;
contourNumber = contourNumber(~ind); % Removing contours that were not found. WARNING: We may never notice that one contour is not used due to a mistake in nameROI as input.

end
