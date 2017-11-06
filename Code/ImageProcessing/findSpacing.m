function sliceSpacing = findSpacing(points,scanType)
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

% README --> This function works for points from at least 2 slices. If only
% one slice is present, the function returns a NaN.


if strcmp(scanType,'MRscan')
    slices = unique(round(points*10)/10); % Rounding to the nearest 0.1 mm, MRI is more problematic due to arbitrary orientations allowed for imaging volumes.
else
    slices = unique(round(points*100)/100); % Rounding to the nearest 0.01 mm 
end
nSlices = numel(slices);
diff = zeros(nSlices - 1,1);
for s = 1:nSlices - 1
    diff(s) = slices(s + 1) - slices(s);
end

[sliceSpacing,nOcc] = mode(diff);
if nOcc == 1
    sliceSpacing = mean(diff);
end

end