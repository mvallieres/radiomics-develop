function  localPeak = getLocPeak(imgObj,roiObj,res)
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

% - imgObj: Image object
% - roiObj: ROI objecy
% - res: [X,Y,Z] resolution vector in mm, e.g. [2,2,2]
% This works only in 3D for now.

% INITIALIZATION
distThresh = (3/(4*pi))^(1/3)*10; % About 6.2 mm, as defined in document

% Insert -inf outside ROI
temp = imgObj;
imgObj(roiObj == 0) = -Inf;

% Find the location(s) of the maximal voxel
maxVal = max(imgObj(:));
indMax = find(imgObj == maxVal);
[I,J,K] = ind2sub(size(imgObj),indMax);
nMax = numel(I);

% Reconverting to full object without -Inf
imgObj = temp;

% Get a meshgrid first
[X,Y,Z] = meshgrid(res(1).*((1:size(imgObj,2))-0.5),res(2).*((1:size(imgObj,1))-0.5),res(3).*((1:size(imgObj,3))-0.5)); % In mm

% Calculate the local peak
maxVal = -Inf;
for n = 1:nMax
    tempX = X - X(I(n),J(n),K(n));
    tempY = Y - Y(I(n),J(n),K(n));
    tempZ = Z - Z(I(n),J(n),K(n));
    tempDistMesh = sqrt(tempX.^2 + tempY.^2 + tempZ.^2);
    val = imgObj(tempDistMesh(:) <= distThresh);
    val(isnan(val(:))) = []; % This is not needed anymore
    if isempty(val)
        tempLocalPeak = imgObj(I(n),J(n),K(n));
    else
        tempLocalPeak = mean(val);
    end
    if tempLocalPeak > maxVal
        maxVal = tempLocalPeak;
    end
end

localPeak = maxVal;
end