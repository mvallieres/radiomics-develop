function [volObjQ] = interpVolume(volObjS,voxDim,interpMet,roundVal,type,boxElements)
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

% --> voxDim: The following format is used [Xin,Yin,Zslice], where Xin and Yin are the X (left to right) and Y (bottom to top) IN-PLANE resolutions, and Zslice is the slice spacing, NO MATTER THE ORIENTATION OF THE VOLUME (i.e. axial , sagittal, coronal).   

% -------------------------------------------------------------------------
% PARSING ARGUMENTS
if isempty(voxDim) || sum(voxDim) == 0
    volObjQ = volObjS;
    return
else
    if numel(voxDim) == 2
        twoD = true;
    else
        twoD = false;
    end
end
if nargin < 5
    error('The number of arguments must at least be 5')
end
if nargin < 6 % Sixth argument is optional
    useBox = false;
else
    useBox = true;
end
if ~strcmp(type,'image') && ~strcmp(type,'roi')
    error('Fifth argument must either be "image" or "roi"')
else
    if strcmp(type,'image')
        if ~strcmp(interpMet,'linear') && ~strcmp(interpMet,'cubic') && ~strcmp(interpMet,'spline')
            error('Third argument must either be "linear" or "cubic" or "spline"')
        else
            val = ~mod(log10(roundVal),1);
            if ~val
                error('~mod(log10(roundVaqueried points ("q" or "Q")l),1) for the fourth argument must be true')
            end
        end
    elseif strcmp(type,'roi')
        if ~strcmp(interpMet,'nearest') && ~strcmp(interpMet,'linear') && ~strcmp(interpMet,'cubic')
            error('Third argument must either be "nearest" or "linear" or "cubic"')
        else
            if roundVal < 0 || roundVal > 1
                error('Fourth argument must be between 0 and 1')
            end
        end
    end
end
% -------------------------------------------------------------------------



% --> QUERIED POINTS: NEW INTERPOLATED VOLUME: "q" or "Q".
% --> SAMPLED POINTS: ORIGINAL VOLUME: "s" or "S".
% --> Always using XYZ coordinates (unless specifically noted), not MATLAB IJK, so beware!

% INITIALIZATION
resQ = voxDim;
if twoD, resQ = [resQ,volObjS.spatialRef.PixelExtentInWorldZ]; end % If 2D, the resolution of the slice dimension of he queried volume is set to the same as the sampled volume.
resS = [volObjS.spatialRef.PixelExtentInWorldX,volObjS.spatialRef.PixelExtentInWorldY,volObjS.spatialRef.PixelExtentInWorldZ];
if isequal(resS,resQ)
    volObjQ = volObjS;
    return
end
spatialRefS = volObjS.spatialRef;
sizeS = spatialRefS.ImageSize;
extentS = [spatialRefS.ImageExtentInWorldX,spatialRefS.ImageExtentInWorldY,spatialRefS.ImageExtentInWorldZ];
lowLimitsS = [spatialRefS.XWorldLimits(1),spatialRefS.YWorldLimits(1),spatialRefS.ZWorldLimits(1)];


% CREATING QUERIED "imref3d" OBJECT CENTERED ON SAMPLED VOLUME
sizeQ = ceil(round(extentS./resQ*1000)/1000); % Rounding up extentS./resQ to 3 decimal places (minimize numerical errors) before applying the IBSI-compliant "ceil" operation (safer than round).
temp = sizeQ(1); sizeQ(1) = sizeQ(2); sizeQ(2) = temp; % Switching to IJK (matlab) reference frame for "imref3d" computation.
if twoD, sizeQ(3) = volObjS.spatialRef.ImageSize(3); end % If 2D, forcing the size of the queried volume in the slice dimension to be the same as the sample volume.
spatialRefQ = imref3d(sizeQ,resQ(1),resQ(2),resQ(3));
extentQ = [spatialRefQ.ImageExtentInWorldX,spatialRefQ.ImageExtentInWorldY,spatialRefQ.ImageExtentInWorldZ];
lowLimitsQ = [spatialRefQ.XWorldLimits(1),spatialRefQ.YWorldLimits(1),spatialRefQ.ZWorldLimits(1)];
diff = extentQ - extentS; newLowLimitsQ = lowLimitsS - diff/2;
spatialRefQ.XWorldLimits = spatialRefQ.XWorldLimits - (lowLimitsQ(1) - newLowLimitsQ(1));
spatialRefQ.YWorldLimits = spatialRefQ.YWorldLimits - (lowLimitsQ(2) - newLowLimitsQ(2));
spatialRefQ.ZWorldLimits = spatialRefQ.ZWorldLimits - (lowLimitsQ(3) - newLowLimitsQ(3));


% REDUCE THE SIZE OF THE VOLUME PRIOR TO INTERPOLATION
if useBox
    boxString = boxElements.boxString;
    roiObjS = boxElements.roiObj;
    [~,~,tempSpatialRef] = computeBox(roiObjS.data,roiObjS.data,roiObjS.spatialRef,boxString); % Computing a new spatialRef for a smaller box as defined by 'boxString'.
    sizeTemp = tempSpatialRef.ImageSize;
    [Xbound,Ybound,Zbound] = intrinsicToWorld(tempSpatialRef,[1;sizeTemp(2)],[1;sizeTemp(1)],[1;sizeTemp(3)]); % Getting world boundaries (center of voxels) of the new box.
    [Xbound,Ybound,Zbound] = worldToIntrinsic(spatialRefQ,Xbound,Ybound,Zbound); % Getting the image positions of the boundaries of the new box, IN THE FULL QUERIED FRAME OF REFERENCE (centered on the sampled frame of reference).
    Xbound = round(Xbound); Ybound = round(Ybound); Zbound = round(Zbound); % Rounding to the nearest image position integer
    sizeQ = [Ybound(2) - Ybound(1) + 1,Xbound(2) - Xbound(1) + 1,Zbound(2) - Zbound(1) + 1];
    [Xbound,Ybound,Zbound] = intrinsicToWorld(spatialRefQ,Xbound,Ybound,Zbound); % Converting back to world positions ion order to correctly define edges of the new box and thus center it onto the full queried reference frame
    newLowLimitsQ(1) = Xbound(1) - resQ(1)/2; newLowLimitsQ(2) = Ybound(1) - resQ(2)/2; newLowLimitsQ(3) = Zbound(1) - resQ(3)/2;
    spatialRefQ = imref3d(sizeQ,resQ(1),resQ(2),resQ(3));
    spatialRefQ.XWorldLimits = spatialRefQ.XWorldLimits - (spatialRefQ.XWorldLimits(1) - newLowLimitsQ(1));
    spatialRefQ.YWorldLimits = spatialRefQ.YWorldLimits - (spatialRefQ.YWorldLimits(1) - newLowLimitsQ(2));
    spatialRefQ.ZWorldLimits = spatialRefQ.ZWorldLimits - (spatialRefQ.ZWorldLimits(1) - newLowLimitsQ(3));
end


% CREATING QUERIED XYZ POINTS
[Xq,Yq,Zq] = meshgrid(1:sizeQ(2),1:sizeQ(1),1:sizeQ(3)); % Since sizeS is in the IJK (matlab) reference frame
[Xq,Yq,Zq] = intrinsicToWorld(spatialRefQ,Xq,Yq,Zq);


% CREATING SAMPLED XYZ POINTS
[Xs,Ys,Zs] = meshgrid(1:sizeS(2),1:sizeS(1),1:sizeS(3)); % Since sizeQ is in the IJK (matlab) reference frame
[Xs,Ys,Zs] = intrinsicToWorld(spatialRefS,Xs,Ys,Zs);


% INTERPOLATING VOLUME
volObjQ = struct;
if strcmp(interpMet,'spline')
    [volObjQ.data] = interp3(Xs,Ys,Zs,volObjS.data,Xq,Yq,Zq,interpMet); % Values outside the domain are interpolated.
else
    [volObjQ.data] = interp3(Xs,Ys,Zs,volObjS.data,Xq,Yq,Zq,interpMet,0); % Values outside the domain are set to 0.
end
volObjQ.spatialRef = spatialRefQ;


% ROUNDING
switch type
    case 'image' % Grey level rounding for 'image' type
        if ~isempty(roundVal)
            volObjQ.data = round(volObjQ.data/roundVal)*roundVal;
        end
    case 'roi' % ROI partial volume for 'roi' type
        volObjQ.data(volObjQ.data >= roundVal) = 1;
        volObjQ.data(volObjQ.data < roundVal) = 0;
end

end