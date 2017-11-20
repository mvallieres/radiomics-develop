function [resultsStruct] = computeNonTextureFeatures(volObjImage,roiObj_Int,roiObj_Morph,scaleNonText,intensity,range,userSetMinVal,IH,IVH,filterType)
% -------------------------------------------------------------------------
% AUTHOR(S): 
% - Martin Vallieres <mart.vallieres@gmail.com>
% -------------------------------------------------------------------------
% HISTORY:
% - Creation: August 2017
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

% NOTE: FOR IH features, we only need userSetMinVal in case a FBS algorithm
% is use for discretisation. This value is set to range(1) if it exists before entering that function.
% However, we also need range(2) for IVH features when using FBS, so we pass on both
% arguments in this function.
% - IH must never be empty.
% - IVH may be empty in the case of CT.


% *************************************************************************
% INITIALIZATIONS (this is a bit complicated, it should be simplifed)
if nargin == 10
    filter = true;
else
    filter = false;
end
if sum(scaleNonText) == 0 % In case the user chose to not interpolate
    scaleNonText = [volObjImage.spatialRef.PixelExtentInWorldX,volObjImage.spatialRef.PixelExtentInWorldY,volObjImage.spatialRef.PixelExtentInWorldZ];
else
    if numel(scaleNonText) == 2 % In case not interpolation is performed in the slice direction (e.g. 2D case)
        scaleNonText = [scaleNonText,volObjImage.spatialRef.PixelExtentInWorldZ];
    end
end

% Scale name
scaleName = num2str(scaleNonText(1)); scaleName = replaceCharacter(scaleName,'.','dot'); % Always isotropic resampling, so the first entry is ok.
scaleName = ['scale',scaleName];

% IH name
IHvalName = num2str(IH.val); IHvalName = replaceCharacter(IHvalName,'.','dot');
IHvalName = ['bin',IHvalName];
if ~isempty(strfind(IH.type,'FBS')) % The minimum value defines the computation.
    if ~isempty(userSetMinVal)
        minValName = num2str(userSetMinVal); minValName = replaceCharacter(minValName,'.','dot'); minValName = replaceCharacter(minValName,'-','M');
        minValName = ['_min',minValName];
    else % Otherwise, minimum value of ROI will be used (not recommended), so no need to report it.
        minValName = [];
    end
else
    minValName = [];
end
IHname = [scaleName,'_algo',IH.type,'_',IHvalName,minValName];

% IVH name
if isempty(IVH)
    IVHAlgoName = 'algoNone'; IVHvalName = ['bin1'];
    if ~isempty(range) % The range defines the computation.
        minValName = num2str(range(1)); minValName = replaceCharacter(minValName,'.','dot'); minValName = replaceCharacter(minValName,'-','M');
        maxValName = num2str(range(2)); maxValName = replaceCharacter(maxValName,'.','dot'); maxValName = replaceCharacter(maxValName,'-','M');
        rangeName = ['_min',minValName,'_max',maxValName];
    else 
        rangeName = [];
    end
else
    IVHAlgoName = ['algo',IVH.type];
    IVHvalName = num2str(IVH.val); IVHvalName = replaceCharacter(IVHvalName,'.','dot');
    IVHvalName = ['bin',IVHvalName];
    if ~isempty(strfind(IH.type,'FBS')) % The range defines the computation.
        if ~isempty(range)
            minValName = num2str(range(1)); minValName = replaceCharacter(minValName,'.','dot'); minValName = replaceCharacter(minValName,'-','M');
            maxValName = num2str(range(2)); maxValName = replaceCharacter(maxValName,'.','dot'); maxValName = replaceCharacter(maxValName,'-','M');
            if strcmp(maxValName,'Inf') % In this case, the maximum value of the ROI is used, so no need to report it.
                rangeName = ['_min',minValName];
            else
                rangeName = ['_min',minValName,'_max',maxValName];
            end
        else % min-max of ROI will be used, no need to report it.
            rangeName = [];
        end
    else % min-max of ROI will be used, no need to report it.
        rangeName  = [];
    end
end
IVHname = [IVHAlgoName,'_',IHvalName,rangeName];
% -------------------------------------------------------------------------



% PREPARATION OF COMPUTATION
if filter && ~isempty(strfind(filterType,'wavelet'))
    ind = strfind(filterType,'_'); waveletName = filterType(ind+1:end);
    [subbands] = getWaveletSubbands(volObjImage.data,waveletName);
    nameTypes = fieldnames(subbands); nType = numel(nameTypes); volObjCell = cell(1,nType);
    for s = 1:nType
        volObjCell{s} = volObjImage; volObjCell{s}.data = subbands.(nameTypes{s});
        subbands.(nameTypes{s}) = []; % Just to clear some memory
    end
else
    volObjCell = cell(1);
    volObjCell{1} = {volObjImage};
    volObjCell = volObjCell{1};
end



% COMPUTATION
resultsStruct = struct;
nObjects = numel(volObjCell);
for o = 1:nObjects
    results = struct;
    volObj = volObjCell{o};
    volInt_RE = roiExtract(volObj.data,roiObj_Int.data);
    try
        % STEP 3: CALCULATION OF MORPHOLOGICAL FEATURES
        results.morph.(scaleName) = getMorphFeatures(volObj.data,roiObj_Int.data,roiObj_Morph.data,scaleNonText,intensity); % For scans with arbitrary units, some features will not be computed.
    catch
        fprintf('PROBLEM WITH COMPUTATION OF MORPHOLOGICAL FEATURES ')
        results.morph.(scaleName).error = 'ERROR_COMPUTATION';
    end

    try
        % STEP 4: CALCULATION OF LOCAL INTENSITY FEATURES
        results.locInt.(scaleName) = getLocIntFeatures(volObj.data,roiObj_Int.data,scaleNonText,intensity); % For scans with arbitrary units, all of these features will not be computed.
    catch
        fprintf('PROBLEM WITH COMPUTATION OF LOCAL INTENSITY FEATURES ')
        results.locInt.(scaleName).error = 'ERROR_COMPUTATION';        
    end

    try
        % STEP 5: CALCULATION OF STATISTICAL FEATURES
        results.stats.(scaleName) = getStatsFeatures(volInt_RE,intensity); % For scans with arbitrary units, some features will not be computed.
    catch
        fprintf('PROBLEM WITH COMPUTATION OF STATISTICAL FEATURES ')
        results.stats.(scaleName).error = 'ERROR_COMPUTATION';           
    end

    try
        % STEP 6: CALCULATION OF INTENSITY HISTOGRAM FEATURES
        [volQuant_RE] = discretisation(volInt_RE,IH.type,IH.val,userSetMinVal); % There would actually be no need to include "userSetMinVal" here as fourth argument, as discretisation for IH features is (logically) always set to "FBN". This value will not be used for "FBN", see discretisation.m code. But this is a safety check in case FBS is used.
        results.intHist.(IHname) = getIntHistFeatures(volQuant_RE);
    catch
        fprintf('PROBLEM WITH COMPUTATION OF INTENSITY HISTOGRAM FEATURES ')
        results.intHist.(IHname).error = 'ERROR_COMPUTATION';           
    end

    try
        % STEP 7: CALCULATION OF INTENSITY-VOLUME HISTOGRAM FEATURES
        if ~isempty(IVH), [volQuant_RE,wd] = discretisation(volInt_RE,IVH.type,IVH.val,userSetMinVal,'ivh'); else volQuant_RE = volInt_RE; wd = 1; end % FOR CT, WE DO NOT WANT TO DISCRETISE. AN EMPTY IVH STRUCT ([]) DEFINES WHAT WE WANT TO USE FOR CT. FOR PET: FBS/0.1; FOR MRI: FBN/1000.
        if ~isempty(IVH)
            if strcmp(IVH.type,'FBS') || strcmp(IVH.type,'FBSequal') % PET example case (definite intensity units -- continuous case)
                rangeFBS = zeros(1,2);
                if isempty(range)
                    rangeFBS(1) = min(volInt_RE(:));
                    rangeFBS(2) = max(volInt_RE(:));
                else
                    if range(1) == -Inf
                        rangeFBS(1) = min(volInt_RE(:));
                    else
                        rangeFBS(1) = range(1);
                    end
                    if range(2) == Inf
                        rangeFBS(2) = max(volInt_RE(:));
                    else
                        rangeFBS(2) = range(2);
                    end
                end
                rangeFBS(1) = rangeFBS(1) + 0.5*wd; % In this case, wd = wb (see discretisation.m)
                rangeFBS(2) = rangeFBS(2) - 0.5*wd; % In this case, wd = wb (see discretisation.m)
                results.intVolHist.(IVHname) = getIntVolHistFeatures(volQuant_RE,wd,rangeFBS);
            else % MRI example case (arbitrary intensity units)
                results.intVolHist.(IVHname) = getIntVolHistFeatures(volQuant_RE,wd); 
            end
        else % CT example case (definite intensity units -- discrete case)
            results.intVolHist.(IVHname) = getIntVolHistFeatures(volQuant_RE,wd,range); 
        end
    catch
        fprintf('PROBLEM WITH COMPUTATION OF INTENSITY-VOLUME HISTOGRAM FEATURES ')
        results.intVolHist.(IVHname).error = 'ERROR_COMPUTATION';         
    end
    
    if nObjects == 1
        resultsStruct = results;
    else
        resultsStruct.(nameTypes{o}) = results;
    end
end
    
end