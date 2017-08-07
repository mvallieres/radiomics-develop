function [resultsStruct] = computeTextureFeatures(volObjImage,roiObj_Int,roiObj_Morph,scaleText,algo,grayLevels,userSetMinVal,filterType)
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



% *************************************************************************
% INITIALIZATIONS
if nargin == 8
    filter = true;
else
    filter = false;
end

% Scale name
scaleName = num2str(scaleText(1)); scaleName = replaceCharacter(scaleName,'.','dot'); % Always isotropic resampling, so the first entry is ok.
scaleName = ['scale',scaleName];

% Discretisation name
grayLevelsName = num2str(grayLevels); grayLevelsName = replaceCharacter(grayLevelsName,'.','dot');
if ~isempty(strfind(algo,'FBS')) % The minimum value defines the computation.
    if ~isempty(userSetMinVal)
        minValName = num2str(userSetMinVal); minValName = replaceCharacter(minValName,'.','dot'); minValName = replaceCharacter(minValName,'-','M');
        minValName = ['_min',minValName];
    else % Otherwise, minimum value of ROI will be used (not recommended), so no need to report it.
        minValName = [];
    end
else
    minValName = [];
end
discretisationName = ['algo',algo,'_bin',grayLevelsName,minValName];

% Processing full name
processingName = [scaleName,'_',discretisationName];

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
    try
        volInt_RE = roiExtract(volObj.data,roiObj_Int.data);

        % DISCRETISATION FOR TEXTURE FEATURES
        [volQuant_RE] = discretisation(volInt_RE,algo,grayLevels,userSetMinVal);

        % COMPUTING ALL TEXTURE FEATURES
        try
            results.glcm.(processingName) = getGLCMfeatures(volQuant_RE);
        catch
            fprintf('PROBLEM WITH COMPUTATION OF GLCM FEATURES ')
            results.glcm.(processingName).error = 'ERROR_COMPUTATION';
        end
        try
            results.glrlm.(processingName) = getGLRLMfeatures(volQuant_RE);
        catch
            fprintf('PROBLEM WITH COMPUTATION OF GLRLM FEATURES ')
            results.glrlm.(processingName).error = 'ERROR_COMPUTATION';                        
        end
        try
            results.glszm.(processingName) = getGLSZMfeatures(volQuant_RE);
        catch
            fprintf('PROBLEM WITH COMPUTATION OF GLSZM FEATURES ')
            results.glszm.(processingName).error = 'ERROR_COMPUTATION';                         
        end
        try
            results.gldzm.(processingName) = getGLDZMfeatures(volQuant_RE,roiObj_Morph.data);
        catch
            fprintf('PROBLEM WITH COMPUTATION OF GLDZM FEATURES ')
            results.gldzm.(processingName).error = 'ERROR_COMPUTATION';                            
        end
        try
            results.ngtdm.(processingName) = getNGTDMfeatures(volQuant_RE);
        catch
            fprintf('PROBLEM WITH COMPUTATION OF NGTDM FEATURES ')
            results.ngtdm.(processingName).error = 'ERROR_COMPUTATION';                           
        end
        try
            results.ngldm.(processingName) = getNGLDMfeatures(volQuant_RE);
        catch
            fprintf('PROBLEM WITH COMPUTATION OF NGLDM FEATURES ')
            results.ngldm.(processingName).error = 'ERROR_COMPUTATION';                          
        end
    catch
        fprintf('PROBLEM WITH DISCRETISATION')
        results.glcm.(processingName).error = 'ERROR_DISCRETISATION';
        results.glrlm.(processingName).error = 'ERROR_DISCRETISATION';
        results.glszm.(processingName).error = 'ERROR_DISCRETISATION';
        results.gldzm.(processingName).error = 'ERROR_DISCRETISATION';
        results.ngtdm.(processingName).error = 'ERROR_DISCRETISATION';
        results.ngldm.(processingName).error = 'ERROR_DISCRETISATION';
    end
        
    if nObjects == 1
        resultsStruct = results;
    else
        resultsStruct.(nameTypes{o}) = results;
    end
end

end