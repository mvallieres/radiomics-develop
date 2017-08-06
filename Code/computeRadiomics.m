function [radiomics] = computeRadiomics(volObjInit,roiObjInit,imParamScan)
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

% - imParamScan.image: Image processing parameters for the image space.
% - imParamScan.filter: Image processing parameters for filtered image space.
%   If this .filter structure is not present in imParamScan, no wavelet
%   features will be computed.
% - Eventually, "imParam" variable could be better organized, this is confusing.



% INITIALIZATION
radiomics.image = struct;
scaleNonText = imParamScan.image.interp.scaleNonText; 
volInterp = imParamScan.image.interp.volInterp; roiInterp = imParamScan.image.interp.roiInterp;
glRound = imParamScan.image.interp.glRound;
roiPV = imParamScan.image.interp.roiPV;
range = imParamScan.image.reSeg.range; if isempty(range), userSetMinVal = imParamScan.image.userSetMinVal; end
outliers = imParamScan.image.reSeg.outliers;
IH = imParamScan.image.discretisation.IH;
IVH = imParamScan.image.discretisation.IVH;
scaleText = imParamScan.image.interp.scaleText; 
algo = imParamScan.image.discretisation.texture.type; 
grayLevels = imParamScan.image.discretisation.texture.val;
nScale = numel(scaleText); nAlgo = numel(algo); nGl = numel(grayLevels{1}); nExp = nScale*nAlgo*nGl;
glcm = cell(nScale,nAlgo,nGl); glrlm = cell(nScale,nAlgo,nGl); glszm = cell(nScale,nAlgo,nGl); ngtdm = cell(nScale,nAlgo,nGl); gldzm = cell(nScale,nAlgo,nGl); ngldm = cell(nScale,nAlgo,nGl);
userSetMinVal = imParamScan.image.userSetMinVal;
type = imParamScan.image.type;
intensity = imParamScan.image.intensity; % Variable used to determine if there is 'arbitrary' (e.g., MRI) or 'definite' (e.g., CT) intensities.

% FILTERS INITIALIZATION
if isfield(imParamScan,'filter')
    filter = true;
    IHfilter = imParamScan.filter.discretisation.IH;
    IVHfilter = imParamScan.filter.discretisation.IVH;
    intensityFilter = imParamScan.filter.intensity;
    filtersType = imParamScan.filter.ToCompute; nFilters = numel(filtersType);
    for f = 1:nFilters
        if ~isempty(strfind(filtersType{f},'wavelet'))
            ind = strfind(filtersType{f},'_'); waveletName = filtersType{f}(ind+1:end); waveletName = replaceCharacter(waveletName,'.','dot');
            waveletType = {'LLL','LLH','LHL','LHH','HLL','HLH','HHL','HHH'}; nWav = numel(waveletType);
            for w = 1:nWav
                radiomics.([waveletType{w},'_',waveletName]) = struct;
            end
        else
            radiomics.(filtersType{f}) = struct;
        end
    end
else
    filter = false;
end


%%%%%%%%%%%%%%%%%%%%%%% COMPUTATION OF NON-TEXTURE FEATURES %%%%%%%%%%%%%%%%%%%%%%%
try
    tic, fprintf('--> Computation of non-texture features in image space: ')
    % STEP 1: INTERPOLATION
    [volObj] = interpVolume(volObjInit,scaleNonText,volInterp,glRound,'image');
    [roiObj_Morph] = interpVolume(roiObjInit,scaleNonText,roiInterp,roiPV,'roi');

    % STEP 2: RE-SEGMENTATION
    roiObj_Int = roiObj_Morph; % Now is the time to create the intensity mask
    [roiObj_Int.data] = rangeReSeg(volObj.data,roiObj_Int.data,range);
    [roiObj_Int.data] = outlierReSeg(volObj.data,roiObj_Int.data,outliers);
    
    % STEP 3: COMPUTE ALL NON-TEXTURE FEATURES IN IMAGE SPACE
    [imageStruct] = computeNonTextureFeatures(volObj,roiObj_Int,roiObj_Morph,scaleNonText,intensity,range,userSetMinVal,IH,IVH);
    [radiomics.image] = concatenateStruct(radiomics.image,imageStruct);
    toc
    
    % STEP 4: COMPUTE ALL NON-TEXTURE FEATURES IN ALL FILTERED SPACES
    if filter
        for f = 1:nFilters
            tic, fprintf('--> Computation of non-texture features in %s space: ',filtersType{f})
            [filterStruct] = computeNonTextureFeatures(volObj,roiObj_Int,roiObj_Morph,scaleNonText,intensityFilter,[],[],IHfilter,IVHfilter,filtersType{f});
            if ~isempty(strfind(filtersType{f},'wavelet'))
                for w = 1:nWav
                    radiomics.([waveletType{w},'_',waveletName]) = concatenateStruct(radiomics.([waveletType{w},'_',waveletName]),filterStruct.([waveletType{w},'_',waveletName]));
                end
            else
                [radiomics.(filtersType{f})] = concatenateStruct(radiomics.(filtersType{f}),filterStruct);
            end
            toc
        end
    end
    
catch
    fprintf('PROBLEM WITH PRE-PROCESSING OF NON-TEXTURE FEATURES')
    radiomics.image = 'ERROR_PROCESSING';
end
% -------------------------------------------------------------------------


% JE SUIS RENDU ICI!!


%%%%%%%%%%%%%%%%%%%%%%% COMPUTATION OF TEXTURE FEATURES %%%%%%%%%%%%%%%%%%%%%%%

% INITIALIZATION


count = 0;
for s = 1:nScale
    
    try
        tic, fprintf(['--> Pre-processing (interp + reSeg) for "Scale=',num2str(scaleText{s}(1)),'": '])
        
        % STEP 1: INTERPOLATION
        [volObj] = interpVolume(volObjInit,scaleText{s},volInterp,glRound,'image');
        [roiObj_Morph] = interpVolume(roiObjInit,scaleText{s},roiInterp,roiPV,'roi');

        % STEP 2: RE-SEGMENTATION
        roiObj_Int = roiObj_Morph; % Now is the time to create the intensity mask
        [roiObj_Int.data] = rangeReSeg(volObj.data,roiObj_Int.data,range);
        [roiObj_Int.data] = outlierReSeg(volObj.data,roiObj_Int.data,outliers);
        
        toc

        for a = 1:nAlgo
            for n = 1:nGl
                
                try
                    count = count + 1;
                    tic, fprintf(['--> Computation of texture features for "Scale=',num2str(scaleText{s}(1)),'", "Algo=',algo{a},'", "GL=',num2str(grayLevels{a}(n)),'" (',num2str(count),'/',num2str(nExp),'): '])
                
                    % Initialization
                    volInt_RE = roiExtract(volObj.data,roiObj_Int.data);
                    
                    % STEP 3: DISCRETISATION FOR TEXTURE FEATURES
                    [volQuant_RE] = discretisation(volInt_RE,algo{a},grayLevels{a}(n),userSetMinVal);

                    % STEP 4: COMPUTING ALL TEXTURE FEATURES
                    try
                        glcm{s,a,n} = getGLCMfeatures(volQuant_RE);
                    catch
                        fprintf('PROBLEM WITH COMPUTATION OF GLCM FEATURES ')
                        glcm{s,a,n} = 'ERROR_COMPUTATION';
                    end
                    try
                        glrlm{s,a,n} = getGLRLMfeatures(volQuant_RE);
                    catch
                        fprintf('PROBLEM WITH COMPUTATION OF GLRLM FEATURES ')
                        glrlm{s,a,n} = 'ERROR_COMPUTATION';                        
                    end
                    try
                        glszm{s,a,n} = getGLSZMfeatures(volQuant_RE);
                    catch
                        fprintf('PROBLEM WITH COMPUTATION OF GLSZM FEATURES ')
                        glszm{s,a,n} = 'ERROR_COMPUTATION';                         
                    end
                    try
                        gldzm{s,a,n} = getGLDZMfeatures(volQuant_RE,roiObj_Morph.data);
                    catch
                        fprintf('PROBLEM WITH COMPUTATION OF GLDZM FEATURES ')
                        gldzm{s,a,n} = 'ERROR_COMPUTATION';                            
                    end
                    try
                        ngtdm{s,a,n} = getNGTDMfeatures(volQuant_RE);
                    catch
                        fprintf('PROBLEM WITH COMPUTATION OF NGTDM FEATURES ')
                        ngtdm{s,a,n} = 'ERROR_COMPUTATION';                           
                    end
                    try
                        ngldm{s,a,n} = getNGLDMfeatures(volQuant_RE);
                    catch
                        fprintf('PROBLEM WITH COMPUTATION OF NGLDM FEATURES ')
                        ngldm{s,a,n} = 'ERROR_COMPUTATION';                          
                    end

                    toc
                catch
                    fprintf('PROBLEM WITH DISCRETISATION')
                    glcm{s,a,n} = 'ERROR_DISCRETISATION'; glrlm{s,a,n} = 'ERROR_DISCRETISATION'; glszm{s,a,n} = 'ERROR_DISCRETISATION'; gldzm{s,a,n} = 'ERROR_DISCRETISATION'; ngtdm{s,a,n} = 'ERROR_DISCRETISATION'; ngldm{s,a,n} = 'ERROR_DISCRETISATION';
                end
            end
        end
    catch 
        fprintf('PROBLEM WITH PRE-PROCESSING OF TEXTURE FEATURES')
        for a = 1:nAlgo
            for n = 1:nGl
                count = count + 1;
                glcm{s,a,n} = 'ERROR_PROCESSING'; glrlm{s,a,n} = 'ERROR_PROCESSING'; glszm{s,a,n} = 'ERROR_PROCESSING'; gldzm{s,a,n} = 'ERROR_PROCESSING'; ngtdm{s,a,n} = 'ERROR_PROCESSING'; ngldm{s,a,n} = 'ERROR_PROCESSING';
            end
        end
    end
end

radiomics.glcm = glcm;
radiomics.glrlm = glrlm;
radiomics.glszm = glszm;
radiomics.gldzm = gldzm;
radiomics.ngtdm = ngtdm;
radiomics.ngldm = ngldm;
% -------------------------------------------------------------------------

radiomics.imParam = imParamScan;
end