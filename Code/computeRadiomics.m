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
radiomics = struct;
scaleNonText = imParamScan.image.interp.scaleNonText; 
volInterp = imParamScan.image.interp.volInterp; roiInterp = imParamScan.image.interp.roiInterp;
glRound = imParamScan.image.interp.glRound;
roiPV = imParamScan.image.interp.roiPV;
range = imParamScan.image.reSeg.range;
outliers = imParamScan.image.reSeg.outliers;
IH = imParamScan.image.discretisation.IH;
IVH = imParamScan.image.discretisation.IVH;
scaleText = imParamScan.image.interp.scaleText; 
algo = imParamScan.image.discretisation.texture.type; 
grayLevels = imParamScan.image.discretisation.texture.val;
nScale = numel(scaleText); nAlgo = numel(algo); nGl = numel(grayLevels{1}); nExp = nScale*nAlgo*nGl;
glcm = cell(nScale,nAlgo,nGl); glrlm = cell(nScale,nAlgo,nGl); glszm = cell(nScale,nAlgo,nGl); ngtdm = cell(nScale,nAlgo,nGl); gldzm = cell(nScale,nAlgo,nGl); ngldm = cell(nScale,nAlgo,nGl);
type = imParamScan.image.type;
intensity = imParamScan.image.intensity; % Variable used to determine if there is 'arbitrary' (e.g., MRI) or 'definite' intensities (

% SETTING UP userSetMinVal
if ~isempty(range)
    userSetMinVal = range(1);
else
    if strcmp(type,'PTscan')
        userSetMinVal = 0; % SUV is defined from 0 to inf.
    elseif strcmp(type,'CTscan')
        userSetMinVal = -1000; % No assumption on range of HU.
    elseif strcmp(type,'MRscan')
        userSetMinVal = 0; % This a dummy value, has no consequence on the rest of the code since FBS discretisation is not used for MRI.
    end
end



%%%%%%%%%%%%%%%%%%%%%%% COMPUTATION OF NON-TEXTURE FEATURES %%%%%%%%%%%%%%%%%%%%%%%
try
    tic, fprintf('--> Computation of non-texture features (1/1): ')
    
    % STEP 1: INTERPOLATION
    [volObj] = interpVolume(volObjInit,scaleNonText,volInterp,glRound,'image');
    [roiObj_Morph] = interpVolume(roiObjInit,scaleNonText,roiInterp,roiPV,'roi');

    % STEP 2: RE-SEGMENTATION
    roiObj_Int = roiObj_Morph; % Now is the time to create the intensity mask
    [roiObj_Int.data] = rangeReSeg(volObj.data,roiObj_Int.data,range);
    [roiObj_Int.data] = outlierReSeg(volObj.data,roiObj_Int.data,outliers);
    
    % Initialization
    volInt_RE = roiExtract(volObj.data,roiObj_Int.data);
    
    try
        % STEP 3: CALCULATION OF MORPHOLOGICAL FEATURES
        radiomics.morph = getMorphFeatures(volObj.data,roiObj_Int.data,roiObj_Morph.data,scaleNonText,intensity); % For scans with arbitrary units, some features will not be computed.
    catch
        fprintf('PROBLEM WITH COMPUTATION OF MORPHOLOGICAL FEATURES ')
        radiomics.morph = 'ERROR_COMPUTATION';
    end
    
    try
        % STEP 4: CALCULATION OF LOCAL INTENSITY FEATURES
        radiomics.locInt = getLocIntFeatures(volObj.data,roiObj_Int.data,scaleNonText,intensity); % For scans with arbitrary units, all of these features will not be computed.
    catch
        fprintf('PROBLEM WITH COMPUTATION OF LOCAL INTENSITY FEATURES ')
        radiomics.locInt = 'ERROR_COMPUTATION';        
    end
    
    try
        % STEP 5: CALCULATION OF STATISTICAL FEATURES
        radiomics.stats = getStatsFeatures(volInt_RE,intensity); % For scans with arbitrary units, some features will not be computed.
    catch
        fprintf('PROBLEM WITH COMPUTATION OF STATISTICAL FEATURES ')
        radiomics.stats = 'ERROR_COMPUTATION';           
    end
    
    try
        % STEP 6: CALCULATION OF INTENSITY HISTOGRAM FEATURES
        [volQuant_RE] = discretisation(volInt_RE,IH.type,IH.val,userSetMinVal); % There would actually be no need to include "userSetMinVal" here as fourth argument, as discretisation for IH features is (logically) always set to "FBN". This value will not be used for "FBN", see discretisation.m code.
        radiomics.intHist = getIntHistFeatures(volQuant_RE);
    catch
        fprintf('PROBLEM WITH COMPUTATION OF INTENSITY HISTOGRAM FEATURES ')
        radiomics.intHist = 'ERROR_COMPUTATION';           
    end
    
    try
        % STEP 7: CALCULATION OF INTENSITY-VOLUME HISTOGRAM FEATURES
        if ~isempty(IVH), [volQuant_RE,wb] = discretisation(volInt_RE,IVH.type,IVH.val,userSetMinVal,'ivh'); else volQuant_RE = volInt_RE; wb = 1; end % FOR CT, WE DO NOT WANT TO DISCRETISE. AN EMPTY IVH STRUCT ([]) DEFINES WHAT WE WANT TO USE FOR CT. FOR PET: FBS/0.1; FOR MRI: FBN/1000.
        if ~isempty(IVH)
            if strcmp(IVH.type,'FBS') || strcmp(IVH.type,'FBSequal')
                radiomics.intVolHist = getIntVolHistFeatures(volQuant_RE,wb,range);
            else
                radiomics.intVolHist = getIntVolHistFeatures(volQuant_RE,wb);
            end
        else
            radiomics.intVolHist = getIntVolHistFeatures(volQuant_RE,wb,range);
        end
    catch
        fprintf('PROBLEM WITH COMPUTATION OF INTENSITY-VOLUME HISTOGRAM FEATURES ')
        radiomics.intVolHist = 'ERROR_COMPUTATION';         
    end
    
    toc
catch
    fprintf('PROBLEM WITH PRE-PROCESSING OF NON-TEXTURE FEATURES')
    radiomics.morph = 'ERROR_PROCESSING'; radiomics.locInt = 'ERROR_PROCESSING'; radiomics.stats = 'ERROR_PROCESSING'; radiomics. intHist = 'ERROR_PROCESSING'; radiomics.intVolHist = 'ERROR_PROCESSING';
end
% -------------------------------------------------------------------------



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