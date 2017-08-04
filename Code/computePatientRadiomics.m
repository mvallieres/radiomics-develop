function  [radiomics] = computePatientRadiomics(sData,nameROI,nameSet,imParam)
% -------------------------------------------------------------------------
% AUTHOR(S): 
% - Martin Vallieres <mart.vallieres@gmail.com>
% -------------------------------------------------------------------------
% HISTORY:
% - Creation: April 2017
% - Revision I: August 2017
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


% INITIALIZATION
radiomics = struct;
scaleNonText = imParam.interp.scaleNonText; 
volInterp = imParam.interp.volInterp; roiInterp = imParam.interp.roiInterp;
glRound = imParam.interp.glRound;
roiPV = imParam.interp.roiPV;
range = imParam.reSeg.range;
outliers = imParam.reSeg.outliers;
IH = imParam.discretisation.IH;
IVH = imParam.discretisation.IVH;
scaleText = imParam.interp.scaleText; 
algo = imParam.discretisation.texture.type; 
grayLevels = imParam.discretisation.texture.val;
nScale = numel(scaleText); nAlgo = numel(algo); nGl = numel(grayLevels{1}); nExp = nScale*nAlgo*nGl;
glcm = cell(nScale,nAlgo,nGl); glrlm = cell(nScale,nAlgo,nGl); glszm = cell(nScale,nAlgo,nGl); ngtdm = cell(nScale,nAlgo,nGl); gldzm = cell(nScale,nAlgo,nGl); ngldm = cell(nScale,nAlgo,nGl);


% SETTING UP userSetMinVal
if ~isempty(range)
    userSetMinVal = range(1);
else
    if strcmp(sData{2}.type,'PTscan')
        userSetMinVal = 0; % SUV is defined from 0 to inf.
    elseif strcmp(sData{2}.type,'CTscan')
        userSetMinVal = -1000; % No assumption on range of HU.
    elseif strcmp(sData{2}.type,'MRscan')
        userSetMinVal = 0; % This a dummy value, has no consequence on the rest of the code since FBS discretisation is not used for MRI.
    end
end

% COMPUTATION OF ROI MASK
tic, fprintf('\n--> Computation of ROI mask: ')
box = 'box10'; % 10 voxels in all three dimensions are added to the smallest bounding box.
try
    contourNumber = findContour(sData,nameROI,nameSet);
    [volObjInit,roiObjInit] = getROI(sData,contourNumber,box); % This takes care of the "Volume resection" step as well using the argument "box". No fourth argument: 'interp' by default.
catch
    fprintf('\nPROBLEM WITH ROI')
    radiomics.morph = []; radiomics.locInt = []; radiomics.stats = []; radiomics. intHist = []; radiomics.intVolHist = []; radiomics.glcm = []; radiomics.glrlm = []; radiomics.glszm = []; radiomics.ngtdm = []; radiomics.gldzm = []; radiomics.ngldm = [];
    radiomics.imParam = imParam;
    return
end
toc


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
        radiomics.morph = getMorphFeatures(volObj.data,roiObj_Int.data,roiObj_Morph.data,scaleNonText);
    catch
        fprintf('PROBLEM WITH COMPUTATION OF MORPHOLOGICAL FEATURES ')
        radiomics.morph = [];
    end
    
    try
        % STEP 4: CALCULATION OF LOCAL INTENSITY FEATURES
        radiomics.locInt = getLocIntFeatures(volObj.data,roiObj_Int.data,scaleNonText);
    catch
        fprintf('PROBLEM WITH COMPUTATION OF LOCAL INTENSITY FEATURES ')
        radiomics.locInt = [];        
    end
    
    try
        % STEP 5: CALCULATION OF STATISTICAL FEATURES
        radiomics.stats = getStatsFeatures(volInt_RE);
    catch
        fprintf('PROBLEM WITH COMPUTATION OF STATISTICAL FEATURES ')
        radiomics.stats = [];           
    end
    
    try
        % STEP 6: CALCULATION OF INTENSITY HISTOGRAM FEATURES
        [volQuant_RE] = discretisation(volInt_RE,IH.type,IH.val,userSetMinVal); % There would actually be no need to include "userSetMinVal" here as fourth argument, as discretisation for IH features is (logically) always set to "FBN". This value will not be used for "FBN", see discretisation.m code.
        radiomics.intHist = getIntHistFeatures(volQuant_RE);
    catch
        fprintf('PROBLEM WITH COMPUTATION OF INTENSITY HISTOGRAM FEATURES ')
        radiomics.intHist = [];           
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
        radiomics.intVolHist = [];         
    end
    
    toc
catch
    fprintf('PROBLEM WITH PRE-PROCESSING OF NON-TEXTURE FEATURES')
    radiomics.morph = []; radiomics.locInt = []; radiomics.stats = []; radiomics. intHist = []; radiomics.intVolHist = [];
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
                        glcm{s,a,n} = [];
                    end
                    try
                        glrlm{s,a,n} = getGLRLMfeatures(volQuant_RE);
                    catch
                        fprintf('PROBLEM WITH COMPUTATION OF GLRLM FEATURES ')
                        glrlm{s,a,n} = [];                        
                    end
                    try
                        glszm{s,a,n} = getGLSZMfeatures(volQuant_RE);
                    catch
                        fprintf('PROBLEM WITH COMPUTATION OF GLSZM FEATURES ')
                        glszm{s,a,n} = [];                         
                    end
                    try
                        gldzm{s,a,n} = getGLDZMfeatures(volQuant_RE,roiObj_Morph.data);
                    catch
                        fprintf('PROBLEM WITH COMPUTATION OF GLDZM FEATURES ')
                        gldzm{s,a,n} = [];                            
                    end
                    try
                        ngtdm{s,a,n} = getNGTDMfeatures(volQuant_RE);
                    catch
                        fprintf('PROBLEM WITH COMPUTATION OF NGTDM FEATURES ')
                        ngtdm{s,a,n} = [];                           
                    end
                    try
                        ngldm{s,a,n} = getNGLDMfeatures(volQuant_RE);
                    catch
                        fprintf('PROBLEM WITH COMPUTATION OF NGLDM FEATURES ')
                        ngldm{s,a,n} = [];                          
                    end

                    toc
                catch
                    fprintf('PROBLEM WITH DISCRETISATION')
                    glcm{s,a,n} = []; glrlm{s,a,n} = []; glszm{s,a,n} = []; ngtdm{s,a,n} = []; gldzm{s,a,n} = []; ngldm{s,a,n} = [];
                end
            end
        end
    catch 
        fprintf('PROBLEM WITH PRE-PROCESSING OF TEXTURE FEATURES')
        for a = 1:nAlgo
            for n = 1:nGl
                count = count + 1;
                glcm{s,a,n} = []; glrlm{s,a,n} = []; glszm{s,a,n} = []; ngtdm{s,a,n} = []; gldzm{s,a,n} = []; ngldm{s,a,n} = [];
            end
        end
    end
end

radiomics.glcm = glcm;
radiomics.glrlm = glrlm;
radiomics.glszm = glszm;
radiomics.ngtdm = ngtdm;
radiomics.gldzm = gldzm;
radiomics.ngldm = ngldm;
% -------------------------------------------------------------------------

radiomics.imParam = imParam;
end