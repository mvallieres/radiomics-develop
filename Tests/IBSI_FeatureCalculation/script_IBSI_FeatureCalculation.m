% ******************************************************************************
% * DESCRIPTION:                                                               *
% * This script computes performs the Image Biomarker Standardisation          * 
% * Initiative (IBSI) calibration test -- Phase 1: Feature calculation using   *
% * the digital phantom.                                                       * 
% * -------------------------------------------------------------------------- *
% * AUTHOR(S):                                                                 *
% * - Martin Vallieres <mart.vallieres@gmail.com>                              *
% * -------------------------------------------------------------------------- *
% * DISCLAIMER:                                                                *                                                             
% * "I'm not a programmer, I'm just a scientist doing stuff!"                  *
% * -------------------------------------------------------------------------  *
% * STATEMENT:                                                                 *
% * This file is part of <https://github.com/mvallieres/radiomics-develop/>,   *  
% * a private repository dedicated to the development of programming code for  *
% * new radiomics applications.                                                * 
% * --> Copyright (C) 2017  Martin Vallieres                                   *
% *     All rights reserved.                                                   *
% *                                                                            *
% * This file is written on the basis of a scientific collaboration for the    * 
% * "radiomics-develop" team.                                                  *
% *                                                                            * 
% * By using this file, all members of the team acknowledge that it is to be   * 
% * kept private until public release. Other scientists willing to join the    * 
% * "radiomics-develop" team is however highly encouraged. Please contact      *
% * Martin Vallieres for this matter.                                          *  
% * -------------------------------------------------------------------------- *
% ******************************************************************************

tic, clc,clear,fprintf('\n'), warning off, pathWORK = pwd;
if isempty(mfilename) % That means we are running the script from the command line using the following format (mandatory): matlab -nodisplay -nodesktop -nosplash -singleCompThread < script_IBSI_FeatureCalculation.m >& script_IBSI_FeatureCalculation.log &
    [scriptFileName,logFileName] = scriptCommandLine_FindNames(pathWORK);
else
    scriptFileName = [mfilename,'.m'];
end
help(scriptFileName)
fprintf('--> COMPUTING IBSI CALIBRATION TEST -- PHASE 1 (Digital Phantom) ... ')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   LOADING VARIABLES USED FOR THE IBSI TESTS USING THE DIGITAL PHANTOM   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd(pathWORK), run testVolume % Variable 'imgObj' and 'roiObj' gets out of there.
maskInt = roiObj; maskMorph = roiObj; % Both "intensity" and "morphological" masks are the same in this example. Both are created for the sake of the exercise.
roiOnlyInt = imgObj; roiOnlyInt(maskInt == 0) = NaN; % Imaging data with NaNs outside the mask.
roiOnlyMorph = imgObj; roiOnlyMorph(maskMorph == 0) = NaN; % Imaging data with NaNs outside the mask.
% ------------------------------------------------------------------------%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              OPTIONS                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
res = [2,2,2]; % XYZ resolution (world), or JIK resolution (intrinsic matlab). The test volume is 2 X 2 X 2 mm^3
distCorrection = false; % Set this variable to true in order to use discretization length difference corrections as used here: https://doi.org/10.1088/0031-9155/60/14/5471. Set this variable to false to replicate IBSI results.
% ------------------------------------------------------------------------%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        FEATURE COMPUTATION                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
features = struct; % Initialization of final structure containing all features.

% CATEGORY 1: MORPHOLOGICAL FEATURES
featuresMorph = getMorphFeatures(imgObj,maskInt,maskMorph,res);
features = concatenateStruct(features,featuresMorph);


% CATEGORY 2: LOCAL INTENSITY FEATURES
featuresLocInt = getLocIntFeatures(imgObj,maskInt,res);
features = concatenateStruct(features,featuresLocInt);


% CATEGORY 3: FOR STATS FEATURES
featuresStats = getStatsFeatures(roiOnlyInt);
features = concatenateStruct(features,featuresStats);


% CATEGORY 4: INTENSITY HISTOGRAM FEATURES
featuresIntHist = getIntHistFeatures(roiOnlyInt);
features = concatenateStruct(features,featuresIntHist);


% CATEGORY 5: INTENSITY-VOLUME HISTOGRAM FEATURES
featuresIntVolHist = getIntVolHistFeatures(roiOnlyInt,1,[1,6]); % Second and third arguments are specific to the IBSI digital phantom.
features = concatenateStruct(features,featuresIntVolHist);


% CATEGORY 6: GLCM FEATURES
featuresGLCM = getGLCMfeatures(roiOnlyInt,distCorrection);
features = concatenateStruct(features,featuresGLCM);


% CATEGORY 7: GLRLM FEATURES
featuresGLRLM = getGLRLMfeatures(roiOnlyInt,distCorrection);
features = concatenateStruct(features,featuresGLRLM);


% CATEGORY 8: GLSZM FEATURES
featuresGLSZM = getGLSZMfeatures(roiOnlyInt);
features = concatenateStruct(features,featuresGLSZM);


% CATEGORY 9: GLDZM FEATURES
featuresGLDZM = getGLDZMfeatures(roiOnlyInt,maskMorph);
features = concatenateStruct(features,featuresGLDZM);


% CATEGORY 10: NGTDM FEATURES
featuresNGTDM = getNGTDMfeatures(roiOnlyInt,distCorrection);
features = concatenateStruct(features,featuresNGTDM);


% CATEGORY 11: NGLDM FEATURES
featuresNGLDM = getNGLDMfeatures(roiOnlyInt);
features = concatenateStruct(features,featuresNGLDM);
% ------------------------------------------------------------------------%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           PRINTING RESULTS                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd(pathWORK), printStructure(features,'IBSIresults_FeatureCalculation');
% ------------------------------------------------------------------------%


fprintf('DONE!\n')
toc
