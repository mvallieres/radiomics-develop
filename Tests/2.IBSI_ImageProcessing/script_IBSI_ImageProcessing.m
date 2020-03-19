% *************************************************************************************
% * DESCRIPTION:                                                                      *
% * This script computes the image processing standardization tests #2 for the IBSI   *
% * --------------------------------------------------------------------------------- *
% * AUTHOR(S): Martin Vallieres <mart.vallieres@gmail.com>                            *
% * --------------------------------------------------------------------------------- *
% * HISTORY:                                                                          *
% * - Creation: June 2018                                                             *
% * --------------------------------------------------------------------------------- *
% * DISCLAIMER:                                                                       *                                                             
% * "I'm not a programmer, I'm just a scientist doing stuff!"                         *
% * --------------------------------------------------------------------------------- *
% * STATEMENT:                                                                        * 
% * --> Copyright (C) 2018  Martin Vallieres                                          *
% *                                                                                   *
% * This program is written on the basis of a scientific collaboration                *
% * between the following parties:                                                    *
% * - Martin Vallieres <mart.vallieres@gmail.com>                                     *
% * - Alex Zwanenburg <alexander.zwanenburg@nct-dresden.de>                           *
% *                                                                                   *
% * By using this file, all parties acknowledge that it is to be kept private         *
% * until publication of a potential scientific study.                                *
% *************************************************************************************

clc,clear,fprintf('\n'), warning off, pathWORK = pwd;
if isempty(mfilename) % That means we are running the script from the command line.
    [scriptFileName,logFileName] = scriptCommandLine_FindNames(pathWORK);
else
    scriptFileName = [mfilename,'.m'];
end
help(scriptFileName)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          PARAMETER OPTIONS                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PARAMETERS - IBSI Case 3
imParamCT3.interp.scaleNonText = [2,2,2];
imParamCT3.interp.scaleText = {[2,2,2]};
imParamCT3.interp.volInterp = 'linear'; 
imParamCT3.interp.glRound = 1; 
imParamCT3.interp.roiInterp = 'linear'; 
imParamCT3.interp.roiPV = 0.5; 
imParamCT3.reSeg.range = [-1000,400];
imParamCT3.reSeg.outliers = ''; 
imParamCT3.discretisation.IH.type = 'FBS'; 
imParamCT3.discretisation.IH.val = 25; 
imParamCT3.discretisation.IVH.type = 'FBS'; 
imParamCT3.discretisation.IVH.val = 2.5; 
imParamCT3.discretisation.texture.type = {'FBS'}; 
imParamCT3.discretisation.texture.val = {[25]}; 
imParamCT3.intensity = 'definite'; imParamCT3.units = 'HU'; imParamCT3.type = 'CTscan';
imParamCT3.distCorrection = false;
imParamCT3.computeDiagFeatures = true;

% PARAMETERS - IBSI Case 4
imParamCT4.interp.scaleNonText = [2,2,2]; 
imParamCT4.interp.scaleText = {[2,2,2]};
imParamCT4.interp.volInterp = 'linear'; 
imParamCT4.interp.glRound = 1; 
imParamCT4.interp.roiInterp = 'linear'; 
imParamCT4.interp.roiPV = 0.5;
imParamCT4.reSeg.range = []; 
imParamCT4.reSeg.outliers = 'Collewet'; 
imParamCT4.discretisation.IH.type = 'FBN'; 
imParamCT4.discretisation.IH.val = 32; 
imParamCT4.discretisation.IVH = []; 
imParamCT4.discretisation.texture.type = {'FBN'}; 
imParamCT4.discretisation.texture.val = {[32]}; 
imParamCT4.intensity = 'definite'; imParamCT4.units = 'HU'; imParamCT4.type = 'CTscan';
imParamCT4.distCorrection = false;
imParamCT4.computeDiagFeatures = true;

% PARAMETERS - IBSI Case 5
imParamCT5.interp.scaleNonText = [2,2,2]; 
imParamCT5.interp.scaleText = {[2,2,2]}; 
imParamCT5.interp.volInterp = 'spline'; 
imParamCT5.interp.glRound = 1; 
imParamCT5.interp.roiInterp = 'linear'; 
imParamCT5.interp.roiPV = 0.5;
imParamCT5.reSeg.range = [-1000,400];
imParamCT5.reSeg.outliers = 'Collewet'; 
imParamCT5.discretisation.IH.type = 'FBN';
imParamCT5.discretisation.IH.val = 32; 
imParamCT5.discretisation.IVH.type = 'FBN'; 
imParamCT5.discretisation.IVH.val = 1000; 
imParamCT5.discretisation.texture.type = {'FBN'};
imParamCT5.discretisation.texture.val = {[32]};
imParamCT5.intensity = 'definite'; imParamCT5.units = 'HU'; imParamCT5.type = 'CTscan';
imParamCT5.distCorrection = false;
imParamCT5.computeDiagFeatures = true;

% ROI OPTIONS
roiTypes = {'GTVp','GTVp','GTVp'};
roiType_labels = {'Case3','Case4','Case5'}; 
nCases = numel(roiTypes);

% PARALLEL OPTIONS
nBatch = 1; 
matlabPATH = ['"',fullfile(matlabroot,'bin','matlab'),'"']; % Full path to the matlab executable on your system.
global codePATH
codePATH = '~/GitHub/radiomics-develop/Code'; % IMPORTANT OPTION: To be adapted to your system. Full path to the "radiomics-develop" code on your system.
addpath(genpath(codePATH));
% -------------------------------------------------------------------------






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          RADIOMICS COMPUTATION                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% STEP 1: READING DATA
cd(pathWORK)
if ispc
    system('move DICOM_IBSI DICOM');
else
    system('mv DICOM_IBSI DICOM');
end
pathDICOM = fullfile(pathWORK,'DICOM'); 
mkdir('DATA'), pathDATA = fullfile(pathWORK,'DATA');
cd(pathWORK)
if ispc
    system('move CSV_IBSI CSV');
else
    system('mv CSV_IBSI CSV');
end 
pathCSV = fullfile(pathWORK,'CSV');
mkdir('FEATURES'), pathFEATURES = fullfile(pathWORK,'FEATURES');
tic, fprintf('\n--> READING DATA USING %u CORES ... ',nBatch)
readAllOrganizedData_batch(pathDICOM,pathDATA,nBatch,matlabPATH) % If data is not organized and only in DICOM format, it is still possible to use the function "readAllDICOM.m'.
fprintf('DONE!\n'), toc


% STEP 2: COMPUTING RADIOMICS FEATURES
allImParams = {imParamCT3,imParamCT4,imParamCT5};
tic, fprintf('\n--> COMPUTING RADIOMIC FEATURES WITH %u CORES ... ',nBatch)
for c = 1:nCases
    roiType = {roiTypes{c}}; roiType_label = {roiType_labels{c}};
    imParams.CTscan.image = allImParams{c};
    computeRadiomics_batchAllPatients(pathDATA,pathCSV,pathFEATURES,imParams,roiType,roiType_label,nBatch,matlabPATH)
end
fprintf('\n'), toc


% STEP3: PRINTING RESULTS
nonText = {'diag','morph_3D','locInt_3D','stats_3D','intHist_3D','intVolHist_3D'}; nNonText = numel(nonText);
text = {'glcm_3Dmrg','glrlm_3Dmrg','glszm_3D','gldzm_3D','ngtdm_3D','ngldm_3D'}; nText = numel(text);
cases = [3,4,5];
for c = 1:nCases
    thisCase = num2str(cases(c));
    features = struct;
    cd(pathFEATURES), load(['Lung-IBSI-001__CT(Case',thisCase,').CTscan.mat'])
    for n = 1:nNonText
        tempStruct = radiomics.image.(nonText{n});
        fields = fieldnames(tempStruct);
        tempStruct = tempStruct.(fields{1});
        features = concatenateStruct(features,tempStruct);
    end
    for n = 1:nText
        tempStruct = radiomics.image.texture.(text{n});
        fields = fieldnames(tempStruct);
        tempStruct = tempStruct.(fields{1});
        features = concatenateStruct(features,tempStruct);
    end
    cd(pathWORK), printStructure(features,['IBSIresults_ImageProcessing_Case',thisCase]);
end
% -------------------------------------------------------------------------
