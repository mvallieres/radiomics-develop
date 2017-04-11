% ******************************************************************************
% * DESCRIPTION:                                                               *
% * This script computes all radiomic features. The script "script_ReadData.m" *
% * must have been run prior to this one.                                      *                    
% * -------------------------------------------------------------------------- *
% * AUTHOR(S):                                                                 *
% * - Martin Vallieres <mart.vallieres@gmail.com>                              *
% * -------------------------------------------------------------------------- *
% * HISTORY:                                                                   *
% * - Creation: April 2017                                                     * 
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

clc,clear,fprintf('\n'), warning off
help script_ComputeRadiomics, pathWORK = pwd;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          PARAMETER OPTIONS                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% RADIOMIC PARAMETERS OPTIONS
imParamMR.interp.scaleNonText = [3,3,3]; % Resolution in mm for the computation of non-texture features. Here, given the large difference between in-plane resolution (~ 1 mm) and slice thickness (~ 5 mm), we take some middle point.
imParamMR.interp.scaleText = {[1,1,1],[2,2,2],[3,3,3],[4,4,4]}; % Different resolutions in mm tested for the computation of texture features. Each cell entry is a [X,Y,Z] resolution in MATLAB world coordinate.
imParamMR.interp.volInterp = 'linear'; % Using linear interpolation for imaging intensities. Conservative choice to prevent interpolation craziness.
imParamMR.interp.glRound = []; % No grey-level rounding
imParamMR.interp.roiInterp = 'linear'; % Using linear interpolation for ROI mask. Conservative choice to prevent interpolation craziness.
imParamMR.interp.roiPV = 0.5; % After interpolation, a value >=0.5 is assigned to a 1, and <0.5 to a 0.
imParamMR.reSeg.range = []; % No range re-segmentation is performed for MR.
imParamMR.reSeg.outliers = 'Collewet'; % Using Collewet normalization to remove outliers for MRI.
imParamMR.discretisation.IH.type = 'FBN'; % Using fixed bin number for intensity-histogram features.
imParamMR.discretisation.IH.val = 64; % Using 64 grey-levels for intensity-histogram features.
imParamMR.discretisation.IVH.type = 'FBN'; % Using fixed bin number for intensity-volume histogram features.
imParamMR.discretisation.IVH.val = 1000; % Using 1000 grey-levels for intensity-volume histogram features.
imParamMR.discretisation.texture.type = {'FBN','FBNequal'}; % Using two different types of quantization algorithms (always fixed bin number for MR for now). If "FBS" or "FBSequal" is used, imParam.reSeg.range must be defined, as userSetMinVal will be assigned to imParam.reSeg.range(1).
imParamMR.discretisation.texture.val = {[8,16,32,64],[8,16,32,64]}; % Gray-levels to test for each algorithm on the above line (definition depends on the algorithm). The total number must be the same in each cell entry.

imParamCT.interp.scaleNonText = [2,2,2]; % Resolution in mm for the computation of non-texture features. Here, given the difference between in-plane resolution (~ 1 mm) and slice thickness (~ 3 mm), we take some middle point.
imParamCT.interp.scaleText = {[1,1,1],[2,2,2],[3,3,3],[4,4,4]}; % Different resolutions in mm tested for the computation of texture features. Each cell entry is a [X,Y,Z] resolution in MATLAB world coordinate.
imParamCT.interp.volInterp = 'linear'; % Using cubic interpolation for imaging intensities. Conservative choice to prevent interpolation craziness.
imParamCT.interp.glRound = 1; % Grey-level rounding to 1 HU.
imParamCT.interp.roiInterp = 'linear'; % Using cubic interpolation for ROI mask. Conservative choice to prevent interpolation craziness.
imParamCT.interp.roiPV = 0.5; % After interpolation, a value >=0.5 is assigned to a 1, and <0.5 to a 0.
imParamCT.reSeg.range = [-500,400]; % Voxels within the ROI with HU outside that range are discarded. This can be adjusted depending on tumour type. For now, it covers a range of HU going from lungs to dense soft-tissues. If you are only working with soft-tissue sarcomas for example, perhaps a range of [-100,400] would be more meaningful. The minimal value is also of special importance for FBS discretisation. 
imParamCT.reSeg.outliers = ''; % Not using outlier removal for CT (for now).
imParamCT.discretisation.IH.type = 'FBN'; % Using fixed bin number for intensity-histogram features.
imParamCT.discretisation.IH.val = 64; % Using 64 grey-levels for intensity-histogram features.
% imParamCT.discretisation.IVH = []; % No need to discretise for CT, natural integer discretisation due to HU. Using "[]" is similar to FBS with a bin width of 1, but without a minimum value set up like in the case of PET (SUV = 0). But shall we also use a set minimal value (e.g. -500 HU as defined in imParam.reSgeg.range(1)) for CT? The answer is yes, so this line is commented for now, and new lines are added below.
imParamCT.discretisation.IVH.type = 'FBS'; % Using fixed bin size for intensity-volume histogram features with CT. imParam.reSeg.range must be defined, as userSetMinVal will be assigned to imParam.reSeg.range(1).
imParamCT.discretisation.IVH.val = 1; % Using 1 HU bin width for intensity-histogram features (natural HU discretisation).
imParamCT.discretisation.texture.type = {'FBN','FBNequal','FBS','FBSequal'}; % Using four different types of quantization algorithms. If "FBS" or "FBSequal" is used, imParam.reSeg.range must be defined, as userSetMinVal will be assigned to imParam.reSeg.range(1).
imParamCT.discretisation.texture.val = {[8,16,32,64],[8,16,32,64],[12.5,25,50,100],[12.5,25,50,100]}; % Gray-levels to test for each algorithm on the above line (definition depends on the algorithm). The total number must be the same in each cell entry.

imParamPET.interp.scaleNonText = [4,4,4]; % Resolution in mm for the computation of non-texture features. 
imParamPET.interp.scaleText = {[1,1,1],[2,2,2],[3,3,3],[4,4,4]}; % Different resolutions in mm tested for the computation of texture features. Each cell entry is a [X,Y,Z] resolution in MATLAB world coordinate.
imParamPET.interp.volInterp = 'linear'; % Using cubic interpolation for imaging intensities. Conservative choice to prevent interpolation craziness.
imParamPET.interp.glRound = []; % No grey-level rounding, as everything is computed from continuous SUV maps.
imParamPET.interp.roiInterp = 'linear'; % Using cubic interpolation for ROI mask. Conservative choice to prevent interpolation craziness.
imParamPET.interp.roiPV = 0.5; % After interpolation, a value >=0.5 is assigned to a 1, and <0.5 to a 0.
imParamPET.reSeg.range = [0,inf]; % We are working with SUV maps, going from 0 to infinity. The minimal value is of special importance for FBS discretisation.
imParamPET.reSeg.outliers = ''; % Not using outlier removal for PET (for now).
imParamPET.discretisation.IH.type = 'FBN'; % Using fixed bin number for intensity-histogram features.
imParamPET.discretisation.IH.val = 64; % Using 64 grey-levels for intensity-histogram features.
imParamPET.discretisation.IVH.type = 'FBS'; % Using fixed bin size for intensity-volume histogram features with PET. imParam.reSeg.range must be defined, as userSetMinVal will be assigned to imParam.reSeg.range(1).
imParamPET.discretisation.IVH.val = 0.1; % Using 0.1 SUV bin width for intensity-histogram features.
imParamPET.discretisation.texture.type = {'FBN','FBNequal','FBS','FBSequal'}; % Using four different types of quantization algorithms. If "FBS" or "FBSequal" is used, imParam.reSeg.range must be defined, as userSetMinVal will be assigned to imParam.reSeg.range(1).
imParamPET.discretisation.texture.val = {[8,16,32,64],[8,16,32,64],[0.25,0.5,1,2],[0.25,0.5,1,2]}; % Gray-levels to test for each algorithm on the above line (definition depends on the algorithm). The total number must be the same in each cell entry.

imParam.MRscan = imParamMR;
imParam.CTscan = imParamCT;
imParam.PTscan = imParamPET;

% ROI OPTIONS
roiType = 'GTV'; % Name you want to give to the type of ROI being analyzed, for file saving purposes. This script currently computes radiomic features for one ROI at a time.
manualROIchoice = false; % If set to true, the user will be prompted to choose which ROI(s) to use for a given imaging scan, so no need to create a roiNames.mat file --> this option is convenient for small datasets. If set to false, the algorithm will use the roiNames.mat file present in WORKSPACE --> this option is convenienent for large datasets, as manual reading and choice for all imaging scans in this script will be long.

% PARALLEL OPTIONS
nBatch = 4; % To compute radiomic features using parallelization. Put this number to the total number of cores that you have on your machine.
matlabPATH = 'matlab'; % IMPORTANT: Full path to the matlab executable on the system. --> Ex: '/home/martin/Programs/MATLAB/R2016a/bin/matlab'. Here, a symbolic link to the full MATLAB path has previously been created on Martin Vallieres' computer. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    RADIOMIC FEATURE COMPUTATION CODE                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pathDATA = fullfile(pathWORK,'DATA');
mkdir('FEATURES'), pathFEATURES = fullfile(pathWORK,'FEATURES'); cd(pathWORK)

tStart = tic;
fprintf('\n\n************************* RADIOMIC FEATURE COMPUTATION *************************')

% 1. ROI CHOICE (optional, see above) OR LOADING PREVIOUSLY CONSTRUCTED roiNames.mat
if manualROIchoice
   roiNames = roiChoice(pathWORK,pathDATA); % User will be prompted to choose each ROI for each scan of each patient. Warning: at the end of the function, a new roiNames.mat files will be saved in the WORKSPACE, and potentially already existing ones will be given another name.
end
cd(pathWORK), load('roiNames') % Variable roiNames is now in MATLAB workspace.

% 2. COMPUTING RADIOMIC FEATURES 
tic, fprintf('\n--> COMPUTING RADIOMIC FEATURES WITH %u CORES ... ',nBatch)
computeAllRadiomics_batch(pathDATA,pathFEATURES,roiNames,imParam,roiType,nBatch,matlabPATH)
fprintf('DONE!\n'), toc

time = toc(tStart);
fprintf('\n\nTOTAL TIME FOR RADIOMIC FEATURE COMPUTATION: %f seconds\n',time)
fprintf('-------------------------------------------------------------------------------------')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         MULTIVARIABLE MODELING CODE IN ANOTHER SCRIPT (to come)         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --> roiNames.mat in /WORKSPACE and getRadiomicNames.m can now be used   %
%     to read radiomics data.                                             %
% --> Subsequent organization of data matrix [nPatient X nFeature] should %
%     not be so difficult.                                                %
% --> With some outcome vector, multivariable analysis can start.         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
