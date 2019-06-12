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
% * - Revision I: August 2017                                                  *                                              
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

clc,clear,fprintf('\n'), warning off, pathWORK = pwd;
if isempty(mfilename) % That means we are running the script from the command line using the following format (mandatory): matlab -nodisplay -nodesktop -nosplash -singleCompThread < script_ComputeRadiomics.m >& script_ComputeRadiomics.log &
    [scriptFileName,logFileName] = scriptCommandLine_FindNames(pathWORK);
else
    scriptFileName = [mfilename,'.m'];
end
help(scriptFileName)
timeStamp = char(datetime('now','Format','yyyy-MM-dd_HH:mm:ss'));
softwareLabel = 'PiCare-Matlab.Radiomics';
softwareVersion = '0.1';
programmingLanguage = 'Matlab2018a';
institution = 'McGill_UCSF_D-Lab_Oncoray';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          PARAMETER OPTIONS                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Change parameters where you see the word "OPTION" in the comment.
% NOTES: - All scale parameters below have the following format: [Xin,Yin,Zslice], where Xin and Yin are the X (left to right) and Y (bottom to top) IN-PLANE resolutions, and Zslice is the slice spacing, NO MATTER THE ORIENTATION OF THE VOLUME (i.e. axial , sagittal, coronal).   
%        - If interpolation in the slice direction (no matter the orientation) is not wanted (i.e. 2D case), please used vectors with only two entries [Xin,Yin] for the scale parameters below (e.g. scale = [2,2]).
%        - If no interpolation at all is wanted, use 0 as the scale parameter (e.g. scale = 0). It will be recorded as such in the final radiomics structure for differentiation purposes.

% RADIOMIC PARAMETERS OPTIONS ---------------------------------------------

% 1. MRI PARAMETERS
imParamMR.interp.scaleNonText = [3,3,3]; % OPTION: Resolution in mm for the computation of non-texture features. Here, given the large difference between commonly seen in-plane resolution (~ 1 mm) and slice thickness (~ 5 mm), we take some middle point.
imParamMR.interp.scaleText = {[1,1,1],[2,2,2],[3,3,3],[4,4,4]}; % OPTION: Different resolutions in mm tested for the computation of texture features. Each cell entry is a [X,Y,Z] resolution in MATLAB world coordinate.
imParamMR.interp.volInterp = 'linear'; % OPTION: Using linear interpolation for imaging intensities. Conservative choice to prevent interpolation craziness. IMPORTANT: ONLY OPTION AVAILABLE FOR NOW.
imParamMR.interp.glRound = []; % OPTION: No grey-level rounding
imParamMR.interp.roiInterp = 'linear'; % OPTION: Using linear interpolation for ROI mask. Conservative choice to prevent interpolation craziness. IMPORTANT: ONLY OPTION AVAILABLE FOR NOW.
imParamMR.interp.roiPV = 0.5; % OPTION: After interpolation, a value >=0.5 is assigned to a 1, and <0.5 to a 0.
imParamMR.reSeg.range = []; % OPTION: No range re-segmentation is performed for MR.
imParamMR.reSeg.outliers = 'Collewet'; % OPTION: Using Collewet normalization to remove outliers for MRI.
imParamMR.discretisation.IH.type = 'FBN'; % OPTION: Using fixed bin number for intensity-histogram features.
imParamMR.discretisation.IH.val = 64; % OPTION: Using 64 grey-levels for intensity-histogram features.
imParamMR.discretisation.IVH.type = 'FBN'; % OPTION: Using fixed bin number for intensity-volume histogram features.
imParamMR.discretisation.IVH.val = 1000; % OPTION: Using 1000 grey-levels for intensity-volume histogram features.
imParamMR.discretisation.texture.type = {'FBN','FBNequal'}; % OPTION: Using two different types of quantization algorithms (always fixed bin number for MR for now). If "FBS" or "FBSequal" is used, the minimum value of the ROI will be used in FBS discretisation since the lower bound of the re-segmentation range cannot be defined for arbitrary MRI intensity units (not recommended to use "FBS" or "FBSequal").
imParamMR.discretisation.texture.val = {[8,16,32,64],[8,16,32,64]}; % OPTION: Gray-levels to test for each algorithm on the above line (definition depends on the algorithm). The total number must be the same in each cell entry.
intensity_MR = 'arbitrary'; % OPTION: MRI is arbitrary intensity units. Use 'arbitrary' to not compute some features that depends on absolute intensity values. Use 'definite' to compute them all anyway.

% 2. CT PARAMETERS
imParamCT.interp.scaleNonText = [2,2,2]; % OPTION: Resolution in mm for the computation of non-texture features. Here, given the difference between commonly seen in-plane resolution (~ 1 mm) and slice thickness (~ 3 mm), we take some middle point.
imParamCT.interp.scaleText = {[1,1,1],[2,2,2],[3,3,3],[4,4,4]}; % OPTION: Different resolutions in mm tested for the computation of texture features. Each cell entry is a [X,Y,Z] resolution in MATLAB world coordinate.
imParamCT.interp.volInterp = 'linear'; % OPTION: Using linear interpolation for imaging intensities. Conservative choice to prevent interpolation craziness. IMPORTANT: ONLY OPTION AVAILABLE FOR NOW.
imParamCT.interp.glRound = 1; % OPTION: Grey-level rounding to 1 HU.
imParamCT.interp.roiInterp = 'linear'; % OPTION: Using linear interpolation for ROI mask. Conservative choice to prevent interpolation craziness. IMPORTANT: ONLY OPTION AVAILABLE FOR NOW.
imParamCT.interp.roiPV = 0.5; % OPTION: After interpolation, a value >=0.5 is assigned to a 1, and <0.5 to a 0.
imParamCT.reSeg.range = [-500,400]; % OPTION: Voxels within the ROI with HU outside that range are discarded. This can be adjusted depending on tumour type. For now, it covers a range of HU going from lungs to dense soft-tissues. If you are only working with soft-tissue sarcomas for example, perhaps a range of [-100,400] would be more meaningful. The minimal value is also of special importance for FBS discretisation. If the lower bound of the re-segmentation range is not defined, the minimum value of the ROI will be used in FBS discretisation (not recommended). 
imParamCT.reSeg.outliers = ''; % OPTION: Not using outlier removal for CT (for now).
imParamCT.discretisation.IH.type = 'FBN'; % OPTION: Using fixed bin number for intensity-histogram features.
imParamCT.discretisation.IH.val = 64; % OPTION: Using 64 grey-levels for intensity-histogram features.
imParamCT.discretisation.IVH = []; % OPTION: No need to discretise for CT, natural integer discretisation due to HU. Using "[]" is similar to FBS with a bin width of 1, but without a minimum value set up like in the case of PET (SUV = 0). But shall we also use a set minimal value (e.g. -500 HU as defined in imParam.reSgeg.range(1)) for CT for IVH features as well? If yes, use the two lines below instead (commented for now).
imParamCT.discretisation.texture.type = {'FBN','FBNequal','FBS','FBSequal'}; % OPTION: Using four different types of quantization algorithms. If "FBS" or "FBSequal" is used, imParam.reSeg.range must be defined, as userSetMinVal will be assigned to imParam.reSeg.range(1).
imParamCT.discretisation.texture.val = {[8,16,32,64],[8,16,32,64],[12.5,25,50,100],[12.5,25,50,100]}; % OPTION: Gray-levels to test for each algorithm on the above line (definition depends on the algorithm). The total number must be the same in each cell entry. For the FBS algorithm, we assume values are in units of HU.
intensity_CT = 'definite'; % OPTION: Definite intensity units (HU). 

% 3. PET PARAMETERS (assuming we always analyze SUV maps)
imParamPET.interp.scaleNonText = [4,4,4]; % OPTION: Resolution in mm for the computation of non-texture features. Common PET resolution is close to that.
imParamPET.interp.scaleText = {[1,1,1],[2,2,2],[3,3,3],[4,4,4]}; % OPTION: Different resolutions in mm tested for the computation of texture features. Each cell entry is a [X,Y,Z] resolution in MATLAB world coordinate.
imParamPET.interp.volInterp = 'linear'; % OPTION: Using linear interpolation for imaging intensities. Conservative choice to prevent interpolation craziness. IMPORTANT: ONLY OPTION AVAILABLE FOR NOW.
imParamPET.interp.glRound = []; % OPTION: No grey-level rounding, as everything is computed from continuous SUV maps.
imParamPET.interp.roiInterp = 'linear'; % OPTION: Using linear interpolation for ROI mask. Conservative choice to prevent interpolation craziness. IMPORTANT: ONLY OPTION AVAILABLE FOR NOW.
imParamPET.interp.roiPV = 0.5; % OPTION: After interpolation, a value >=0.5 is assigned to a 1, and <0.5 to a 0.
imParamPET.reSeg.range = [0,inf]; % OPTION: We are working with SUV maps, going from 0 to infinity. The minimal value is of special importance for FBS discretisation. If the lower bound of the re-segmentation range is not defined, the minimum value of the ROI will be used in FBS discretisation (not recommended).
imParamPET.reSeg.outliers = ''; % OPTION: Not using outlier removal for PET (for now).
imParamPET.discretisation.IH.type = 'FBN'; % OPTION: Using fixed bin number for intensity-histogram features.
imParamPET.discretisation.IH.val = 64; % OPTION: Using 64 grey-levels for intensity-histogram features.
imParamPET.discretisation.IVH.type = 'FBS'; % OPTION: Using fixed bin size for intensity-volume histogram features with PET. imParam.reSeg.range must be defined, as userSetMinVal will be assigned to imParam.reSeg.range(1).
imParamPET.discretisation.IVH.val = 0.1; % OPTION: Using 0.1 SUV bin width for intensity-histogram features.
imParamPET.discretisation.texture.type = {'FBN','FBNequal','FBS','FBSequal'}; % OPTION: Using four different types of quantization algorithms. If "FBS" or "FBSequal" is used, imParam.reSeg.range must be defined, as userSetMinVal will be assigned to imParam.reSeg.range(1).
imParamPET.discretisation.texture.val = {[8,16,32,64],[8,16,32,64],[0.25,0.5,1,2],[0.25,0.5,1,2]}; % OPTION: Gray-levels to test for each algorithm on the above line (definition depends on the algorithm). The total number must be the same in each cell entry. For the FBS algorithm, we assume values are in units of SUV.
intensity_PET = 'definite'; % OPTION: Definite intensity units (SUV).

% FILTER PARAMETERS (considered as arbitrary intensities). No FBS algorithm should ever be used at the moment, or otherwise the code will fail (TO SOLVE TO ALLOW THE USE OF THE MINIMUM VALUE OF THE ROI?).
imParamFilter.discretisation.IH.type = 'FBN'; % Using fixed bin number for intensity-histogram features.
imParamFilter.discretisation.IH.val = 64; % OPTION: Using 64 grey-levels for intensity-histogram features.
imParamFilter.discretisation.IVH.type = 'FBN'; % Using fixed bin number for intensity-volume histogram features.
imParamFilter.discretisation.IVH.val = 1000; % OPTION: Using 1000 grey-levels for intensity-volume histogram features.
imParamFilter.discretisation.texture.type = {'FBN','FBNequal'}; % Using two different types of quantization algorithms (always fixed bin number for MR for now).
imParamFilter.discretisation.texture.val = {[8,16,32,64],[8,16,32,64]}; % OPTION: Gray-levels to test for each algorithm on the above line (definition depends on the algorithm). The total number must be the same in each cell entry.
intensity_filter = 'filter'; % OPTION: Arbitrary intensity units. Use 'filter' to not compute some features that depends on absolute intensity values. Use 'definite' to compute them all anyway.

% WAVELET FILTERS OPTIONS 
computeWavelet = true; % OPTION: Set to true if you want to compute radiomic features on wavelet-flitered images. This significantly increases computation time. Set to false otherwise.
waveletName = 'coif1'; % OPTION: MATLAB name for the wavelet basis function used in wavelet decomposition.
% -------------------------------------------------------------------------

% ROI OPTIONS
roiTypes = {'undefined','GTVp','GTVp'}; % Name specifying the ROI type as defined as defined by the Radiomics Ontology (RO) effort (see the RO instructions for more details). For the moment, must be set to either one of these character arrays: 'GTVp', 'GTVn', 'CTVp', 'CTVn', 'PTVp', 'PTVn', 'phantom' or 'undefined'.
roiType_labels = {'rings','tumour','tumourAndEdema'}; % User-defined name of radiomics experiments to run. Each entry in the cell must correspond to one of the $nameExperiment$ in the "roiNames_$nameExperiment$.csv"files /WORKSPACE/CSV. Thus, this cell could have multiple entries.
% THIS OPTION BELOW IS DISABLED. AS OF MARCH 2018, THE USER MUST PROVIDE THE .csv FILES. THESE FILES MUST BE THOROUGHLY CHECKED BY THE USER BEFORE RUNNING THIS SCRIPT. THIS OPTION WILL SOON BE REPLACED BY AN OPTION ALLOWING TO MANUALLY VERIFY EACH ROI.
%manualROIchoice = false; % OPTION: % If set to true, the user will be prompted to choose which ROI(s) to use for a given imaging scan (and verify the chosen ROI definition on the imagin scans), so no need to create a roiNames.mat file --> this option is convenient for small datasets. If set to false, the algorithm will use the roiNames.mat file present in WORKSPACE --> this option is convenienent for large datasets, as manual reading and choice for all imaging scans in this script will be long.
segmentationMethod = 'InPolygon'; % Name describing the overall segmentation method, as defined by the Radiomics Ontology (RO) effort (see the RO instructions for more details). For the moment, must be set to either one of these character arrays: 'Poly2Mask', 'InPolygon', 'RegionGrowing', 'DeepLearning', 'CMeansFuzzy' or 'AtlasModel'.  

% POST-ACQUISITION IMAGE PROCESSING OPTIONS
partialVolumeEffectCorrection_MR = ''; % Value specifying if partial volume effect corrections were applied to MRI data prior to radiomics analysis. Insert 'Y' if it was performed, or '' if it was not. See the instructions of the Radiomics Ontology effort for more details. In future release of the code, this process will be integrated into script_ReadData.
partialVolumeEffectCorrection_CT = ''; % Value specifying if partial volume effect corrections were applied to CT data prior to radiomics analysis. Insert 'Y' if it was performed, or '' if it was not. See the instructions of the Radiomics Ontology effort for more details. In future release of the code, this process will be integrated into script_ReadData.
partialVolumeEffectCorrection_PET = ''; % Value specifying if partial volume effect corrections were applied to PET data prior to radiomics analysis. Insert 'Y' if it was performed, or '' if it was not. See the instructions of the Radiomics Ontology effort for more details. In future release of the code, this process will be integrated into script_ReadData.
noiseReduction_MR = ''; % Value specifying if noise reduction was applied to MRI data prior to radiomics analysis. Insert 'Y' if it was performed, or '' if it was not. See the instructions of the Radiomics Ontology effort for more details. In future release of the code, this process will be integrated into script_ReadData.
noiseReduction_CT = ''; % Value specifying if noise reduction was applied to CT data prior to radiomics analysis. Insert 'Y' if it was performed, or '' if it was not. See the instructions of the Radiomics Ontology effort for more details. In future release of the code, this process will be integrated into script_ReadData.
noiseReduction_PET = ''; % Value specifying if noise reduction was applied to PET data prior to radiomics analysis. Insert 'Y' if it was performed, or '' if it was not. See the instructions of the Radiomics Ontology effort for more details. In future release of the code, this process will be integrated into script_ReadData.
imageNonUniformityCorrection_MR = ''; % Value specifying if image nonuniformity corrections were applied to MRI data prior to radiomics analysis. Insert 'Y' if it was performed, or '' if it was not. See the instructions of the Radiomics Ontology effort for more details. In future release of the code, this process will be integrated into script_ReadData.
imageNonUniformityCorrection_CT = ''; % Value specifying if image nonuniformity corrections were applied to CT data prior to radiomics analysis. Insert 'Y' if it was performed, or '' if it was not. See the instructions of the Radiomics Ontology effort for more details. In future release of the code, this process will be integrated into script_ReadData.
imageNonUniformityCorrection_PET = ''; % Value specifying if image nonuniformity corrections were applied to PET data prior to radiomics analysis. Insert 'Y' if it was performed, or '' if it was not. See the instructions of the Radiomics Ontology effort for more details. In future release of the code, this process will be integrated into script_ReadData.

% PARALLEL OPTIONS
nBatch = 4; % OPTION: To compute radiomic features using parallelization. Put this number to the total number of cores that you have on your machine, but beware of RAM usage.
matlabPATH = 'C:\"Program Files"\MATLAB\R2016b\bin\matlab.exe'; % OPTION: IMPORTANT --> Full path to the matlab executable on the system. --> Ex: '/home/martin/Programs/MATLAB/R2017b/bin/matlab'. Here, a symbolic link to the full MATLAB path has previously been created on Martin Vallieres' computer. 
global codePATH
codePATH = 'C:\wrk\radiomics\Code';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    RADIOMIC FEATURE COMPUTATION CODE                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pathDATA = fullfile(pathWORK,'DATA'); pathCSV = fullfile(pathWORK,'CSV');
mkdir('FEATURES','RAW'), pathFEATURES = fullfile(pathWORK,'FEATURES','RAW'); 
pathTABLES = fullfile(pathWORK,'FEATURES'); cd(pathWORK)

tStart = tic;
fprintf('\n\n************************* RADIOMICS FEATURE COMPUTATION *************************')

% DISABLED AS OF MARCH 2018
% % 1. ROI CHOICE (optional, see above) OR LOADING PREVIOUSLY CONSTRUCTED roiNames.mat
% if manualROIchoice
%    flagRadiomics = roiChoice(pathWORK,pathDATA); % User will be prompted to choose each ROI for each scan of each patient and to verify if the chosen ROI read by the program is correctly positioned onto the imaging volume. Warning: at the end of the function, a new roiNames.mat files will be saved in the WORKSPACE, and potentially already existing ones will be given another name.
% else
%     flagRadiomics = true;
% end
% cd(pathWORK), load('roiNames') % Variable roiNames is now in MATLAB workspace.


% 1. MANUAL CHECK OF ROIs (optional, but recommended)
% --> To come soon in future release.



% 2. INITIALIZATIONS (in future release of the code, some parameters will be integrated to radiomics options)
imParamMR.morph.method = 'ISO'; imParamCT.morph.method = 'ISO'; imParamPET.morph.method = 'ISO'; imParamFilter.morph.method = 'ISO';
imParamMR.morph.value = 0.5; imParamCT.morph.value = 0.5; imParamPET.morph.value = 0.5; imParamFilter.morph.value = 0.5;
imParamMR.glcm.symmetry = 'SYM'; imParamCT.glcm.symmetry = 'SYM'; imParamPET.glcm.symmetry = 'SYM'; imParamFilter.glcm.symmetry = 'SYM';
imParamMR.glcm.distanceNorm.method = 'Chebyshev'; imParamCT.glcm.distanceNorm.method = 'Chebyshev'; imParamPET.glcm.distanceNorm.method = 'Chebyshev'; imParamFilter.glcm.distanceNorm.method = 'Chebyshev';
imParamMR.glcm.distanceNorm.value = 1; imParamCT.glcm.distanceNorm.value = 1; imParamPET.glcm.distanceNorm.value = 1; imParamFilter.glcm.distanceNorm.value = 1;
imParamMR.glcm.distanceNorm.unit = ''; imParamCT.glcm.distanceNorm.unit = ''; imParamPET.glcm.distanceNorm.unit = ''; imParamFilter.glcm.distanceNorm.unit = '';
imParamMR.glcm.distanceWeightingFunction = 'Inverse'; imParamCT.glcm.distanceWeightingFunction = 'Inverse'; imParamPET.glcm.distanceWeightingFunction = 'Inverse'; imParamFilter.glcm.distanceWeightingFunction = 'Inverse';
imParamMR.glrlm.distanceWeightingFunction = 'Inverse'; imParamCT.glrlm.distanceWeightingFunction = 'Inverse'; imParamPET.glrlm.distanceWeightingFunction = 'Inverse'; imParamFilter.glrlm.distanceWeightingFunction = 'Inverse';
imParamMR.gldzm.distanceNorm.method = 'Chebyshev'; imParamCT.gldzm.distanceNorm.method = 'Chebyshev'; imParamPET.gldzm.distanceNorm.method = 'Chebyshev'; imParamFilter.gldzm.distanceNorm.method = 'Chebyshev';
imParamMR.gldzm.distanceNorm.value = 1; imParamCT.gldzm.distanceNorm.value = 1; imParamPET.gldzm.distanceNorm.value = 1; imParamFilter.gldzm.distanceNorm.value = 1;
imParamMR.gldzm.distanceNorm.unit = ''; imParamCT.gldzm.distanceNorm.unit = ''; imParamPET.gldzm.distanceNorm.unit = ''; imParamFilter.gldzm.distanceNorm.unit = '';
imParamMR.ngtdm.distanceNorm.method = 'Chebyshev'; imParamCT.ngtdm.distanceNorm.method = 'Chebyshev'; imParamPET.ngtdm.distanceNorm.method = 'Chebyshev'; imParamFilter.ngtdm.distanceNorm.method = 'Chebyshev';
imParamMR.ngtdm.distanceNorm.value = 1; imParamCT.ngtdm.distanceNorm.value = 1; imParamPET.ngtdm.distanceNorm.value = 1; imParamFilter.ngtdm.distanceNorm.value = 1;
imParamMR.ngtdm.distanceNorm.unit = ''; imParamCT.ngtdm.distanceNorm.unit = ''; imParamPET.ngtdm.distanceNorm.unit = ''; imParamFilter.ngtdm.distanceNorm.unit = '';
imParamMR.ngtdm.distanceWeightingFunction = 'Inverse'; imParamCT.ngtdm.distanceWeightingFunction = 'Inverse'; imParamPET.ngtdm.distanceWeightingFunction = 'Inverse'; imParamFilter.ngtdm.distanceWeightingFunction = 'Inverse';
imParamMR.ngldm.coarseness = 0; imParamCT.ngldm.coarseness = 0; imParamPET.ngldm.coarseness = 0; imParamFilter.ngldm.coarseness = 0;
imParamMR.ngldm.distanceNorm.method = 'Chebyshev'; imParamCT.ngldm.distanceNorm.method = 'Chebyshev'; imParamPET.ngldm.distanceNorm.method = 'Chebyshev'; imParamFilter.ngldm.distanceNorm.method = 'Chebyshev';
imParamMR.ngldm.distanceNorm.value = 1; imParamCT.ngldm.distanceNorm.value = 1; imParamPET.ngldm.distanceNorm.value = 1; imParamFilter.ngldm.distanceNorm.value = 1;
imParamMR.ngldm.distanceNorm.unit = ''; imParamCT.ngldm.distanceNorm.unit = ''; imParamPET.ngldm.distanceNorm.unit = ''; imParamFilter.ngldm.distanceNorm.unit = '';
imParamMR.intensity = intensity_MR; imParamMR.units = ''; imParamMR.type = 'MRscan'; imParamCT.intensity = intensity_CT; imParamCT.units = 'HU'; imParamCT.type = 'CTscan'; imParamPET.intensity = intensity_PET; imParamPET.units = 'SUV'; imParamPET.type = 'PTscan'; imParamFilter.intensity = intensity_filter; imParamFilter.units = ''; imParamFilter.type = 'filter';
imParams.MRscan.image = imParamMR; imParams.CTscan.image = imParamCT; imParams.PTscan.image = imParamPET;
imParamFilterMR = imParamFilter; imParamFilterCT = imParamFilter; imParamFilterPT = imParamFilter;
imParamFilterMR.interp = imParamMR.interp; imParamFilterCT.interp = imParamCT.interp; imParamFilterPT.interp = imParamPET.interp;
imParamFilterMR.reSeg = imParamMR.reSeg; imParamFilterCT.reSeg = imParamCT.reSeg; imParamFilterPT.reSeg = imParamPET.reSeg;
imParamFilterMR = orderfields(imParamFilterMR,imParamMR); imParamFilterCT = orderfields(imParamFilterCT,imParamCT); imParamFilterPT = orderfields(imParamFilterPT,imParamPET);
if computeWavelet
    imParams.MRscan.LLL_coif1 = imParamFilterMR; imParams.MRscan.LHL_coif1 = imParamFilterMR; imParams.MRscan.LLH_coif1 = imParamFilterMR; imParams.MRscan.LHH_coif1 = imParamFilterMR; imParams.MRscan.HLL_coif1 = imParamFilterMR; imParams.MRscan.HHL_coif1 = imParamFilterMR; imParams.MRscan.HLH_coif1 = imParamFilterMR; imParams.MRscan.HHH_coif1 = imParamFilterMR;
    imParams.CTscan.LLL_coif1 = imParamFilterCT; imParams.CTscan.LHL_coif1 = imParamFilterCT; imParams.CTscan.LLH_coif1 = imParamFilterCT; imParams.CTscan.LHH_coif1 = imParamFilterCT; imParams.CTscan.HLL_coif1 = imParamFilterCT; imParams.CTscan.HHL_coif1 = imParamFilterCT; imParams.CTscan.HLH_coif1 = imParamFilterCT; imParams.CTscan.HHH_coif1 = imParamFilterCT;
    imParams.PTscan.LLL_coif1 = imParamFilterPT; imParams.PTscan.LHL_coif1 = imParamFilterPT; imParams.PTscan.LLH_coif1 = imParamFilterPT; imParams.PTscan.LHH_coif1 = imParamFilterPT; imParams.PTscan.HLL_coif1 = imParamFilterPT; imParams.PTscan.HHL_coif1 = imParamFilterPT; imParams.PTscan.HLH_coif1 = imParamFilterPT; imParams.PTscan.HHH_coif1 = imParamFilterPT;
    imParams.MRscan.filter.ToCompute = {['wavelet_',waveletName]}; imParams.CTscan.filter.ToCompute = {['wavelet_',waveletName]}; imParams.PTscan.filter.ToCompute = {['wavelet_',waveletName]};
    imParams.MRscan.filter.discretisation = imParamFilter.discretisation; imParams.CTscan.filter.discretisation = imParamFilter.discretisation; imParams.PTscan.filter.discretisation = imParamFilter.discretisation;
    imParams.MRscan.filter.intensity = imParamFilter.intensity; imParams.CTscan.filter.intensity = imParamFilter.intensity; imParams.PTscan.filter.intensity = imParamFilter.intensity;
    imParams.MRscan.filter.units = imParamFilter.units; imParams.CTscan.filter.units = imParamFilter.units; imParams.PTscan.filter.units = imParamFilter.units;
    imParams.MRscan.filter.type = imParamFilter.type; imParams.CTscan.filter.type = imParamFilter.type; imParams.PTscan.filter.type = imParamFilter.type;
end
imParams.MRscan.postAcquisitionProcessing.partialVolumeEffectCorrection = partialVolumeEffectCorrection_MR; imParams.CTscan.postAcquisitionProcessing.partialVolumeEffectCorrection = partialVolumeEffectCorrection_CT; imParams.PTscan.postAcquisitionProcessing.partialVolumeEffectCorrection = partialVolumeEffectCorrection_PET;
imParams.MRscan.postAcquisitionProcessing.noiseReduction = noiseReduction_MR; imParams.CTscan.postAcquisitionProcessing.noiseReduction = noiseReduction_CT; imParams.PTscan.postAcquisitionProcessing.noiseReduction = noiseReduction_PET;
imParams.MRscan.postAcquisitionProcessing.imageNonUniformityCorrection = imageNonUniformityCorrection_MR; imParams.CTscan.postAcquisitionProcessing.imageNonUniformityCorrection = imageNonUniformityCorrection_CT; imParams.PTscan.postAcquisitionProcessing.imageNonUniformityCorrection = imageNonUniformityCorrection_PET;
imParams.MRscan.segmentationMethod = segmentationMethod; imParams.CTscan.segmentationMethod = segmentationMethod; imParams.PTscan.segmentationMethod = segmentationMethod;
imParams.MRscan.calculationRun.dateTime = timeStamp; imParams.CTscan.calculationRun.dateTime = timeStamp; imParams.PTscan.calculationRun.dateTime = timeStamp;
imParams.MRscan.calculationRun.software.label = softwareLabel; imParams.CTscan.calculationRun.software.label = softwareLabel; imParams.PTscan.calculationRun.software.label = softwareLabel;
imParams.MRscan.calculationRun.software.version = softwareVersion; imParams.CTscan.calculationRun.software.version = softwareVersion; imParams.PTscan.calculationRun.software.version = softwareVersion;
imParams.MRscan.calculationRun.software.programmingLanguage = programmingLanguage; imParams.CTscan.calculationRun.software.programmingLanguage = programmingLanguage; imParams.PTscan.calculationRun.software.programmingLanguage = programmingLanguage;
imParams.MRscan.calculationRun.software.institution = institution; imParams.CTscan.calculationRun.software.institution = institution; imParams.PTscan.calculationRun.software.institution = institution;



% 3. RADIOMICS COMPUTATIONS
tic, fprintf('\n--> COMPUTING RADIOMICS FEATURES WITH %u CORES ... ',nBatch)
computeRadiomics_batchAllPatients(pathDATA,pathCSV,pathFEATURES,imParams,roiTypes,roiType_labels,nBatch,matlabPATH)
fprintf('\n    --> ALL ROI TYPES DONE!\n'), toc

tic, fprintf('\n--> COMPUTING RADIOMICS TABLES WITH %u CORES ... ',nBatch)
computeRadiomics_batchAllTables(pathFEATURES,pathTABLES,roiType_labels,nBatch,matlabPATH)
fprintf('DONE!\n'), toc

time = toc(tStart);
fprintf('\n\nTOTAL TIME FOR RADIOMICS COMPUTATIONS: %f seconds\n',time)
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


% DISABLED FOR NOW: MORE CHECKS ABOUT THE IMAGING DATA RANGE OF ALL PATIENTS FIRST NEED TO BE DEVELOPED
% % ADC MAP PARAMETERS: To be done specifically for that type of maps. Later, we may add KTrans as well. 
% imParamADC.interp.scaleNonText = [3,3,3]; % OPTION: Resolution in mm for the computation of non-texture features. Here, given the large difference between commonly seen in-plane resolution (~ 1 mm) and slice thickness (~ 5 mm), we take some middle point.
% imParamADC.interp.scaleText = {[1,1,1],[2,2,2],[3,3,3],[4,4,4]}; % OPTION: Different resolutions in mm tested for the computation of texture features. Each cell entry is a [X,Y,Z] resolution in MATLAB world coordinate.
% imParamADC.interp.volInterp = 'linear'; % OPTION: Using cubic interpolation for imaging intensities. Conservative choice to prevent interpolation craziness.
% imParamADC.interp.glRound = []; % OPTION: No grey-level rounding for ADC maps.
% imParamADC.interp.roiInterp = 'linear'; % OPTION: Using cubic interpolation for ROI mask. Conservative choice to prevent interpolation craziness.
% imParamADC.interp.roiPV = 0.5; % OPTION: After interpolation, a value >=0.5 is assigned to a 1, and <0.5 to a 0.
% imParamADC.reSeg.range = [0,inf]; % OPTION: We are working with ADC maps, going from 0 to infinity. The minimal value is of special importance for FBS discretisation. If the lower bound of the re-segmentation range is not defined, the minimum value of the ROI will be used in FBS discretisation (not recommended).
% imParamADC.reSeg.outliers = 'Collewet'; % OPTION: Using Collewet normalization to remove outliers for MRI.
% imParamADC.discretisation.IH.type = 'FBN'; % OPTION: Using fixed bin number for intensity-histogram features.
% imParamADC.discretisation.IH.val = 64; % OPTION: Using 64 grey-levels for intensity-histogram features.
% imParamADC.discretisation.IVH.type = 'FBS'; % OPTION: Using fixed bin size for intensity-volume histogram features with PET. imParam.reSeg.range must be defined, as userSetMinVal will be assigned to imParam.reSeg.range(1).
% imParamADC.discretisation.IVH.val = 1; % OPTION: Using 1 um^2/s bin width for intensity-histogram features.
% imParamADC.discretisation.texture.type = {'FBN','FBNequal','FBS','FBSequal'}; % OPTION: Using four different types of quantization algorithms. If "FBS" or "FBSequal" is used, imParam.reSeg.range must be defined, as userSetMinVal will be assigned to imParam.reSeg.range(1).
% imParamADC.discretisation.texture.val = {[8,16,32,64],[8,16,32,64],[12.5,25,50,100],[12.5,25,50,100]}; % OPTION: Gray-levels to test for each algorithm on the above line (definition depends on the algorithm). The total number must be the same in each cell entry. For the FBS algorithm, we assume values are in units of um^2/s.
% imParamADC.type = 'ADCscan';
% imParamADC.intensity = 'definite'; % OPTION: Definite intensity units (SUV).
% imParams.ADCscan.image = imParamADC; 
