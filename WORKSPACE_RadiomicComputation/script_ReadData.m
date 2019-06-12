% ******************************************************************************
% * DESCRIPTION:                                                               *
% * This script reads all organized DICOM data.                                *
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

clc,clear,fprintf('\n'), warning off, pathWORK = pwd;
if isempty(mfilename) % That means we are running the script from the command line using the following format (mandatory): matlab -nodisplay -nodesktop -nosplash -singleCompThread < script_ReadData.m >& script_ReadData.log &
    [scriptFileName,logFileName] = scriptCommandLine_FindNames(pathWORK);
else
    scriptFileName = [mfilename,'.m'];
end
help(scriptFileName)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          PARAMETER OPTIONS                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PARALLEL OPTIONS
nBatch_Read = 1; % For the initial reading of the data using parallelization. Beware: RAM usage limitations. In doubt, just put 1.
matlabPATH = 'matlab'; % IMPORTANT: Full path to the matlab executable on the system. --> Ex: '/home/martin/Programs/MATLAB/R2016a/bin/matlab'. Here, a symbolic link to the full MATLAB path has previously been created on Martin Vallieres' computer. 
global codePATH
codePATH = '~/GitHub/radiomics-develop/Code'; % To be adapted to your system. Full path to the "radiomics-develop" code on your system.
addpath(genpath(codePATH));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      READING DICOM DATA CODE                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pathDICOM = fullfile(pathWORK,'DICOM'); 
mkdir('DATA'), pathDATA = fullfile(pathWORK,'DATA'); 
cd(pathWORK)

tStart = tic;
fprintf('\n\n************************* READING ALL DICOM DATA *************************')

% 1. READING DATA (organization of patient and scan folders as defined in the instructions is mandatory)
tic, fprintf('\n--> READING DATA USING %u CORES ... ',nBatch_Read)
readAllOrganizedData_batch(pathDICOM,pathDATA,nBatch_Read,matlabPATH) % If data is not organized and only in DICOM format, it is still possible to use the function "readAllDICOM.m'.
fprintf('DONE!\n'), toc

time = toc(tStart);
fprintf('\n\nTOTAL TIME FOR READING ALL DATA: %f seconds\n',time)
fprintf('-------------------------------------------------------------------------------------')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
