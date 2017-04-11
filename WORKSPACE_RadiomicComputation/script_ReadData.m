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

clc,clear,fprintf('\n'), warning off
help script_ReadData, pathWORK = pwd;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          PARAMETER OPTIONS                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PARALLEL OPTIONS
nBatch_Read = 2; % For the initial reading of the data using parallelization. Beware: RAM usage limitations. In doubt, just put 1.
matlabPATH = 'matlab'; % IMPORTANT: Full path to the matlab executable on the system. --> Ex: '/home/martin/Programs/MATLAB/R2016a/bin/matlab'. Here, a symbolic link to the full MATLAB path has previously been created on Martin Vallieres' computer. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      READING DICOM DATA CODE                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pathDICOM = fullfile(pathWORK,'DICOM'); 
mkdir('DATA'), pathDATA = fullfile(pathWORK,'DATA'); 
cd(pathWORK)

tStart = tic;
fprintf('\n\n************************* READING ALL DICOM DATA *************************')

% 1 PRE-ANONYMIZING THE DATA USING "anonymize_dicom.py"
tic, fprintf('\n--> PRE-ANONYMIZING THE DATA ... ')
preAnonymize(pathDICOM)
fprintf('DONE!\n'), toc

% 2 READING DATA (organization of patient and scans folders is crucial here)
tic, fprintf('\n--> READING DATA USING %u CORES ... ',nBatch_Read)
readAllDICOM(pathDICOM,pathDATA,nBatch_Read,matlabPATH,'folder')
fprintf('DONE!\n'), toc

time = toc(tStart);
fprintf('\n\nTOTAL TIME FOR READING ALL DICOM DATA: %f seconds\n',time)
fprintf('-------------------------------------------------------------------------------------')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
