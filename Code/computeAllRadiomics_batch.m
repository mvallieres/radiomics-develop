function computeAllRadiomics_batch(pathRead,pathSave,roiNames,imParams,roiType,nBatch,matlabPATH)
% -------------------------------------------------------------------------
% AUTHOR(S): 
% - Martin Vallieres <mart.vallieres@gmail.com>
% -------------------------------------------------------------------------
% HISTORY:
% - Creation: April 2017
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


startpath = pwd;

cd(pathSave)
nameBatchLog = ['batchLog_',roiType];
if exist(nameBatchLog,'dir')
    info = dir([nameBatchLog,'*']);
    date = info.date; ind = strfind(date,' '); date(ind) = '_';
    newName = [nameBatchLog,'_',date];
    system(['mv ',nameBatchLog,' ',newName]);    
end
mkdir(nameBatchLog), cd(nameBatchLog), pathBatch = pwd;
time = 60; % Time spent in seconds between checks by the master process to verify the end of parallel computations.
patientNames = getPatientNames(roiNames);
nameROI = roiNames(:,4);
nameSet = roiNames(:,5);

% PRODUCE BATCH COMPUTATIONS
if numel(patientNames) < nBatch
    nBatch = numel(patientNames);
end
[patients] = batchPatients(numel(patientNames),nBatch);
cd(pathBatch), save('workspace','pathRead','pathSave','patients','patientNames','nameROI','nameSet','imParams','roiType'), pause(2);
for i = 1:nBatch
    nameScript = ['batch',num2str(i),'_script.m'];
    fid = fopen(nameScript,'w');
    fprintf(fid,'load(''workspace'')\n');
    fprintf(fid,['computeAllRadiomics(pathRead,pathSave,patientNames(patients{',num2str(i),'}),nameROI(patients{',num2str(i),'}),nameSet(patients{',num2str(i),'}),imParams,roiType)\n']);
    fprintf(fid,['system(''touch batch',num2str(i),'_end'');\n']);
    fprintf(fid,'clear all');
    fclose(fid);
    system([matlabPATH,' -nojvm -nodisplay -nodesktop -nosplash < ',nameScript,' >& ',nameScript(1:end-1),'log &']);
end

% WAITING LOOP
waitBatch(pathBatch,time,nBatch)
delete('workspace.mat')

cd(startpath)
end