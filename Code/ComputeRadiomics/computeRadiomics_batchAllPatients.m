function computeRadiomics_batchAllPatients(pathRead,pathCSV,pathSave,imParams,roiTypes,roiType_labels,nBatchInit,matlabPATH)
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

global codePATH

startpath = pwd;

time = 10; % Time spent in seconds between checks by the master process to verify the end of parallel computations.
nROItypes = numel(roiType_labels);

for r = 1:nROItypes
    roiType = roiTypes{r};
    roiType_label = roiType_labels{r};
    fprintf('\n    --> Computing features for the "%s" roi type ... ',roiType_label);
    
    % READING CSV EXPERIMENT TABLE
    cd(pathCSV), tableROI = readtable(['roiNames_',roiType_label,'.csv']);
    patientNames = getPatientNames([tableROI.PatientID,tableROI.ImagingScanName,tableROI.ImagingModality]);
    nameROI = tableROI.ROIname;
    if sum(contains(tableROI.Properties.VariableNames,'StructureSetName'))
        nameSet = tableROI.StructureSetName;
    else
        nameSet = cell(size(tableROI,1),1); % Creation of cell with empty entries.
    end
    
    
    % INITIALIZATION
    cd(pathSave), nameBatchLog = ['batchLog_',roiType_label];
    if exist(nameBatchLog,'dir')
        info = dir([nameBatchLog,'*']);
        date = info.date; ind = strfind(date,' '); date(ind) = '_';
        newName = [nameBatchLog,'_',date];
        if ispc
            system(['move ', nameBatchLog, ' ', newName]);
        else
            system(['mv ',nameBatchLog,' ',newName]); 
        end	
    end
    mkdir(nameBatchLog), cd(nameBatchLog), pathBatch = pwd;

    
    % PRODUCE BATCH COMPUTATIONS
    nBatch = nBatchInit;
    if numel(patientNames) < nBatch
        nBatch = numel(patientNames);
    end
    [patients] = batchPatients(numel(patientNames),nBatch);
    cd(pathBatch), save('workspace','pathRead','pathSave','patients','patientNames','nameROI','nameSet','imParams','roiType','roiType_label'), pause(2);
    for i = 1:nBatch
        nameScript = ['batch',num2str(i),'_script.m'];
        fid = fopen(nameScript,'w');
        fprintf(fid,'load(''workspace'')\n');
        addpathline = ['addpath(genpath(''', codePATH,'''));'];
        fprintf(fid,'%s\n',addpathline); %Windows paths contain \ which are interpreted as escape characters and cannot be used in the fprintf format string
        fprintf(fid,['computeRadiomics_AllPatients(pathRead,pathSave,patientNames(patients{',num2str(i),'}),nameROI(patients{',num2str(i),'}),nameSet(patients{',num2str(i),'}),imParams,roiType,roiType_label)\n']);
        if ispc
            fprintf(fid,['system(''type nul > batch', num2str(i),'_end'');\n']);
        else
            fprintf(fid,['system(''touch batch',num2str(i),'_end'');\n']);
        end
        fprintf(fid,'clear all');
        fclose(fid);
        if ispc
            system(['start /B ',matlabPATH,' -nodisplay -nodesktop -nosplash -singleCompThread -r "diary ',nameScript(1:end-1),'log;',nameScript(1:end-2),';diary off;exit" ']);
        else
            system([matlabPATH,' -nodisplay -nodesktop -nosplash -singleCompThread -r "diary ',nameScript(1:end-1),'log;',nameScript(1:end-2),';diary off;exit" & ']);
        end
    end

    
    % WAITING LOOP
    waitBatch(pathBatch,time,nBatch)
    delete('workspace.mat')

    fprintf('DONE');
end

cd(startpath)
end
