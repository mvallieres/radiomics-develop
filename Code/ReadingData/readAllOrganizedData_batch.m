function readAllOrganizedData_batch(pathRead,pathSave,nBatch,matlabPATH,codePATH)
% -------------------------------------------------------------------------
% AUTHOR(S): 
% - Martin Vallieres <mart.vallieres@gmail.com>
% -------------------------------------------------------------------------
% HISTORY:
% - Creation: March 2018
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

warning off
startpath = pwd;
time = 10; % Number of seconds to wait before checking if parallel computations are done

% INITIALIZATION
cd(pathRead), listPatients = dir;
listPatients = listPatients(~ismember({listPatients.name},{'.','..','.DS_Store','._.DS_Store'}));
indRemove = [];
for p = 1:numel(listPatients)
    indHyphen = strfind(listPatients(p).name,'-');
    if numel(indHyphen) ~= 2
        indRemove = [indRemove,p];
    end
end
listPatients(indRemove) = []; % We remove any folder that do not follow conventions: $cancer$-$institution$-000
nPatients = numel(listPatients); namePatients = cell(1,nPatients);
for p = 1:nPatients
    namePatients{p} = listPatients(p).name;
end


% STARTING BATCH COMPUTATIONS FOR READING EACH PATIENT FOLDER
cd(pathSave), mkdir('batchLog_ReadData'), cd('batchLog_ReadData'), pathBatch = pwd;
if nPatients < nBatch
    nBatch = nPatients;
end
[indPatient] = batchPatients(nPatients,nBatch);
cd(pathBatch), save('workspace','pathRead','pathSave','namePatients','indPatient')
for b = 1:nBatch
    nameScript = ['batch',num2str(b),'_script.m'];
    fid = fopen(nameScript,'w');
    fprintf(fid,'load(''workspace'')\n');
    addpathline = ['addpath(genpath(''', codePATH,'''));'];
    fprintf(fid,'%s\n',addpathline); %Windows paths contain \ which are interpreted as escape characters and cannot be used in the fprintf format string
    fprintf(fid,['readAllOrganizedData(pathRead,pathSave,namePatients(indPatient{',num2str(b),'}))\n']);
    if ispc
        fprintf(fid,['system(''type nul > batch', num2str(b),'_end'');\n']);
    else
    	fprintf(fid,['system(''touch batch',num2str(b),'_end'');\n']);
    end
    fprintf(fid,'clear all');
    fclose(fid);
    if ispc
        system(['start /B ',matlabPATH,' -nodisplay -nodesktop -nosplash -singleCompThread -r "diary ',nameScript(1:end-1),'log;',nameScript(1:end-2),';diary off;exit" ']);
    else
    	system([matlabPATH,' -nodisplay -nodesktop -nosplash -singleCompThread < ',nameScript,' >& ',nameScript(1:end-1),'log &']);
    end
end

% WAITING LOOP
waitBatch(pathBatch,time,nBatch)


% DISABLED FOR NOW: MORE CHECKS ABOUT THE IMAGING DATA RANGE OF ALL PATIENTS FIRST NEED TO BE DEVELOPED
% % CHECKING AND MODIFYING ADC DATA
% modifADCdata(pathSave)

cd(startpath)
end
