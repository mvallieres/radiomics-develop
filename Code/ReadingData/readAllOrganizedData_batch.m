function readAllOrganizedData_batch(pathRead,pathSave,nBatch,pathMATLAB)
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
time = 60; % Number of seconds to wait before checking if parallel computations are done

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
    [~,~] = system(['rm batch',num2str(b),'_script.m']); % TO REMOVE AFTER DEBUGGING
    [~,~] = system(['rm batch',num2str(b),'_script.log']); % TO REMOVE AFTER DEBUGGING
    [~,~] = system(['rm batch',num2str(b),'_end']); % TO REMOVE AFTER DEBUGGING
    nameScript = ['batch',num2str(b),'_script.m'];
    fid = fopen(nameScript,'w');
    fprintf(fid,'load(''workspace'')\n');
    fprintf(fid,['readAllOrganizedData(pathRead,pathSave,namePatients(indPatient{',num2str(b),'}))\n']);
    fprintf(fid,['system(''touch batch',num2str(b),'_end'');\n']);
    fprintf(fid,'clear all');
    fclose(fid);
    system([pathMATLAB,' -nodisplay -nodesktop -nosplash -singleCompThread < ',nameScript,' >& ',nameScript(1:end-1),'log &']);
end

% WAITING LOOP
waitBatch(pathBatch,time,nBatch)

cd(startpath)
end
