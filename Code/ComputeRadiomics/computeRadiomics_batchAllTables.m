function computeRadiomics_batchAllTables(pathFEATURES,pathTABLES,roiType_labels,nBatchInit,matlabPATH,codePATH)
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
time = 10; % Time spent in seconds between checks by the master process to verify the end of parallel computations.



% GETTING COMBINATIONS OF scan,  roiType and imageSpaces
nROItypes = numel(roiType_labels);
tableTags = {};
for r = 1:nROItypes
    label = roiType_labels{r};
    
    % Get all scan names present for that given roiType label
    filePaths = getFilePaths(pathFEATURES,['*',label,'*']); nFiles = numel(filePaths);
    scans = cell(nFiles,1);
    for f = 1:nFiles
        [~,radFileName,~] = fileparts(filePaths{f}); 
        scans{f} = getScanName_fromRadName(radFileName); 
    end
    scans = unique(scans); nScans = numel(scans);
    
    for s = 1:nScans
        scan = scans{s};
        
        % Get all scan names present for that given roiType label and scans
        filePaths = getFilePaths(pathFEATURES,['*',scan,'(',label,')*']); nFiles = numel(filePaths);

        % Finding the images spaces for a test file (assuming that all
        % files for a given scan and roiType label have the same image
        % spaces
        radiomics = load(filePaths{ceil(nFiles*rand(1,1))}); radiomics = struct2cell(radiomics); radiomics = radiomics{1};
        imSpaces = fieldnames(radiomics); imSpaces(end) = []; nImSpaces = numel(imSpaces); % The last field of radiomics is "imParam" --> Not an image space.

        % Constructing the tableTags variable
        for i = 1:nImSpaces
            imSpace = imSpaces{i};
            tableTags = [tableTags;{scan,label,imSpace}];
        end
    
    end
end
nTables = size(tableTags,1);



% INITIALIZATION
cd(pathTABLES), nameBatchLog = 'batchLog_tables';
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
if nTables < nBatch
    nBatch = nTables;
end
[tables] = batchPatients(nTables,nBatch);
cd(pathBatch), save('workspace','pathFEATURES','pathTABLES','tables','tableTags'), pause(2);
for i = 1:nBatch
    nameScript = ['batch',num2str(i),'_script.m'];
    fid = fopen(nameScript,'w');
    fprintf(fid,'load(''workspace'')\n');
    addpathline = ['addpath(genpath(''', codePATH,'''));'];
    fprintf(fid,'%s\n',addpathline); %Windows paths contain \ which are interpreted as escape characters and cannot be used in the fprintf format string
    fprintf(fid,['computeRadiomics_AllTables(pathFEATURES,pathTABLES,tableTags(tables{',num2str(i),'},:));\n']);
    if ispc
        fprintf(fid,['system(''type nul > batch', num2str(i),'_end'');\n']);
    else
    	fprintf(fid,['system(''touch batch',num2str(i),'_end'');\n']);
    end
    fprintf(fid,'clear all');
    fclose(fid);
    if ispc
        system(['start /B ',matlabPATH,' -nodisplay -nodesktop -nosplash -singleCompThread -r "diary ',nameScript(1:end-1),'log;',nameScript(1:end-2),';diary off;exit" ']);
        %system(['start /B ',matlabPATH,' -nodisplay -nodesktop -nosplash -singleCompThread < ',nameScript,' > ',nameScript(1:end-1),'log 2>&1']);
    else
    	system([matlabPATH,' -nodisplay -nodesktop -nosplash -singleCompThread < ',nameScript,' >& ',nameScript(1:end-1),'log &']);
    end
end


% WAITING LOOP
waitBatch(pathBatch,time,nBatch)
%delete('workspace.mat')

cd(startpath)
end
