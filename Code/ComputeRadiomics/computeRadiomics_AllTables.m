function computeRadiomics_AllTables(pathFEATURES,pathTABLES,tableTags)
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
nTables = size(tableTags,1);

fprintf('\n')
for t = 1:nTables
    scan= tableTags{t,1};
    roiType = tableTags{t,2};
    imSpace = tableTags{t,3};
    nameTable = ['radiomics__',scan,'(',roiType,')__',imSpace];
    wildcard = ['*',scan,'(',roiType,')*']; % Wildcard used to look only in the parent folder (pathFEATURES), no need to recursively look into sub-folders using '**/'.
    tic, fprintf('\n-->  Computing Radiomics table:"%s" ... ',nameTable);
    
    % Create radiomics table
    radiomicsTable = createRadiomicsTable(getFilePaths(pathFEATURES,wildcard),imSpace);
    radiomicsTable.Properties.Description = nameTable;
    
    % Save radiomics table
    save(fullfile(pathTABLES,nameTable),'radiomicsTable');
    
    % Create CSV table and Definitions
    writeRadiomicsCSV(fullfile(pathTABLES,nameTable))
    
    fprintf('DONE\n'), toc
end

cd(startpath)
end