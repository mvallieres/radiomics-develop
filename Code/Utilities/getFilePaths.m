function filePaths = getFilePaths(pathToParentFolder,wildcard)
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

% IMPORTANT
% - pathFiles: Full path to where the files are located 
%
% - wildcard: (optional). String specifying which type of files to locate in the parent folder.
%    Ex: '*.dcm*'. If not present --> same as wildcard '*' (look for all files, but only in the parent folder, i.e. non-recursively).
%    --> Hint: To recursively look for all .dcm files in all sub-folders of the parent folder, use
%    fullfile('**','*.dcm') as the wildcard (i.e. '**/*.dcm' on linux)
%
% --> We only want names of files. We thus remove directories with potentially the
%         same wildcard name


if nargin < 2
    wildcard = '*';
end

% Getting the list of all files
list = dir(fullfile(pathToParentFolder,wildcard));

% Removing directories 
listCell = struct2cell(list);
keepFlag = ~cell2mat(listCell (5,:));  % Fifth row tells us if each entry is a directory or not.
list = list(keepFlag); 

% Getting filePaths
filePaths = fullfile({list.folder},{list.name})';

end