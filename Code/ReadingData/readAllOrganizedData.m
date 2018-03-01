function readAllOrganizedData(pathRead,pathSave,namePatients)
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

startpath = pwd;

nPatients = numel(namePatients);

fprintf('\n');
tStart = tic;
for p = 1:nPatients
    namePatient = namePatients{p};
    fprintf('\n--> Reading "%s" ... ',namePatient);
    readOrganizedData(fullfile(pathRead,namePatient),pathSave,namePatient);
    fprintf('DONE');
end
time = toc(tStart);
fprintf('\n\n TOTAL TIME FOR READING %u PATIENTS: %f min',nPatients,time/60);

cd(startpath)
end
