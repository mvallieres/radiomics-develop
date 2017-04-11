function [structureArray] = appendStruct(structureArray,newStructure)
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


if isempty(structureArray) || isempty(fieldnames(structureArray))
    structureArray = newStructure;
    return
end

structLength = length(structureArray);
fields = fieldnames(structureArray(1));
nFields = length(fields);

for i = 1:nFields
    try
        structureArray(structLength + 1).(fields{i}) = newStructure.(fields{i});
    catch
        structureArray(structLength + 1).(fields{i}) = 'FIELD NOT PRESENT';
    end
end
