function printStructure(structure,nameCSV)
% -------------------------------------------------------------------------
% AUTHOR(S): 
% - Martin Vallieres <mart.vallieres@gmail.com>
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

% THIS FUNCTION WORKS ONLY A ONE LEVEL STRUCTURE WITH NUMERICAL VALUES


% CREATE THE VECTOR OF VALUES
nameFields = fieldnames(structure); nFields = numel(nameFields);
vector = zeros(nFields,1);
for f = 1:nFields
    value = structure.(nameFields{f});
    if ~isnumeric(value) || isempty(value)
        vector(f) = NaN;
    else
        vector(f) = value;
    end
end


% CREATE TABLE
VALUE = vector;
tableValue = table(VALUE,'RowNames',nameFields);
tableValue.Properties.DimensionNames{1} = 'FEATURE';

% WRITE TABLE
writetable(tableValue,[nameCSV,'.csv'],'WriteRowNames',true);

end