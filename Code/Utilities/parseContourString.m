function [contourNumber,operations] = parseContourString(contourString)
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

if isnumeric(contourString)
    contourNumber = contourString;
    operations = [];
    return
end

indPlus = strfind(contourString,'+');
indMinus = strfind(contourString,'-');
indOperations = sort([indPlus,indMinus]); 

% Parsing operations
if isempty(indOperations)
    operations = [];
else
    nOp = numel(indOperations); operations = cell(1,nOp);
    for o = 1:nOp
        operations{o} = contourString(indOperations(o));
    end
end

% Parsing contour numbers
if isempty(indOperations)
    contourNumber = str2num(contourString);
else
    contourNumber = zeros(1,nOp+1);
    contourNumber(1) = str2num(contourString(1:indOperations(1)-1));
    for c = 2:nOp
        contourNumber(c) = str2num(contourString((indOperations(c-1)+1):(indOperations(c)-1)));
    end
    contourNumber(end) = str2num(contourString((indOperations(end)+1):end));
end

end
