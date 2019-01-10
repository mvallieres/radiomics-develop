function [ID,scan,type] = findIDscanType(name)
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


% ARGUMENT NAME MUST BE IN THE FOLLOWING FORMAT: ID_scan.modality.mat

indUnderscore = strfind(name,'_');
indDot = strfind(name,'.');

ID = name(1:(indUnderscore(1)-1));
scan = name((indUnderscore(1)+1):(indDot(1)-1));
type = name((indDot(1)+1):(indDot(2)-1));
end