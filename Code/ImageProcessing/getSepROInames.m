function nameROIout = getSepROInames(nameROIin,delimiter)
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


ind = strfind(nameROIin,delimiter);
if isempty(ind)
    nameROIout = {nameROIin};
else
    nInd = numel(ind);
    nameROIout = cell(1,nInd+1);
    nameROIout{1} = nameROIin(1:(ind(1)-1));
    for i = 2:nInd
        nameROIout{i} = nameROIin((ind(i-1)+1):(ind(i)-1));
    end
    nameROIout{nInd+1} = nameROIin((ind(end)+1):end);
end

end