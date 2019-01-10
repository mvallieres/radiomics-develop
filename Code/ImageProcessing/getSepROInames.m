function [nameROIout,vectPlusMinus] = getSepROInames(nameROIin,delimiters)
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

% WORKS ONLY FOR DELIMITERS "+" and "-"


% FINDING "+" and "-"
indPlus = strfind(nameROIin,delimiters{1}); vectPlus = ones(1,numel(indPlus));
indMinus = strfind(nameROIin,delimiters{2}); vectMinus = ones(1,numel(indMinus))*-1;
[ind,seq] = sort([indPlus,indMinus]);
vectPlusMinus = [vectPlus,vectMinus]; vectPlusMinus = vectPlusMinus(seq); nDelim = numel(vectPlusMinus);


% MAKING SURE "+" and "-" ARE NOT INSIDE A ROIname
indStart = strfind(nameROIin,'{'); nROI = numel(indStart);
indStop = strfind(nameROIin,'}');
indKeep = true(1,nDelim);
for d = 1:nDelim
    for r = 1:nROI
        if (indStop(r) - ind(d)) > 0 && (ind(d) - indStart(r)) > 0 % Thus inside a ROI name
            indKeep(d) = false;
            break
        end
    end
end
ind = ind(indKeep);
vectPlusMinus = vectPlusMinus(indKeep);


% PARSING ROI NAMES
if isempty(ind)
    nameROIout = {nameROIin(2:end-1)}; % Excluding the "{" and "}" at the start and end of the ROIname
else
    nInd = numel(ind);
    nameROIout = cell(1,nInd+1);
    nameROIout{1} = nameROIin(2:(ind(1)-2)); % Excluding the "{" and "}" at the start and end of the ROIname
    for i = 2:nInd
        nameROIout{i} = nameROIin((ind(i-1)+2):(ind(i)-2)); % Excluding the "{" and "}" at the start and end of the ROIname
    end
    nameROIout{nInd+1} = nameROIin((ind(end)+2):end-1); % Excluding the "{" and "}" at the start and end of the ROIname
end

end