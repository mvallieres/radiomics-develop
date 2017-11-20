function [volQuant_RE,wd] = discretisation(vol_RE,type,nQ,userSetMinVal,ivh)
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

% --> Use useSetValue as the minimum of range re-segmentation. This is the
% only way to create comparable textures using FBS discretisation. For FBN
% discretisation, this value has no importance as an argument to the
% function and will not be used.
% --> Last argument is optional and MUST BE SET TO 'ivh' FOR IVH FEATURES!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IMPORTANT --> FOR 'FBS' TYPE, IT IS ASSUMED THAT RE-SEGMENTATION WITH   %
%               PROPER RANGE WAS ALREADY PERFORMED.                       % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% -------------------------------------------------------------------------
% PARSING ARGUMENTS
volQuant_RE = vol_RE;
if isempty(nQ)
    return
end
if nargin < 3
    error('The total number of arguments must be at least be 3')
else
    if ~isnumeric(nQ)
        error('Third argument must be a numerical number')
    end
    if ~strcmp(type,'FBS') && ~strcmp(type,'FBN') && ~strcmp(type,'FBSequal') && ~strcmp(type,'FBNequal')
        error('The second argument must either be "FBS" or "FBN" or "FBSequal" or "FBNequal"')
    end
end
% -------------------------------------------------------------------------


% DISCRETISATION
if strcmp(type,'FBS') || strcmp(type,'FBSequal')
    if ~isempty(userSetMinVal)
        minVal = userSetMinVal;
    else
        minVal = min(volQuant_RE(:));
    end
else
    minVal = min(volQuant_RE(:));
end
maxVal = max(volQuant_RE(:));
switch type
    case 'FBS'
        wb = nQ; wd = wb;
        volQuant_RE = ceil((volQuant_RE - minVal)/nQ);
        if min(volQuant_RE(:)) == 0
            volQuant_RE(volQuant_RE == min(volQuant_RE(:))) = 1;
        end
    case 'FBN'
        wb = (maxVal - minVal)/nQ; wd = 1;
        volQuant_RE = ceil(nQ * ((volQuant_RE - minVal)/(maxVal - minVal)));
        volQuant_RE(volQuant_RE == min(volQuant_RE(:))) = 1;
    case 'FBSequal'
        wb = nQ; wd = wb;
        volQuant_RE = equalization(volQuant_RE);
        volQuant_RE = ceil((volQuant_RE - minVal)/nQ);
        if min(volQuant_RE(:)) == 0
            volQuant_RE(volQuant_RE == min(volQuant_RE(:))) = 1;
        end
    case 'FBNequal'
        wb = (maxVal - minVal)/nQ; wd = 1;
        volQuant_RE = equalization(volQuant_RE);
        volQuant_RE = ceil(nQ * ((volQuant_RE - minVal)/(maxVal - minVal)));
        volQuant_RE(volQuant_RE == min(volQuant_RE(:))) = 1;
end

if nargin == 5 && strcmp(ivh,'ivh') && (strcmp(type,'FBS') || strcmp(type,'FBSequal'))
    volQuant_RE = minVal + (volQuant_RE - 0.5)*wb;
end

end