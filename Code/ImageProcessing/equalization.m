function volEqual_RE = equalization(vol_RE)
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

% THIS IS A PURE "WHAT IS CONTAINED WITHIN THE ROI" EQUALIZATION. THIS IS
% NOT INFLUENCED BY THE "userSetMinVal" USED FOR FBS DISCRESTISATION.

Ng = 64; % This is the default we will use. It means that when using 'FBS', nQ should be chosen wisely such that the total number of grey levels does not exceed 64, for all patients (recommended). This choice was amde by considering that the best equalization performance for "histeq.m" is obtained with low Ng. WARNING: The effective number of grey levels coming out of "histeq.m" may be lower than Ng.


% CONSERVE THE INDICES OF THE ROI
Xgl = vol_RE(:);
indROI = find(~isnan(Xgl));
Xgl = Xgl(~isnan(Xgl));

% ADJUST RANGE BETWEEN 0 and 1
minVal = min(Xgl);
maxVal = max(Xgl);
Xgl01 = (Xgl - minVal)/(maxVal - minVal);

% EQUALIZATION
Xgl_equal = histeq(Xgl01,Ng);

% RE-ADJUST TO CORRECT RANGE
Xgl_equal = (Xgl_equal - min(Xgl_equal))/(max(Xgl_equal) - min(Xgl_equal));
Xgl_equal = Xgl_equal*(maxVal - minVal);
Xgl_equal = Xgl_equal + minVal;

% RECONSTRUCT THE VOLUME WITH EQUALIZED VALUES
volEqual_RE = vol_RE;
volEqual_RE(indROI) = Xgl_equal;

end
