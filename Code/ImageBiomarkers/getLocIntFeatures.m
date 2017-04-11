function locInt = getLocIntFeatures(imgObj,roiObj,res)
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

% - imgObj: Continous image intentisity distribution, with no NaNs outside the ROI
% - roiObj: Mask defining the ROI
% - res: [a,b,c] vector specfying the resolution of the volume in mm.  % XYZ resolution (world), or JIK resolution (intrinsic matlab).


% Local grey level peak
locInt.Floc_peak_loc = getLocPeak(imgObj,roiObj,res);

% Global grey level peak
% locInt.Floc_peak_glob = getGlobPeak(imgObj,roiObj,res); % NEEDS TO BE VECTORIZED FOR FASTER CALCULATION!

end
