function ngldm = getNGLDMfeatures(vol)
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

% - vol: 3D volume, isotropically resampled, quantized, (e.g. Ng = 32, levels = [1, ..., Ng]) with NaNs outside the region of interest


% GET THE NGLDM MATRIX
levels = 1:max(vol(~isnan(vol(:)))); % Correct definition, without any assumption
[NGLDM] = getNGLDMmatrix(vol,levels); Ns = sum(NGLDM(:));
NGLDM = NGLDM./sum(NGLDM(:)); % Normalization of NGLDM
sz = size(NGLDM); % Size of NGLDM
cVect = 1:sz(2); rVect = 1:sz(1);% Row and column vectors
[cMat,rMat] = meshgrid(cVect,rVect); % Column and row indicators for each entry of the NGLDM
pg = sum(NGLDM,2)'; % Gray-Level Vector
pd = sum(NGLDM); % Dependence Count Vector



% COMPUTING TEXTURES

% Low dependence emphasis
ngldm.Fngl_lde = pd*(cVect.^(-2))';

% High dependence emphasis
ngldm.Fngl_hde = pd*(cVect.^2)';

% Low grey level count emphasis
ngldm.Fngl_lgce = pg*(rVect.^(-2))';

% High grey level count emphasis
ngldm.Fngl_hgce = pg*(rVect.^2)';

% Low dependence low grey level emphasis
ngldm.Fngl_ldlge = sum(sum(NGLDM.*(rMat.^(-2)).*(cMat.^(-2))));

% Low dependence high grey level emphasis
ngldm.Fngl_ldhge = sum(sum(NGLDM.*(rMat.^2).*(cMat.^(-2))));

% High dependence low grey levels emphasis
ngldm.Fngl_hdlge = sum(sum(NGLDM.*(rMat.^(-2)).*(cMat.^2)));

% High dependence high grey level emphasis
ngldm.Fngl_hdhge = sum(sum(NGLDM.*(rMat.^2).*(cMat.^2)));

% Gray level non-uniformity
ngldm.Fngl_glnu = sum(pg.^2) * Ns;

% Gray level non-uniformity normalised
ngldm.Fngl_glnu_norm = sum(pg.^2);

% Dependence count non-uniformity
ngldm.Fngl_dcnu = sum(pd.^2) * Ns;

% Dependence count non-uniformity normalised
ngldm.Fngl_dcnu_norm = sum(pd.^2);

% Dependence count percentage
% Omitted, always evaluates to 1.

% Grey level variance
temp = rMat .* NGLDM;
u = sum(temp(:));
temp = ((rMat - u).^2) .* NGLDM;
ngldm.Fngl_gl_var = sum(temp(:));

% Dependence count variance
temp = cMat .* NGLDM;
u = sum(temp(:));
temp = ((cMat - u).^2) .* NGLDM;
ngldm.Fngl_dc_var = sum(temp(:));

% Dependence count entropy
valPos = NGLDM(find(NGLDM(:)));
temp = valPos .* log2(valPos);
ngldm.Fngl_dc_entr = -sum(temp(:));

% Dependence count energy
temp = NGLDM.^2;
ngldm.Fngl_dc_energy = sum(temp(:));

end