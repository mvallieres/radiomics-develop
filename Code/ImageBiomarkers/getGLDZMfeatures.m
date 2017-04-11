function gldzm = getGLDZMfeatures(vol)
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

% - vol: 3D volume, isotropically resampled, quantized (e.g. Ng = 32, levels = [1, ..., Ng]), with NaNs outside the region of interest


% GET THE GLDZM MATRIX
levels = 1:max(vol(~isnan(vol(:)))); % Correct definition, without any assumption
[GLDZM] = getGLDZMmatrix(vol,levels); Ns = sum(GLDZM(:));
GLDZM = GLDZM./sum(GLDZM(:)); % Normalization of GLDZM
sz = size(GLDZM); % Size of GLSZM
cVect = 1:sz(2); rVect = 1:sz(1);% Row and column vectors
[cMat,rMat] = meshgrid(cVect,rVect); % Column and row indicators for each entry of the GLDZM
pg = sum(GLDZM,2)'; % Gray-Level Vector
pd = sum(GLDZM); % Distance Zone Vector



% COMPUTING TEXTURES

% Small distance emphasis
gldzm.Fdzm_sde = pd*(cVect.^(-2))';

% Large distance emphasis
gldzm.Fdzm_lde = pd*(cVect.^2)';

% Low grey level zone emphasis
gldzm.Fdzm_lgze = pg*(rVect.^(-2))';

% High grey level zone emphasis
gldzm.Fdzm_hgze = pg*(rVect.^2)';

% Small distance low grey level emphasis
gldzm.Fdzm_sdlge = sum(sum(GLDZM.*(rMat.^(-2)).*(cMat.^(-2))));

% Small distance high grey level emphasis
gldzm.Fdzm_sdhge = sum(sum(GLDZM.*(rMat.^2).*(cMat.^(-2))));

% Large distance low grey levels emphasis
gldzm.Fdzm_ldlge = sum(sum(GLDZM.*(rMat.^(-2)).*(cMat.^2)));

% Large distance high grey level emphasis
gldzm.Fdzm_ldhge = sum(sum(GLDZM.*(rMat.^2).*(cMat.^2)));

% Gray level non-uniformity
gldzm.Fdzm_glnu = sum(pg.^2) * Ns;

% Gray level non-uniformity normalised
gldzm.Fdzm_glnu_norm = sum(pg.^2);

% Zone distance non-uniformity
gldzm.Fdzm_zdnu = sum(pd.^2) * Ns;

% Zone distance non-uniformity normalised
gldzm.Fdzm_zdnu_norm = sum(pd.^2);

% Zone percentage
gldzm.Fdzm_z_perc = Ns/sum(~isnan(vol(:))); % Must change the original definition here.

% Grey level variance
temp = rMat .* GLDZM;
u = sum(temp(:));
temp = ((rMat - u).^2) .* GLDZM;
gldzm.Fdzm_gl_var = sum(temp(:));

% Zone distance variance
temp = cMat .* GLDZM;
u = sum(temp(:));
temp = ((cMat - u).^2) .* GLDZM;
gldzm.Fdzm_zd_var = sum(temp(:));

% Zone distance entropy
valPos = GLDZM(find(GLDZM(:)));
temp = valPos .* log2(valPos);
gldzm.Fdzm_zd_entr = -sum(temp(:));

end