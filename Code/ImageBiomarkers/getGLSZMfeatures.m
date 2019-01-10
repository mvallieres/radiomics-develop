function glszm = getGLSZMfeatures(vol)
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


% GET THE GLSZM MATRIX
levels = 1:max(vol(~isnan(vol(:)))); % Correct definition, without any assumption
[GLSZM] = getGLSZMmatrix(vol,levels); Ns = sum(GLSZM(:));
GLSZM = GLSZM./sum(GLSZM(:)); % Normalization of GLSZM
sz = size(GLSZM); % Size of GLSZM
cVect = 1:sz(2); rVect = 1:sz(1);% Row and column vectors
[cMat,rMat] = meshgrid(cVect,rVect); % Column and row indicators for each entry of the GLSZM
pg = sum(GLSZM,2)'; % Gray-Level Vector
pz = sum(GLSZM); % Zone Size Vector



% COMPUTING TEXTURES

% Small zone emphasis
glszm.Fszm_sze = pz*(cVect.^(-2))';

% Large zone emphasis
glszm.Fszm_lze = pz*(cVect.^2)';

% Low grey level zone emphasis
glszm.Fszm_lgze = pg*(rVect.^(-2))';

% High grey level zone emphasis
glszm.Fszm_hgze = pg*(rVect.^2)';

% Small zone low grey level emphasis
glszm.Fszm_szlge = sum(sum(GLSZM.*(rMat.^(-2)).*(cMat.^(-2))));

% Small zone high grey level emphasis
glszm.Fszm_szhge = sum(sum(GLSZM.*(rMat.^2).*(cMat.^(-2))));

% Large zone low grey levels emphasis
glszm.Fszm_lzlge = sum(sum(GLSZM.*(rMat.^(-2)).*(cMat.^2)));

% Large zone high grey level emphasis
glszm.Fszm_lzhge = sum(sum(GLSZM.*(rMat.^2).*(cMat.^2)));

% Gray level non-uniformity
glszm.Fszm_glnu = sum(pg.^2) * Ns;

% Gray level non-uniformity normalised
glszm.Fszm_glnu_norm = sum(pg.^2);

% Zone size non-uniformity
glszm.Fszm_zsnu = sum(pz.^2) * Ns;

% Zone size non-uniformity normalised
glszm.Fszm_zsnu_norm = sum(pz.^2);

% Zone percentage
glszm.Fszm_z_perc = sum(pg)/(pz*cVect');

% Grey level variance
temp = rMat .* GLSZM;
u = sum(temp(:));
temp = ((rMat - u).^2) .* GLSZM;
glszm.Fszm_gl_var = sum(temp(:));

% Zone size variance
temp = cMat .* GLSZM;
u = sum(temp(:));
temp = ((cMat - u).^2) .* GLSZM;
glszm.Fszm_zs_var = sum(temp(:));

% Zone size entropy
valPos = GLSZM(find(GLSZM(:)));
temp = valPos .* log2(valPos);
glszm.Fszm_zs_entr = -sum(temp(:));

end