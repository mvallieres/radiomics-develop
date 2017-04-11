function glrlm = getGLRLMfeatures(vol)
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


% GET THE GLRLM MATRIX
levels = 1:max(vol(~isnan(vol(:)))); % Correct definition, without any assumption
[GLRLM] = getGLRLMmatrix(vol,levels); Ns = sum(GLRLM(:));
GLRLM = GLRLM./sum(GLRLM(:)); % Normalization of GLRLM
sz = size(GLRLM); % Size of GLRLM
cVect = 1:sz(2); rVect = 1:sz(1);% Row and column vectors
[cMat,rMat] = meshgrid(cVect,rVect); % Column and row indicators for each entry of the GLRLM
pg = sum(GLRLM,2)'; % Gray-Level Run-Number Vector
pr = sum(GLRLM); % Run-Length Run-Number Vector



% COMPUTING TEXTURES

% Short runs emphasis
glrlm.Frlm_sre = pr*(cVect.^(-2))';

% Long runs emphasis
glrlm.Frlm_lre = pr*(cVect.^2)';

% Low grey level run emphasis
glrlm.Frlm_lgre = pg*(rVect.^(-2))';

% High grey level run emphasis
glrlm.Frlm_hgre = pg*(rVect.^2)';

% Short run low grey level emphasis
glrlm.Frlm_srlge = sum(sum(GLRLM.*(rMat.^(-2)).*(cMat.^(-2))));

% Short run high grey level emphasis
glrlm.Frlm_srhge = sum(sum(GLRLM.*(rMat.^2).*(cMat.^(-2))));

% Long run low grey levels emphasis
glrlm.Frlm_lrlge = sum(sum(GLRLM.*(rMat.^(-2)).*(cMat.^2)));

% Long run high grey level emphasis
glrlm.Frlm_lrhge = sum(sum(GLRLM.*(rMat.^2).*(cMat.^2)));

% Gray level non-uniformity
glrlm.Frlm_glnu = sum(pg.^2) * Ns;

% Gray level non-uniformity normalised
glrlm.Frlm_glnu_norm = sum(pg.^2);

% Run length non-uniformity
glrlm.Frlm_rlnu = sum(pr.^2) * Ns;

% Run length non-uniformity normalised
glrlm.Frlm_rlnu_norm = sum(pr.^2);

% Run percentage
glrlm.Frlm_r_perc = sum(pg)/(pr*cVect');

% Grey level variance
temp = rMat .* GLRLM;
u = sum(temp(:));
temp = ((rMat - u).^2) .* GLRLM;
glrlm.Frlm_gl_var = sum(temp(:));

% Run length variance
temp = cMat .* GLRLM;
u = sum(temp(:));
temp = ((cMat - u).^2) .* GLRLM;
glrlm.Frlm_rl_var = sum(temp(:));

% Run entropy
valPos = GLRLM(find(GLRLM(:)));
temp = valPos .* log2(valPos);
glrlm.Frlm_rl_entr = -sum(temp(:));

end