function glcm = getGLCMfeatures(vol)
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


% GET THE GLCM MATRIX
levels = 1:max(vol(~isnan(vol(:)))); % Correct definition, without any assumption
[GLCM] = getGLCMmatrix(vol,levels);
p_ij = GLCM./(sum(GLCM(:))); % Normalization of GLCM
p_i = sum(p_ij,2); p_j = sum(p_ij);
p_iminusj = gclm_DiagProb(p_ij);
p_iplusj = gclm_CrossDiagProb(p_ij);
Ng = max(size(GLCM));
vectNg = 1:Ng;
[colGrid,rowGrid] = meshgrid(vectNg,vectNg);



% COMPUTING TEXTURES

% Joint maximum
glcm.Fcm_joint_max = max(p_ij(:));

% Joint average
temp = rowGrid.*p_ij;
u = sum(temp(:));
glcm.Fcm_joint_avg = u; 

% Joint variance
temp = ((rowGrid - u).^2) .* p_ij;
var = sum(temp(:));
glcm.Fcm_joint_var = var;

% Joint entropy
pPos = p_ij(find(p_ij(:))); % Exclusing those with 0 probability
temp = pPos.*log2(pPos);
glcm.Fcm_joint_entr = -sum(temp(:));

% Difference average
k = 0:(Ng-1);
u = k * p_iminusj;
glcm.Fcm_diff_avg = u;

% Difference variance
var = ((k - u).^2) * p_iminusj;
glcm.Fcm_diff_var = var;

% Difference entropy
kPos = p_iminusj(find(p_iminusj(:)));
glcm.Fcm_diff_entr = - (kPos' * log2(kPos));

% Sum average
k = 2:(Ng*2);
u = k * p_iplusj;
glcm.Fcm_sum_avg = u;

% Sum variance
var = ((k - u).^2) * p_iplusj;
glcm.Fcm_sum_var = var;

% Sum entropy
kPos = p_iplusj(find(p_iplusj(:)));
glcm.Fcm_sum_entr = - (kPos' * log2(kPos));

% Angular second moment
temp = p_ij.^2;
glcm.Fcm_energy = sum(temp(:));

% Contrast
temp = (rowGrid - colGrid).^2 .* p_ij;
glcm.Fcm_contrast = sum(temp(:));

% Dissimilarity
temp = abs(rowGrid - colGrid) .* p_ij;
glcm.Fcm_dissimilarity = sum(temp(:));

% Inverse difference
temp = p_ij./(1 + abs(rowGrid - colGrid));
glcm.Fcm_inv_diff = sum(temp(:));

% Inverse difference normalised
temp = p_ij./(1 + abs(rowGrid - colGrid)/Ng);
glcm.Fcm_inv_diff_norm = sum(temp(:));

% Inverse difference moment
temp = p_ij./(1 + (rowGrid - colGrid).^2);
glcm.Fcm_inv_diff_mom = sum(temp(:));

% Inverse difference moment normalised
temp = p_ij./(1 + ((rowGrid - colGrid).^2)/Ng^2);
glcm.Fcm_inv_diff_mom_norm = sum(temp(:));

% Inverse variance
p = 0;
for i = 1:Ng
    for j = (i+1):Ng
        p = p + p_ij(i,j)/((i-j)^2);
    end
end
glcm.Fcm_inv_var = 2*p;

% Correlation
u_i = vectNg*p_i;
u_j = vectNg*p_j';
std_i = sqrt(((vectNg - u_i).^2) *  p_i);
std_j = sqrt(((vectNg - u_j).^2) *  p_j');
temp = rowGrid .* colGrid .* p_ij;
glcm.Fcm_corr = (1/(std_i*std_j)) * (-u_i*u_j + sum(temp(:)));

% Autocorrelation
temp = rowGrid .* colGrid .* p_ij;
glcm.Fcm_auto_corr = sum(temp(:));

% Cluster tendency
temp = ((rowGrid + colGrid - u_i - u_j).^2) .* p_ij;
glcm.Fcm_clust_tend = sum(temp(:));

% Cluster shade
temp = ((rowGrid + colGrid - u_i - u_j).^3) .* p_ij;
glcm.Fcm_clust_shade = sum(temp(:));

% Cluster prominence
temp = ((rowGrid + colGrid - u_i - u_j).^4) .* p_ij;
glcm.Fcm_clust_prom = sum(temp(:));

% First measure of information correlation
pPos = p_ij(find(p_ij(:))); temp = pPos.*log2(pPos); HXY = -sum(temp(:));
pPos = p_i(find(p_i(:))); temp = pPos.*log2(pPos); HX = -sum(temp(:));
p_i_temp = repmat(p_i,[1,Ng]); p_j_temp = repmat(p_j,[Ng,1]);
p_temp = p_i_temp .* p_j_temp;
pPos = p_ij(find(p_temp(:))); pPos_temp = p_temp(find(p_temp(:)));
temp = pPos .* log2(pPos_temp);
HXY1 = -sum(temp(:));
glcm.Fcm_info_corr1 = (HXY - HXY1)/HX;

% Second measure of information correlation
temp = pPos_temp .* log2(pPos_temp);
HXY2 = -sum(temp(:));
if HXY > HXY2
    glcm.Fcm_info_corr2 = 0;
else
    glcm.Fcm_info_corr2 = sqrt(1 - exp(-2*(HXY2 - HXY)));
end

end