function morph = getMorphFeatures(vol,maskInt,maskMorph,res,intensity)
% -------------------------------------------------------------------------
% AUTHOR(S): 
% - Martin Vallieres <mart.vallieres@gmail.com>
% -------------------------------------------------------------------------
% HISTORY:
% - Creation: April 2017
% - Revision I: August 2017
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

% - vol: 3D volume, NON-QUANTIZED, continous imaging intensity distribution
% - maskInt: Intensity mask
% - maskMorph: Morphological mask
% - res: [a,b,c] vector specfying the resolution of the volume in mm.  % XYZ resolution (world), or JIK resolution (intrinsic matlab).
% - intensity (optional): If 'arbitrary', some feature will not be computed. If
%   'definite', all feature will be computed. If not present as an argument,
%   all features will be computed. If 'filter', most features will not be
%   computed, except some. The 'filter' option encompasses 'arbitrary',
%   must is even more stringent. Please see below.
% REFERENCES
% [1] https://arxiv.org/abs/1612.07003 (make it formal with authors, etc.)


% INTIALIZATION
if nargin < 5
    definite = true;
    filter = false;
else
    if strcmp(intensity,'arbitrary')
        definite = false;
        filter = false;
    elseif strcmp(intensity,'definite')
        definite = true;
        filter = false;
    elseif strcmp(intensity,'filter')
        definite = false;
        filter = true;
    else
        error('Fifth argument must either be "arbitrary" or "definite" or "filter"')
    end 
end


% PADDING THE VOLUME WITH A LAYER OF NaNs (reduce mesh computation errors of associated mask)
vol = padarray(vol,[1,1,1],NaN);
% PADDING THE MASKS WITH A LAYER OF 0's (reduce mesh computation errors of associated mask)
maskInt = padarray(maskInt,[1,1,1],0);
maskMorph = padarray(maskMorph,[1,1,1],0);

% GETTING IMPORTANT VARIABLES
Xgl_int = vol(maskInt == 1);
Xgl_morph = vol(maskMorph == 1);
[XYZ_int,~,~] = getMesh(maskInt,res); % XYZ refers to [Xc,Yc,Zc] in ref. [1].
[XYZ_morph,faces,vertices] = getMesh(maskMorph,res); % XYZ refers to [Xc,Yc,Zc] in ref. [1].
convHull = convhull(vertices(:,1),vertices(:,2),vertices(:,3)); % [X,Y,Z] points of the convex hull.


% STARTING COMPUTATION

% Volume
if ~filter
    volume = getMeshVolume(faces,vertices); % In mm^3
    morph.Fmorph_volume = volume;
else
    morph.Fmorph_volume = [];
end

% Approximate Volume
if ~filter
    morph.Fmorph_approx_volume = sum(maskMorph(:)) * prod(res);
else
    morph.Fmorph_approx_volume = [];
end

% Surface area
if ~filter
    area = getMeshArea(faces,vertices); % In mm^2
    morph.Fmorph_area = area;
else
    morph.Fmorph_area = [];
end

% Surface to volume ratio
if ~filter
    morph.Fmorph_av = area/volume;
else
    morph.Fmorph_av = [];
end

% Compactness 1
if ~filter
    morph.Fmorph_comp1 = volume/((pi^(1/2))*(area^(3/2)));
else
    morph.Fmorph_comp1 = [];
end

% Compactness 2
if ~filter
    morph.Fmorph_comp2 = 36*pi*(volume^2)/(area^3);
else
    morph.Fmorph_comp2 = [];
end

% Spherical disproportion
if ~filter
    morph.Fmorph_sph_dispr = area/(36*pi*volume^2)^(1/3);
else
    morph.Fmorph_sph_dispr = [];
end

% Sphericity
if ~filter
    morph.Fmorph_sphericity = ((36*pi*volume^2)^(1/3))/area;
else
    morph.Fmorph_sphericity = [];
end

% Asphericity
if ~filter
    morph.Fmorph_asphericity = ((area^3)/(36*pi*volume^2))^(1/3) - 1;
else
    morph.Fmorph_asphericity = [];
end

% Centre of mass shift
morph.Fmorph_com = getCOM(Xgl_int,Xgl_morph,XYZ_int,XYZ_morph);

% Maximum 3D diameter
if ~filter
    morph.Fmorph_diam = getMax3Ddiam(convHull,vertices);
else
    morph.Fmorph_diam = [];
end

% Major axis length
if ~filter
    [major,minor,least] = getAxisLengths(XYZ_morph);
    morph.Fmorph_pca_major = 4*sqrt(major);
else
    morph.Fmorph_pca_major = [];
end

% Minor axis length
if ~filter
    morph.Fmorph_pca_minor = 4*sqrt(minor);
else
    morph.Fmorph_pca_minor = [];
end

% Least axis length
if ~filter
    morph.Fmorph_pca_least = 4*sqrt(least);
else
    morph.Fmorph_pca_least = [];
end

% Elongation
if ~filter
    morph.Fmorph_pca_elongation = sqrt(minor/major);
else
    morph.Fmorph_pca_elongation = [];
end

% Flatness
if ~filter
    morph.Fmorph_pca_flatness = sqrt(least/major);
else
    morph.Fmorph_pca_flatness = [];
end

% Volume density - axis-aligned bounding box
if ~filter
    Xc_aabb = max(vertices(:,1)) - min(vertices(:,1)); Yc_aabb = max(vertices(:,2)) - min(vertices(:,2)); Zc_aabb = max(vertices(:,3)) - min(vertices(:,3));
    Vaabb = Xc_aabb * Yc_aabb * Zc_aabb; 
    morph.Fmorph_v_dens_aabb = volume / Vaabb;
else
    morph.Fmorph_v_dens_aabb = [];
end

% Area density - axis-aligned bounding box
if ~filter
    Aaabb = 2*Xc_aabb*Yc_aabb + 2*Xc_aabb*Zc_aabb + 2*Yc_aabb*Zc_aabb;
    morph.Fmorph_a_dens_aabb = area / Aaabb;
else
    morph.Fmorph_a_dens_aabb = [];
end

% Volume density - oriented minimum bounding box
% Copyright (c) 2015: Johannes Korsawe, 2012: John D'Errico
% https://www.mathworks.com/matlabcentral/fileexchange/18264-minimal-bounding-box
if ~filter
    [~,~,Vombb,Aombb,~] = minboundbox(vertices(:,1),vertices(:,2),vertices(:,3),'v',3);
    morph.Fmorph_v_dens_ombb = volume / Vombb;
else
    morph.Fmorph_v_dens_ombb = [];
end

% Area density - oriented minimum bounding box
if ~filter
    morph.Fmorph_a_dens_ombb = area / Aombb;
else
    morph.Fmorph_a_dens_ombb = [];
end

% Volume density - approximate enclosing ellipsoid
if ~filter
    a = 2*sqrt(major); b = 2*sqrt(minor); c = 2*sqrt(least);
    Vaee = (4*pi*a*b*c)/3;
    morph.Fmorph_v_dens_aee = volume / Vaee;
else
    morph.Fmorph_v_dens_aee = [];
end

% Area density - approximate enclosing ellipsoid
if ~filter
    Aaee = getAreaDensApprox(a,b,c,20);
    morph.Fmorph_a_dens_aee = area / Aaee;
else
    morph.Fmorph_a_dens_aee = [];
end

% Volume density - minimum volume enclosing ellipsoid (Rotate the volume
% first??)
% Copyright (c) 2009, Nima Moshtagh
% http://www.mathworks.com/matlabcentral/fileexchange/9542-minimum-volume-enclosing-ellipsoid
%[A,~] = MinVolEllipse([vertices(:,1),vertices(:,2),vertices(:,3)]',0.01); % Subsequent singular value decomposition of matrix A and and taking the inverse of the square root of the diagonal of the sigma matrix will produce respective semi-axis lengths.
if ~filter
    [A,~] = MinVolEllipse([vertices(convHull(:,1),1),vertices(convHull(:,2),2),vertices(convHull(:,3),3)]',0.01); % Subsequent singular value decomposition of matrix A and and taking the inverse of the square root of the diagonal of the sigma matrix will produce respective semi-axis lengths.
    [~,Q,~] = svd(A); a = 1/sqrt(Q(3,3)); b = 1/sqrt(Q(2,2)); c = 1/sqrt(Q(1,1)); % New semi-axis lengths
    Vmvee = (4*pi*a*b*c)/3;
    morph.Fmorph_v_dens_mvee = volume / Vmvee;
else
    morph.Fmorph_v_dens_mvee = [];
end

% Area density - minimum volume enclosing ellipsoid
if ~filter
    Amvee = getAreaDensApprox(a,b,c,20); % Using a new set of (a,b,c), see Volume density - minimum volume enclosing ellipsoid
    morph.Fmorph_a_dens_mvee = area / Amvee;
else
    morph.Fmorph_a_dens_mvee = [];
end

% Volume density - convex hull
if ~filter
    Vconvex = getMeshVolume(convHull,vertices);
    morph.Fmorph_v_dens_conv_hull = volume / Vconvex;
else
    morph.Fmorph_v_dens_conv_hull = [];
end

% Area density - convex hull
if ~filter
    Aconvex = getMeshArea(convHull,vertices);
    morph.Fmorph_a_dens_conv_hull = area / Aconvex;
else
    morph.Fmorph_a_dens_conv_hull = [];
end

% Integrated intensity
if definite
    morph.Fmorph_integ_int = mean(Xgl_int) * volume;
else
    morph.Fmorph_integ_int = [];
end

% % Moran's I index
% vol(maskInt == 0) = NaN;
% morph.Fmorph_moran_i = getMoranI(vol,res); % NEEDS TO BE VECTORIZED FOR FASTER CALCULATION!
% 
% % Geary's C measure
% morph.Fmorph_geary_c = getGearyC(vol,res); % NEEDS TO BE VECTORIZED FOR FASTER CALCULATION!
morph.Fmorph_moran_i = [];
morph.Fmorph_geary_c = [];


end