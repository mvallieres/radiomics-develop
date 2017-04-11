function morph = getMorphFeatures(vol,res)
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

% - vol: 3D volume, NON-QUANTIZED, with NaNs outside the region of interest
% --> vol: continous imaging intensity distribution
% - res: [a,b,c] vector specfying the resolution of the volume in mm.  % XYZ resolution (world), or JIK resolution (intrinsic matlab).
%
% REFERENCES
% [1] https://arxiv.org/abs/1612.07003 (make it formal with authors, etc.)


% PADDING THE VOLUME WITH A LAYER OF NaNs (reduce mesh computation errors of associated mask)
vol = padarray(vol,[1,1,1],NaN);

% GETTING IMPORTANT VARIABLES
Xgl = vol(~isnan(vol(:)));
mask = vol; mask(~isnan(mask)) = 1; mask(isnan(mask)) = 0;
[XYZ,faces,vertices] = getMesh(mask,res); % XYZ refers to [Xc,Yc,Zc] in ref. [1].
convHull = convhull(vertices(:,1),vertices(:,2),vertices(:,3)); % [X,Y,Z] points of the convex hull.


% STARTING COMPUTATION

% Volume
volume = getMeshVolume(faces,vertices); % In mm^3
morph.Fmorph_volume = volume;

% Approximate Volume
morph.Fmorph_approx_volume = sum(mask(:)) * prod(res);

% Surface area
area = getMeshArea(faces,vertices); % In mm^2
morph.Fmorph_area = area;

% Surface to volume ratio
morph.Fmorph_av = area/volume;

% Compactness 1
morph.Fmorph_comp1 = volume/((pi^(1/2))*(area^(3/2)));

% Compactness 2
morph.Fmorph_comp2 = 36*pi*(volume^2)/(area^3);

% Spherical disproportion
morph.Fmorph_sph_dispr = area/(36*pi*volume^2)^(1/3);

% Sphericity
morph.Fmorph_sphericity = ((36*pi*volume^2)^(1/3))/area;

% Asphericity
morph.Fmorph_asphericity = ((area^3)/(36*pi*volume^2))^(1/3) - 1;

% Centre of mass shift
morph.Fmorph_com = getCOM(Xgl,XYZ);

% Maximum 3D diameter
morph.Fmorph_diam = getMax3Ddiam(convHull,vertices);

% Major axis length
[major,minor,least] = getAxisLengths(XYZ);
morph.Fmorph_pca_major = 4*sqrt(major);

% Minor axis length
morph.Fmorph_pca_minor = 4*sqrt(minor);

% Least axis length
morph.Fmorph_pca_least = 4*sqrt(least);

% Elongation
morph.Fmorph_pca_elongation = sqrt(minor/major);

% Flatness
morph.Fmorph_pca_flatness = sqrt(least/major);

% Volume density - axis-aligned bounding box
Xc_aabb = max(vertices(:,1)) - min(vertices(:,1)); Yc_aabb = max(vertices(:,2)) - min(vertices(:,2)); Zc_aabb = max(vertices(:,3)) - min(vertices(:,3));
Vaabb = Xc_aabb * Yc_aabb * Zc_aabb; 
morph.Fmorph_v_dens_aabb = volume / Vaabb;

% Area density - axis-aligned bounding box
Aaabb = 2*Xc_aabb*Yc_aabb + 2*Xc_aabb*Zc_aabb + 2*Yc_aabb*Zc_aabb;
morph.Fmorph_a_dens_aabb = area / Aaabb;

% Volume density - oriented minimum bounding box
% Copyright (c) 2015: Johannes Korsawe, 2012: John D'Errico
% https://www.mathworks.com/matlabcentral/fileexchange/18264-minimal-bounding-box
[~,~,Vombb,Aombb,~] = minboundbox(vertices(:,1),vertices(:,2),vertices(:,3),'v',3);
morph.Fmorph_v_dens_ombb = volume / Vombb;

% Area density - oriented minimum bounding box
morph.Fmorph_a_dens_ombb = area / Aombb;

% Volume density - approximate enclosing ellipsoid
a = 2*sqrt(major); b = 2*sqrt(minor); c = 2*sqrt(least);
Vaee = (4*pi*a*b*c)/3;
morph.Fmorph_v_dens_aee = volume / Vaee;

% Area density - approximate enclosing ellipsoid
Aaee = getAreaDensApprox(a,b,c,20);
morph.Fmorph_a_dens_aee = area / Aaee;

% Volume density - minimum volume enclosing ellipsoid
% Copyright (c) 2009, Nima Moshtagh
% http://www.mathworks.com/matlabcentral/fileexchange/9542-minimum-volume-enclosing-ellipsoid
%[A,~] = MinVolEllipse([vertices(:,1),vertices(:,2),vertices(:,3)]',0.01); % Subsequent singular value decomposition of matrix A and and taking the inverse of the square root of the diagonal of the sigma matrix will produce respective semi-axis lengths.
[A,~] = MinVolEllipse([vertices(convHull(:,1),1),vertices(convHull(:,2),2),vertices(convHull(:,3),3)]',0.01); % Subsequent singular value decomposition of matrix A and and taking the inverse of the square root of the diagonal of the sigma matrix will produce respective semi-axis lengths.
[~,Q,~] = svd(A); a = 1/sqrt(Q(3,3)); b = 1/sqrt(Q(2,2)); c = 1/sqrt(Q(1,1)); % New semi-axis lengths
Vmvee = (4*pi*a*b*c)/3;
morph.Fmorph_v_dens_mvee = volume / Vmvee;

% Area density - minimum volume enclosing ellipsoid
Amvee = getAreaDensApprox(a,b,c,20); % Using a new set of (a,b,c), see Volume density - minimum volume enclosing ellipsoid
morph.Fmorph_a_dens_mvee = area / Amvee;

% Volume density - convex hull
Vconvex = getMeshVolume(convHull,vertices);
morph.Fmorph_v_dens_conv_hull = volume / Vconvex;

% Area density - convex hull
Aconvex = getMeshArea(convHull,vertices);
morph.Fmorph_a_dens_conv_hull = area / Aconvex;

% Integrated intensity
morph.Fmorph_integ_int = mean(Xgl) * volume;

% Moran's I index
%morph.Fmorph_moran_i = getMoranI(vol,res); % NEEDS TO BE VECTORIZED FOR FASTER CALCULATION!

% Geary's C measure
%morph.Fmorph_geary_c = getGearyC(vol,res); % NEEDS TO BE VECTORIZED FOR FASTER CALCULATION!

end