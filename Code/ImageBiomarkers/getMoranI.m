function moranI= getMoranI(vol,res)
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


% Find the location(s) of all non NaNs voxels
ind = find(~isnan(vol));
[I,J,K] = ind2sub(size(vol),ind);
nVox = numel(I);

% Get the mean
u = mean(vol(ind));
volMmean = vol - u; % (Xgl,i - u)
volMmeanS = (vol - u).^2; % % (Xgl,i - u).^2
sumS = sum(volMmeanS(~isnan(volMmeanS(:)))); % Sum of (Xgl,i - u).^2 over all i

% Get a meshgrid first
[X,Y,Z] = meshgrid(res(1).*((1:size(vol,2))-0.5),res(2).*((1:size(vol,1))-0.5),res(3).*((1:size(vol,3))-0.5));


temp = 0;
sumW = 0;
for i = 1:nVox
   
    % Distance mesh
    tempX = X - X(I(i),J(i),K(i));
    tempY = Y - Y(I(i),J(i),K(i));
    tempZ = Z - Z(I(i),J(i),K(i));
    tempDistMesh = 1./sqrt(tempX.^2 + tempY.^2 + tempZ.^2); % meshgrid of weigths
    tempDistMesh(isnan(vol)) = NaN; tempDistMesh(I(i),J(i),K(i)) = NaN; % Removing NaNs
    sumW = sumW + sum(tempDistMesh(~isnan(tempDistMesh(:)))); % Running sum of weights
    
    % Inside sum calculation
    tempVol = volMmean; tempVol(I(i),J(i),K(i)) = NaN; % Removing NaNs
    tempVol = tempDistMesh .* tempVol; % (wij .* (Xgl,j - u))
    sumVal = sum(tempVol(~isnan(tempVol(:)))); % Summing (wij .* (Xgl,j - u)) over all j
    
    temp = temp + volMmean(I(i),J(i),K(i)) * sumVal; % Running sum of (Xgl,i - u)*(wij .* (Xgl,j - u)) over all i
end

moranI = temp*nVox/sumS/sumW;
end