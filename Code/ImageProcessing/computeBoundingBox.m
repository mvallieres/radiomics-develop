function [boxBound] = computeBoundingBox(mask)
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


[iV,jV,kV] = find3d(mask);
boxBound(1,1) = min(iV);
boxBound(1,2) = max(iV);
boxBound(2,1) = min(jV);
boxBound(2,2) = max(jV);
boxBound(3,1) = min(kV);
boxBound(3,2) = max(kV);

end


% CERR UTILITY FUNCTIONS (can be found at: https://github.com/adityaapte/CERR)
function [iV,jV,kV] = find3d(mask3M)
indV = find(mask3M(:));
[iV,jV,kV] = fastind2sub(size(mask3M),indV);
iV = iV';
jV = jV';
kV = kV';
end

function varargout = fastind2sub(siz,ndx)
nout = max(nargout,1);
if length(siz)<=nout,
  siz = [siz ones(1,nout-length(siz))];
else
  siz = [siz(1:nout-1) prod(siz(nout:end))];
end
n = length(siz);
k = [1 cumprod(siz(1:end-1))];
ndx = ndx - 1;
for i = n:-1:1,
  varargout{i} = floor(ndx/k(i)) + 1;
  ndx = ndx - (varargout{i}-1) * k(i);
end
end