function Aell = getAreaDensApprox(a,b,c,n)
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

% - a: Major semi-axis length
% - b: Minor semi-axis length
% - c: Least semi-axis length
% - n: Number of iterations

alpha = sqrt(1-b^2/a^2); beta = sqrt(1-c^2/a^2); 
AB = alpha*beta; point = (alpha^2 + beta^2) / (2*AB);
Aell = 0;
for v = 0:n
    Aell = Aell + AB^v / (1-4*v^2) * legendreP(v,point);
end
Aell = Aell * 4 * pi * a * b;

end