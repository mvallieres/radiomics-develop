function p_iminusj = gclm_DiagProb(p_ij)
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

Ng = size(p_ij,1);
valK = 0:(Ng-1); nK = numel(valK);
p_iminusj = zeros(nK,1);

for iterationK = 1:nK
    k = valK(iterationK);
    p = 0;
    for i = 1:Ng
        for j = 1:Ng
            if (k - abs(i-j)) == 0
                p = p + p_ij(i,j);
            end
        end
    end
    p_iminusj(iterationK) = p;
end

end