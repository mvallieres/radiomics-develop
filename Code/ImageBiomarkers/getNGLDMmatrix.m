function [NGLDM] = getNGLDMmatrix(ROIonly,levels)
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

% FOR a = 0, d = 1.


% PRELIMINARY
nLevel = length(levels);
if nLevel > 100
    adjust = 10000;
else
    adjust = 1000;
end
levelTemp = max(levels)+1;
ROIonly(isnan(ROIonly)) = levelTemp;
levels = [levels,levelTemp];

dim = size(ROIonly);
if ndims(ROIonly) == 2
	dim(3) = 1;
end
q2 = reshape(ROIonly,1,prod(dim));


% % QUANTIZATION EFFECTS CORRECTION (M. Vallieres)
% % In case (for example) we initially wanted to have 64 levels, but due to
% % quantization, only 60 resulted.
% qs = round(levels*adjust)/adjust;
% q2 = round(q2*adjust)/adjust;
qs = levels;


% EL NAQA CODE
q3 = q2*0;
for k = 1:length(qs)
	q3(q2==qs(k)) = k;
end
q3 = reshape(q3,dim);
lqs = length(qs);
NGLDM = double(zeros(lqs,27)); % Min dependence = 0, Max dependence = 26; So 27 columns
for i = 1:dim(1)
    i_min = max(1,i-1);
    i_max = min(i+1,dim(1));
	for j = 1:dim(2)
        j_min = max(1,j-1);
		j_max = min(j+1,dim(2));
		for k = 1:dim(3)
			k_min = max(1,k-1);
			k_max = min(k+1,dim(3));
            val_q3 = q3(i,j,k);
            count = 0;
			for I2 = i_min:i_max
				for J2 = j_min:j_max
					for K2 = k_min:k_max
						if I2 == i && J2 == j && K2 == k
							continue;
                        else
                            val_neighbor = q3(I2,J2,K2);
                            if (val_q3 - val_neighbor == 0) % a = 0
                                count = count + 1;
                            end
						end
					end
				end
            end
            NGLDM(val_q3,count + 1) = NGLDM(val_q3,count + 1) + 1;
		end
	end
end

NGLDM(end,:) = []; % Last column was for the NaN voxels, to be removed
stop = find(sum(NGLDM),1,'last');
NGLDM(:,(stop+1):end) = [];

end
