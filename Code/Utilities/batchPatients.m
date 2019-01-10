function [patients] = batchPatients(nPatient,nBatch)
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


% FIND THE NUMBER OF PATIENTS IN EACH BATCH
patients = cell(1,nBatch);
patientVect = randperm(nPatient); % To randomize stuff a bit.
if nBatch
    nP = nPatient / nBatch; nSup = ceil(nP); nInf = floor(nP);
    if nSup ~= nInf
        nSubInf = nBatch - 1; nSubSup = 1; total = nSubInf*nInf + nSubSup*nSup;
        while total ~= nPatient
            nSubInf = nSubInf - 1; nSubSup = nSubSup + 1;
            total = nSubInf*nInf + nSubSup*nSup;
        end
        nP = [repmat(nInf,[1,nSubInf]),repmat(nSup,[1,nSubSup])];
    else % The number of patients in all batches will be the same
        nP = repmat(nSup,[1,nBatch]);
    end
    start = 1;
    for i = 1:nBatch
        patients{i} = patientVect(start:(start+nP(i)-1));
        start = start+nP(i);
    end
end

end
