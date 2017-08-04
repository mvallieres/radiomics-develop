function computeAllRadiomics(pathRead,pathSave,nameRead,nameROI,nameSet,imParams,roiType)
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

startpath = pwd;

nPatient = numel(nameRead);

% COMPUTATION
fprintf('\n')
for p = 1:nPatient
    fprintf(['\n*********************** COMPUTING FEATURES: %s, ***********************'],nameRead{p});
    tStart = tic;
    cd(pathRead)
    load(nameRead{p}) % Variable 'sData' now in MATLAB Workspace;
    scanType = sData{2}.type;
    imParamScan = imParams.(scanType);
    if strcmp(scanType,'PTscan') % MOVE THIS IN THE READIND DATA PART AND ALSO APPLY PVE CORRECTIONS AND DENOISING (for MRI(N3/N4 method) and CT too)? YES, THIS IS NEEDED ESPECIALLY FOR MRI AND PET. FUTURE WORK ON ITS WAY (image post-processing).
        sData{2}.scan.volume.data = computeSUVmap(sData{2}.scan.volume.data,sData{3}(1));
    end
    [radiomics] = computePatientRadiomics(sData,nameROI{p},nameSet{p},imParamScan);
    indDot = strfind(nameRead{p},'.');
    nameSave = [nameRead{p}(1:(indDot(1)-1)),'(',roiType,')',nameRead{p}(indDot(1):end)];
    cd(pathSave), save(nameSave,'radiomics')
    time = toc(tStart);
    fprintf('TOTAL TIME: %.2f seconds\n',time)
end

cd(startpath)
end