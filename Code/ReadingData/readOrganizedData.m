function readOrganizedData(pathPatient,pathSave,namePatient)
% -------------------------------------------------------------------------
% AUTHOR(S): 
% - Martin Vallieres <mart.vallieres@gmail.com>
% -------------------------------------------------------------------------
% HISTORY:
% - Creation: March 2018
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

% INITIALIZATION
cd(pathPatient), listScans = dir;
listScans = listScans(~ismember({listScans.name},{'.','..','.DS_Store','._.DS_Store'}));
nScans = numel(listScans);


for s = 1:nScans
    scan = listScans(s).name;
    cd(fullfile(pathPatient,scan))
    ok = zeros(1,3);
    
    % STEP 1: READING DICOM DATA (there must be DICOM data present in the directory for now --> minimim requirement is to report on image acquisition parameters).
    try
        readAllDICOM(pwd,pwd,0,'matlab');
        listMat = dir('*.mat'); % There must be only one scan series in the organized patient-scan folder
        sData = load(listMat(1).name); sData = struct2cell(sData); sData = sData{1};
        delete('*.mat')
        sData{2}.scan.volume.data = uint16(sData{2}.scan.volume.data); % Usually the native format of DICOM. Could be int16 though, need to check on how to automatically detect that.
        nameSave = [namePatient,'_',scan,'.',sData{2}.type,'.mat'];
        ok(1) = 1;
    catch % PROBLEM WITH DICOM DATA. CODE DOWN THE ROAD WILL FAIL. TO MANUALLY CORRECT. DICOM DATA MUST BE PRESENT AND READABLE.
        sData = cell(1,7); % MINIMAL sData
        sData{2} = struct;
        nameSave = [namePatient,'_',scan,'.Unknown.mat'];
    end
    
    
    
    % STEP 2: VERIFY FOR THE PRESENCE OF .nnrd files
    listNRRD = dir('*.nrrd');
    if ~isempty(listNRRD)
        if exist('imagingVolume.nrrd') % We recommend to always provide "imagingVolume.nrrd" if .nrrd segmentation is performed. However, it may still happen that the only imagin volume present comes from the DICOM data.
            ok(2) = 1;
            [sData{2}.nrrd.volume.data,sData{2}.nrrd.volume.header] = nrrdread('imagingVolume.nrrd');
        else
            sData{2}.nrrd.volume.data = sData{2}.scan.volume.data; % Copied from DICOM data for now.
            sData{2}.nrrd.volume.header = 'None -- Imaging volume = DICOM data';
        end
        listMask = dir('segMask*.nrrd'); nMask = numel(listMask);
        sData{2}.nrrd.mask = struct;
        for m = 1:nMask
            nameMaskFile = listMask(m).name;
            indUnderScore = strfind(nameMaskFile,'_');
            indDot = strfind(nameMaskFile,'.');
            if numel(indUnderScore) == 2
                nameROI = nameMaskFile((indUnderScore(1)+1):(indUnderScore(2)-1));
                labelROI = str2num(nameMaskFile((indUnderScore(2)+6):(indDot(1)-1)));
            else % Then there must be only one underscore
                nameROI = nameMaskFile((indUnderScore(1)+1):(indDot(1)-1));
                labelROI = 1;
            end
            sData{2}.nrrd.mask(m).name = nameROI;
            [sData{2}.nrrd.mask(m).data,sData{2}.nrrd.mask(m).header] = nrrdread(nameMaskFile);
            sData{2}.nrrd.mask(m).data(sData{2}.nrrd.mask(m).data ~= labelROI) = NaN;
            sData{2}.nrrd.mask(m).data(sData{2}.nrrd.mask(m).data == labelROI) = 1;
            sData{2}.nrrd.mask(m).data(sData{2}.nrrd.mask(m).data ~= labelROI) = 0;
            sData{2}.nrrd.mask(m).data = uint16(sData{2}.nrrd.mask(m).data); % Apparently, this is the native format of .nrrd
       end
    end
   
    % STEP 3: VERIFY FOR THE PRESENCE OF .img files
    listIMG = dir('*.img');
    if ~isempty(listIMG)
        if exist('imagingVolume.img')
            ok(3) = 1;
            sData{2}.img.volume.data = niftiread('imagingVolume');
            sData{2}.img.volume.header = niftiinfo('imagingVolume');
            listMask = dir('segMask*.img'); nMask = numel(listMask);
            sData{2}.img.mask = struct;
            for m = 1:nMask
                nameMaskFile = listMask(m).name;
                indUnderScore = strfind(nameMaskFile,'_');
                indDot = strfind(nameMaskFile,'.');
                if numel(indUnderScore) == 2
                    nameROI = nameMaskFile((indUnderScore(1)+1):(indUnderScore(2)-1));
                    labelROI = str2num(nameMaskFile((indUnderScore(2)+6):(indDot(1)-1)));
                else % Then there must be only one underscore
                    nameROI = nameMaskFile((indUnderScore(1)+1):(indDot(1)-1));
                    labelROI = 1;
                end
                sData{2}.img.mask(m).name = nameROI;
                sData{2}.img.mask(m).data = niftiread(nameMaskFile(1:end-4));
                sData{2}.img.mask(m).header = niftiinfo(nameMaskFile(1:end-4));
                sData{2}.img.mask(m).data(sData{2}.img.mask(m).data ~= labelROI) = NaN;
                sData{2}.img.mask(m).data(sData{2}.img.mask(m).data == labelROI) = 1;
                sData{2}.img.mask(m).data(isnan(sData{2}.img.mask(m).data)) = 0;
                sData{2}.img.mask(m).data = single(sData{2}.img.mask(m).data); % % Apparently, this is the native format of .img
            end
        end
    end
   
    % STEP 4: SAVE FINAL "sData" FILE
    if sum(ok) % Otherwise, no imaging volume was present or without error, so there is nothing to save.
        cd(pathSave), save(nameSave,'sData','-v7.3')
    end
end

cd(startpath)
end
