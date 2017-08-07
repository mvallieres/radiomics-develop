function flagRadiomics = roiChoice(pathWORK,pathDATA)
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


startpath = pwd;

% List of MATLAB imaging data files
cd(pathDATA), list = dir('*.mat'); nFile = numel(list); roiCheck = {};
close all

% MANUAL CHOICE OF ROI
fprintf('\n\n%%%%%%%%%%%% MANUAL CHOICE OF ROIs %%%%%%%%%%%%')
for f = 1:nFile
    roiOK = 0;
    name = list(f).name;
    [ID,scan,type] = findIDscanType(name); roiNames{f,1} = ID; roiNames{f,2} = scan; roiNames{f,3} = type;
    sData = load(name); sData = struct2cell(sData); sData = sData{1};
    
    while ~roiOK
        nROI = numel(sData{2}.scan.contour);
        
        % CHOICE IN THE LIST
        fprintf('\n----- LIST OF ROIs FOR "%s" -----',name);
        for r = 1:nROI
            fprintf(['\n',num2str(r),'. ',sData{2}.scan.contour(r).name,' -- ',sData{2}.scan.contour(r).nameSet])
        end
        choiceOK = 0;
        fprintf('\n')
        while ~choiceOK
            roi = input('--> WHICH ROI DO YOU CHOOSE FOR THAT SCAN? Example --> Enter 1 for roi #1, [2,4] for roi #2 and #4: ');
            if isnumeric(roi) && ~sum((roi - nROI) > 0) && ~sum((1 - roi) > 0)
                choiceOK = 1;
            end
        end
        nROI = numel(roi);
        nameROI = []; nameSet = [];
        for r = 1:nROI
            nameROI = [nameROI,sData{2}.scan.contour(roi(r)).name,','];
            nameSet = [nameSet,sData{2}.scan.contour(roi(r)).nameSet,','];
        end
        nameROI = nameROI(1:end-1); nameSet = nameSet(1:end-1); % Removing the comma at the end
        roiNames{f,4} = nameROI; roiNames{f,5} = nameSet;

        % ROI COMPUTATION
        fprintf('Computing the ROI ... ')
        contourNumber = findContour(sData,nameROI,nameSet);
        [volObjInit,roiObjInit] = getROI(sData,contourNumber,'full');

        % FINDING MIDDLE SLICES (axial, sagittal, coronal)
        sizeM = size(roiObjInit.data);
        [J,I,K] = meshgrid(1:sizeM(2),1:sizeM(1),1:sizeM(3)); % IJK in MATLAB coordinates
        indMask = find(roiObjInit.data);
        I = I(indMask); J = J(indMask); K = K(indMask);
        com = round(sum([I,J,K])/numel(indMask)); % Center of mass in MATLAB coordinates;
        slices.axial.vol = volObjInit.data(:,:,com(3));
        slices.axial.roi = roiObjInit.data(:,:,com(3));
        slices.sagittal.vol = squeeze(volObjInit.data(:,com(2),:));
        slices.sagittal.roi = squeeze(roiObjInit.data(:,com(2),:));
        slices.coronal.vol = squeeze(volObjInit.data(com(1),:,:))';
        slices.coronal.roi = squeeze(roiObjInit.data(com(1),:,:))';
        fprintf('DONE')
        
        % DISPLAYING THE SLICES IN TWO PLOTS
         figure
        subplot(2,3,1), imagesc(slices.axial.vol), title('AXIAL IMAGE'), colormap gray
        subplot(2,3,2), imagesc(slices.sagittal.vol), title('SAGITTAL IMAGE'), colormap gray
        subplot(2,3,3), imagesc(slices.coronal.vol), title('CORONAL IMAGE'), colormap gray
        subplot(2,3,4), imagesc(slices.axial.roi), title('AXIAL MASK'), colormap gray
        subplot(2,3,5), imagesc(slices.sagittal.roi), title('SAGITTAL MASK'), colormap gray
        subplot(2,3,6), imagesc(slices.coronal.roi), title('CORONAL MASK'), colormap gray
        figure    
        subplot(1,3,1), imshowpair(slices.axial.vol,bwperim(slices.axial.roi,8)), axis square, title('AXIAL')
        subplot(1,3,2), imshowpair(slices.sagittal.vol,bwperim(slices.sagittal.roi,8)), axis square, title('SAGITTAL')
        subplot(1,3,3), imshowpair(slices.coronal.vol,bwperim(slices.coronal.roi,8)), axis square, title('CORONAL')
        
        % IMAGE OK?
        choiceOK = 0;
        fprintf('\n')
        while ~choiceOK
            okImage = input('--> IS THE ROI CORRECTLY DEFINED IN THE IMAGE? --> Enter 1 if yes, 0 if no: ');
            if okImage == 1 || okImage == 0
                choiceOK = 1;
            end
        end
        if ~okImage
            choiceOK = 0;
            while ~choiceOK
                okAnother = input('--> DO YOU WANT TO CHOOSE ANOTHER ROI (otherwise you still save that roi name)? --> Enter 1 if yes, 0 if no: ');
                if okAnother == 1 || okAnother == 0
                    choiceOK = 1;
                end
            end
            if ~okAnother
                roiOK = 1; % Simply to get out of the loop.
                roiCheck = [roiCheck;[name,': ERROR']];
            end
        else
            roiOK = 1;
            roiCheck = [roiCheck;[name,': OK']];
        end
        close all
    end
    
end
fprintf('-------------------------------------------------------------------------------------\n')

% Saving new roiNames.mat
cd(pathWORK)
if exist('roiNames.mat')
    info = dir('roiNames.mat');
    date = info.date; ind = strfind(date,' '); date(ind) = '_';
    newName = ['roiNames_',date,'.mat'];
    system(['mv roiNames.mat ',newName]);
end
save('roiNames','roiNames')

% Saving roiCheck.may
cd(pathWORK)
if exist('roiCheck.mat')
    info = dir('roiCheck.mat');
    date = info.date; ind = strfind(date,' '); date(ind) = '_';
    newName = ['roiCheck_',date,'.mat'];
    system(['mv roiCheck.mat ',newName]);
end
save('roiCheck','roiCheck')


% Start radiomic feature computation?
fprintf('\n')
choiceOK = 0;
while ~choiceOK
    okRadiomics = input('--> DO YOU WANT TO START THE COMPUTATION OF RADIOMICS FEATURES? --> Enter 1 if yes, 0 if no: ');
    if okRadiomics == 1 || okRadiomics == 0
        choiceOK = 1;
    end
end
if okRadiomics
    flagRadiomics = true;
else
    flagRadiomics = false;
end

cd(startpath)
end