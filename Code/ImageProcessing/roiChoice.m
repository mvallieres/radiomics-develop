function roiNames = roiChoice(pathWORK,pathDATA)
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

% List of MATLAB imaging data files
cd(pathDATA), list = dir('*.mat'); nFile = numel(list);


% MANUAL CHOICE OF ROI
fprintf('\n\n%%%%%%%%%%%% MANUAL CHOICE OF ROIs %%%%%%%%%%%%')
for f = 1:nFile
    name = list(f).name;
    [ID,scan,type] = findIDscanType(name); roiNames{f,1} = ID; roiNames{f,2} = scan; roiNames{f,3} = type;
    sData = load(name); sData = struct2cell(sData); sData = sData{1};
    nROI = numel(sData{2}.scan.contour);
    fprintf('\n----- LIST OF ROIs FOR "%s" -----',name);
    for r = 1:nROI
        fprintf(['\n',num2str(r),'. ',sData{2}.scan.contour(r).name,' -- ',sData{2}.scan.contour(r).nameSet])
    end
    ok = 0;
    fprintf('\n')
    while ~ok
        roi = input('--> WHICH ROI DO YOU CHOOSE FOR THAT SCAN? Example --> Enter 1 for roi #1, [2,4] for roi #2 and #4: ');
        if isnumeric(roi) && ~sum((roi - nROI) > 0) && ~sum((1 - roi) > 0)
            ok = 1;
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
end
fprintf('-------------------------------------------------------------------------------------')


% Saving new roiNames.mat
cd(pathWORK)
if exist('roiNames.mat')
    info = dir('roiNames.mat');
    date = info.date; ind = strfind(date,' '); date(ind) = '_';
    newName = ['roiNames_',date,'.mat'];
    system(['mv roiNames.mat ',newName]);
end
save('roiNames','roiNames')

cd(startpath)
end