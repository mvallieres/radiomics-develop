function radiomicsTable = createRadiomicsTable(radiomicsFilePaths,imageSpace);
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


% INITIALIZATIONS OF RADIOMICS STRUCTURES
nFiles = numel(radiomicsFilePaths);
patientID = cell(nFiles,1);
radStructs = cell(nFiles,1); fileOpen = false(nFiles,1);
for f = 1:nFiles
    try
        radStruct = load(radiomicsFilePaths{f}); radStruct = struct2cell(radStruct);  radStruct = radStruct{1};
        radStructs{f} = radStruct;
        fileOpen(f) = true;
    end
    [~,nameFile,~] = fileparts(radiomicsFilePaths{f}); 
    patientID{f} = getPatientID_fromRadName(nameFile);
end



% INITIALIZE FEATURE NAMES
patientRandom = ceil(nFiles * rand(1,1));
radiomicsStruct = load(radiomicsFilePaths{patientRandom}); radiomicsStruct = struct2cell(radiomicsStruct); radiomicsStruct = radiomicsStruct{1};
imageSpaceStruct = radiomicsStruct.(imageSpace); %IMAGE SPACE STRUCTURE --> .morph, .locInt, ...,  .texture
[nonTextCell,textCell] = initializeFeatureNames(imageSpaceStruct);




% CREATE TABLE DATA
radiomicsTable = table; strTable = []; strNames = '||'; countVar = 0;

% Non-texture features
for type = 1:numel(nonTextCell{1})
    for param = 1:numel(nonTextCell{3}{type})
        for feat = 1:numel(nonTextCell{2}{type})
            countVar = countVar + 1; nameF = ['radVar',num2str(countVar)];
            eval([nameF,' = zeros(nFiles,1);'])
            realNameF = [nonTextCell{1}{type},'__',nonTextCell{2}{type}{feat},'__',nonTextCell{3}{type}{param}];
            strTable = [strTable,nameF,',']; strNames = [strNames,nameF,':',realNameF,'||'];
            for f = 1:nFiles
                if fileOpen(f)
                    try
                        val = radStructs{f}.(imageSpace).(nonTextCell{1}{type}).(nonTextCell{3}{type}{param}).(nonTextCell{2}{type}{feat});
                    catch
                        val = NaN;
                    end
                    if ischar(val) || isempty(val)
                        val = NaN;
                    end
                else
                    val = NaN;
                end
                eval([nameF,'(f) = val;'])
            end
        end
    end
end

% Texture features
for type = 1:numel(textCell{1})
    for param = 1:numel(textCell{3}{type})
        for feat = 1:numel(textCell{2}{type})
            countVar = countVar + 1; nameF = ['radVar',num2str(countVar)];
            eval([nameF,' = zeros(nFiles,1);'])
            realNameF = [textCell{1}{type},'__',textCell{2}{type}{feat},'__',textCell{3}{type}{param}];
            strTable = [strTable,nameF,',']; strNames = [strNames,nameF,':',realNameF,'||'];
            for f = 1:nFiles
                if fileOpen(f)
                    try
                        val = radStructs{f}.(imageSpace).texture.(textCell{1}{type}).(textCell{3}{type}{param}).(textCell{2}{type}{feat});
                    catch
                        val = NaN;
                    end
                    if ischar(val) || isempty(val)
                        val = NaN;
                    end
                else
                    val = NaN;
                end
                eval([nameF,'(f) = val;'])
            end
        end
    end
end
strTable(end) = []; % Removing the last comma
eval(['radiomicsTable = table(',strTable,');'])
radiomicsTable.Properties.UserData = strNames;
radiomicsTable.Properties.RowNames = patientID;
radiomicsTable.Properties.DimensionNames{1} = 'PatientID';

end