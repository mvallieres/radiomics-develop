function [nonTextCell,textCell] = initializeFeatureNames(imageSpaceStruct)
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

nonTextCell = cell(1,3); % First entry is the names of feature types. Second entry is the name of the features for a given feature type. Third entry is the name of the extraction parameters for all features of a given feature type.
textCell = cell(1,3); % First entry is the names of feature types. Second entry is the name of the features for a given feature type. Third entry is the name of the extraction parameters for all features of a given feature type.


% NON-TEXTURE FEATURES
fieldNonText = fieldnames(imageSpaceStruct); fieldNonText = fieldNonText(1:end-1); % Removing 'texture' at the end.
nNonTextType = numel(fieldNonText); nonTextCell{1} = fieldNonText;
nonTextCell{2} = cell(nNonTextType,1);
nonTextCell{3} = cell(nNonTextType,1);
for t = 1:nNonTextType
    fieldParams = fieldnames(imageSpaceStruct.(nonTextCell{1}{t}));
    fieldFeat = fieldnames(imageSpaceStruct.(nonTextCell{1}{t}).(fieldParams{1}));
    nonTextCell{2}{t} = fieldFeat;
    nonTextCell{3}{t} = fieldParams;
end

% TEXTURE FEATURES
fieldText = fieldnames(imageSpaceStruct.texture);
nTextType = numel(fieldText); textCell{1} = fieldText;
textCell{2} = cell(nTextType,1);
textCell{3} = cell(nTextType,1);
for t = 1:nTextType
    fieldParams = fieldnames(imageSpaceStruct.texture.(textCell{1}{t}));
    fieldFeat = fieldnames(imageSpaceStruct.texture.(textCell{1}{t}).(fieldParams{1}));
    textCell{2}{t} = fieldFeat;
    textCell{3}{t} = fieldParams;
end

end