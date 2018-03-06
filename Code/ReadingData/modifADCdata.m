function modifADCdata(pathRead)
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

cd(pathRead), listADC = dir('*ADC.MRscan.mat'); nADC = numel(listADC);
for a = 1:nADC
    name = listADC(a).name;
    sData = load(name); sData = struct2cell(sData); sData = sData{1};
    sData{2}.type = 'ADCscan'; save(name,'sData','-v7.3')
    indDot = strfind(name,'.');
    system(['mv ',name,' ',name(1:indDot(1)),'ADCscan.mat']);
end

cd(startpath)
end