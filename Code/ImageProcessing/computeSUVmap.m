function [SUVmap] = computeSUVmap(rawPET,dicomH)
% -------------------------------------------------------------------------
% function [SUVmap] = computeSUVmap(rawPET,dicomH)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function computes the SUVmap of a raw input PET volume. It is 
% assumed that the calibration factor was applied beforehand to the PET 
% volume (e.g., rawPET = rawPET*RescaleSlope + RescaleIntercept).
% -------------------------------------------------------------------------
% INPUTS:
% - rawPET: 3D array representing the PET volume in raw format.
% - dicomH: DICOM header of one of the corresponding slice of 'rawPET'.
% -------------------------------------------------------------------------
% OUTPUTS:
% - SUVmap: 'rawPET' converted to SUVs (standard uptake values).
% -------------------------------------------------------------------------
% AUTHOR(S): 
% - Martin Vallieres <mart.vallieres@gmail.com>
% - Issam El Naqa <ielnaqa@med.umich.edu>
% - CERR development team <http://www.cerr.info/>
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
% --> Copyright (C) 2015  Martin Vallieres, Issam El Naqa
% --> Copyright 2010, Joseph O. Deasy, on behalf of the CERR development team
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

% Get patient weight
if isfield(dicomH,'PatientWeight')
    weight = dicomH.PatientWeight*1000;  % in grams
elseif isfield(dicomH,'PatientsWeight')
    weight = dicomH.PatientsWeight*1000;  % in grams
else
    weight = [];
end
if isempty(weight) || weight == 0
    weight = 75000; % Estimation
end

try
    % Get Scan time
    scantime = dcm_hhmmss(dicomH.AcquisitionTime);
    % Start Time for the Radiopharmaceutical Injection
    injection_time = dcm_hhmmss(dicomH.RadiopharmaceuticalInformationSequence.Item_1.RadiopharmaceuticalStartTime);
    % Half Life for Radionuclide
    half_life = dicomH.RadiopharmaceuticalInformationSequence.Item_1.RadionuclideHalfLife;
    % Total dose injected for Radionuclide
    injected_dose = dicomH.RadiopharmaceuticalInformationSequence.Item_1.RadionuclideTotalDose;

    % Calculate decay
    decay = exp(-log(2)*(scantime-injection_time)/half_life);
    % Calculate the dose decayed during procedure
    injected_dose_decay = injected_dose*decay; % in Bq

catch % Estimation
    decay = exp(-log(2)*(1.75*3600)/6588); % 90 min waiting time, 15 min preparation
    injected_dose_decay = 420000000 * decay; % 420 MBq
end

% Calculate SUV
SUVmap = rawPET*weight/injected_dose_decay;

end

% CERR UTILITY FUNCTION (can be found at: https://github.com/adityaapte/CERR)
function [totSec] = dcm_hhmmss(dateStr)
if ~ischar(dateStr)
    dateStr = num2str(dateStr);
end
hh = str2double(dateStr(1,1:2));
mm = str2double(dateStr(1,3:4));
ss = str2double(dateStr(1,5:6));
totSec = hh*60*60 + mm*60 + ss;
end