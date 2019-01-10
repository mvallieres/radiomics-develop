function patientID = getPatientID_fromRadName(radName)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% radName example: 'STS-McGill-001__T1(tumourAndEdema).MRscan';
%                                           --> patientID is 'STS-McGill-001'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

indDoubleUnder = strfind(radName,'__');
patientID = radName(1:(indDoubleUnder(1)-1));

end