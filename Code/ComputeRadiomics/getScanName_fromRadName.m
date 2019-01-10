function scanName = getScanName_fromRadName(radName)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% radName example: 'STS-McGill-001__T1(tumourAndEdema).MRscan';
%                                           --> scanName  is 'T1'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

indDoubleUnder = strfind(radName,'__');
indOpenPar = strfind(radName,'(');
scanName = radName((indDoubleUnder(1)+2):(indOpenPar(end)-1));

end