function [subbands] = getWaveletSubbandsLevel_rotInv(vol,level,wavelet_name,L_filt,H_filt,pooling)

% INPUTS
% - vol: 3D array to filter
% - wavelet_name: Char array defining the 'waveletname'
% - level: level of decomposition
% - L_filt: low-pass decomposition filter (mandatory)
% - H_filt: low-pass decomposition filter (mandatory)
% - pooling: either 'max' or 'average'


% GETTING THE DIFFERENT SUBBANDS NEEDED (FLIPS)
sub_noflip = getWaveletSubbandsLevel(vol,level,wavelet_name,L_filt,H_filt);
sub_flipX = getWaveletSubbandsLevel(flip(vol,2),level,wavelet_name,L_filt,H_filt); sub_flipX = unflip(sub_flipX,2);
sub_flipY = getWaveletSubbandsLevel(flip(vol,1),level,wavelet_name,L_filt,H_filt); sub_flipY = unflip(sub_flipY,1);
sub_flipZ = getWaveletSubbandsLevel(flip(vol,3),level,wavelet_name,L_filt,H_filt); sub_flipZ = unflip(sub_flipZ,3);
sub_flipXY = getWaveletSubbandsLevel(flip(flip(vol,2),1),level,wavelet_name,L_filt,H_filt); sub_flipXY = unflip(sub_flipXY,[1,2]);
sub_flipXZ = getWaveletSubbandsLevel(flip(flip(vol,2),3),level,wavelet_name,L_filt,H_filt); sub_flipXZ = unflip(sub_flipXZ,[3,2]);
sub_flipYZ = getWaveletSubbandsLevel(flip(flip(vol,1),3),level,wavelet_name,L_filt,H_filt); sub_flipYZ = unflip(sub_flipYZ,[3,1]);
sub_flipXYZ = getWaveletSubbandsLevel(flip(flip(flip(vol,2),1),3),level,wavelet_name,L_filt,H_filt); sub_flipXYZ = unflip(sub_flipXYZ,[3,1,2]);


% INITIALIZATION
% --> For reference
% subbands.(sub_names{1}) = LLL
% subbands.(sub_names{2}) = LLH
% subbands.(sub_names{3}) = LHL
% subbands.(sub_names{4}) = LHH
% subbands.(sub_names{5}) = HLL
% subbands.(sub_names{6}) = HLH
% subbands.(sub_names{7}) = HHL
% subbands.(sub_names{8}) = HHH
subbands = sub_noflip; % INITIALIZATION --> g_0_0
sub_names = fieldnames(subbands);
szVol = size(subbands.(sub_names{1}));

% SUBBAND LLL
if strcmp(pooling,'average')
    subbands.(sub_names{1}) = subbands.(sub_names{1}) + sub_flipX.(sub_names{1});
    subbands.(sub_names{1}) = subbands.(sub_names{1}) + sub_flipXZ.(sub_names{1});
    subbands.(sub_names{1}) = subbands.(sub_names{1}) + sub_flipZ.(sub_names{1});
    subbands.(sub_names{1}) = subbands.(sub_names{1}) + sub_flipYZ.(sub_names{1});
    subbands.(sub_names{1}) = subbands.(sub_names{1}) + sub_flipY.(sub_names{1});
    subbands.(sub_names{1}) = subbands.(sub_names{1}) + sub_flipXY.(sub_names{1});
    subbands.(sub_names{1}) = subbands.(sub_names{1}) + sub_flipXYZ.(sub_names{1});
    subbands.(sub_names{1}) = subbands.(sub_names{1})/8;
elseif strcmp(pooling,'max')
    temp = [subbands.(sub_names{1})(:)';...
                    sub_flipX.(sub_names{1})(:)';...
                    sub_flipXZ.(sub_names{1})(:)';...
                    sub_flipZ.(sub_names{1})(:)';...
                    sub_flipYZ.(sub_names{1})(:)';...
                    sub_flipY.(sub_names{1})(:)';...
                    sub_flipXY.(sub_names{1})(:)';...
                    sub_flipXYZ.(sub_names{1})(:)'];
    subbands.(sub_names{1}) = nanmax(temp);
    subbands.(sub_names{1}) = reshape(subbands.(sub_names{1})',szVol(1),szVol(2),szVol(3));
end

% SUBBAND LLH
if strcmp(pooling,'average')
    subbands.(sub_names{2}) = subbands.(sub_names{2}) + sub_flipX.(sub_names{5});
    subbands.(sub_names{2}) = subbands.(sub_names{2}) + sub_flipXZ.(sub_names{2});
    subbands.(sub_names{2}) = subbands.(sub_names{2}) + sub_flipZ.(sub_names{5});
    subbands.(sub_names{2}) = subbands.(sub_names{2}) + sub_noflip.(sub_names{3});
    subbands.(sub_names{2}) = subbands.(sub_names{2}) + sub_flipYZ.(sub_names{3});
    subbands.(sub_names{2}) = subbands.(sub_names{2}) + sub_flipY.(sub_names{2});
    subbands.(sub_names{2}) = subbands.(sub_names{2}) + sub_flipXY.(sub_names{2});
    subbands.(sub_names{2}) = subbands.(sub_names{2}) + sub_flipX.(sub_names{2});
    subbands.(sub_names{2}) = subbands.(sub_names{2}) + sub_flipXY.(sub_names{5});
    subbands.(sub_names{2}) = subbands.(sub_names{2}) + sub_flipXYZ.(sub_names{5});
    subbands.(sub_names{2}) = subbands.(sub_names{2}) + sub_flipXZ.(sub_names{5});
    subbands.(sub_names{2}) = subbands.(sub_names{2}) + sub_flipXYZ.(sub_names{2});
    subbands.(sub_names{2}) = subbands.(sub_names{2}) + sub_flipYZ.(sub_names{2});
    subbands.(sub_names{2}) = subbands.(sub_names{2}) + sub_flipZ.(sub_names{2});
    subbands.(sub_names{2}) = subbands.(sub_names{2}) + sub_flipYZ.(sub_names{5});
    subbands.(sub_names{2}) = subbands.(sub_names{2}) + sub_flipY.(sub_names{5});
    subbands.(sub_names{2}) = subbands.(sub_names{2}) + sub_noflip.(sub_names{5});
    subbands.(sub_names{2}) = subbands.(sub_names{2}) + sub_flipX.(sub_names{3});
    subbands.(sub_names{2}) = subbands.(sub_names{2}) + sub_flipXZ.(sub_names{3});
    subbands.(sub_names{2}) = subbands.(sub_names{2}) + sub_flipZ.(sub_names{3});
    subbands.(sub_names{2}) = subbands.(sub_names{2}) + sub_flipXYZ.(sub_names{3});
    subbands.(sub_names{2}) = subbands.(sub_names{2}) + sub_flipXY.(sub_names{3});
    subbands.(sub_names{2}) = subbands.(sub_names{2}) + sub_flipY.(sub_names{3});
    subbands.(sub_names{2}) = subbands.(sub_names{2})/24;
elseif strcmp(pooling,'max')
    temp = [subbands.(sub_names{2})(:)';...
                    sub_flipX.(sub_names{5})(:)';...
                    sub_flipXZ.(sub_names{2})(:)';...
                    sub_flipZ.(sub_names{5})(:)';...
                    sub_noflip.(sub_names{3})(:)';...
                    sub_flipYZ.(sub_names{3})(:)';...
                    sub_flipY.(sub_names{2})(:)';...
                    sub_flipXY.(sub_names{2})(:)';...
                    sub_flipX.(sub_names{2})(:)';...
                    sub_flipXY.(sub_names{5})(:)';...
                    sub_flipXYZ.(sub_names{5})(:)';...
                    sub_flipXZ.(sub_names{5})(:)';...
                    sub_flipXYZ.(sub_names{2})(:)';...
                    sub_flipYZ.(sub_names{2})(:)';...
                    sub_flipZ.(sub_names{2})(:)';...
                    sub_flipYZ.(sub_names{5})(:)';...
                    sub_flipY.(sub_names{5})(:)';...
                    sub_noflip.(sub_names{5})(:)';...
                    sub_flipX.(sub_names{3})(:)';...
                    sub_flipXZ.(sub_names{3})(:)';...
                    sub_flipZ.(sub_names{3})(:)';...
                    sub_flipXYZ.(sub_names{3})(:)';...
                    sub_flipXY.(sub_names{3})(:)';...
                    sub_flipY.(sub_names{3})(:)'];
    subbands.(sub_names{2}) = nanmax(temp);
    subbands.(sub_names{2}) = reshape(subbands.(sub_names{2})',szVol(1),szVol(2),szVol(3));
end

% SUBBAND LHL
if strcmp(pooling,'average')
    subbands.(sub_names{3}) = subbands.(sub_names{3}) + sub_flipX.(sub_names{3});
    subbands.(sub_names{3}) = subbands.(sub_names{3}) + sub_flipXZ.(sub_names{3});
    subbands.(sub_names{3}) = subbands.(sub_names{3}) + sub_flipZ.(sub_names{3});
    subbands.(sub_names{3}) = subbands.(sub_names{3}) + sub_noflip.(sub_names{5});
    subbands.(sub_names{3}) = subbands.(sub_names{3}) + sub_flipYZ.(sub_names{5});
    subbands.(sub_names{3}) = subbands.(sub_names{3}) + sub_flipY.(sub_names{5});
    subbands.(sub_names{3}) = subbands.(sub_names{3}) + sub_flipXY.(sub_names{3});
    subbands.(sub_names{3}) = subbands.(sub_names{3}) + sub_flipX.(sub_names{5});
    subbands.(sub_names{3}) = subbands.(sub_names{3}) + sub_flipXY.(sub_names{2});
    subbands.(sub_names{3}) = subbands.(sub_names{3}) + sub_flipXYZ.(sub_names{3});
    subbands.(sub_names{3}) = subbands.(sub_names{3}) + sub_flipXZ.(sub_names{2});
    subbands.(sub_names{3}) = subbands.(sub_names{3}) + sub_flipXYZ.(sub_names{5});
    subbands.(sub_names{3}) = subbands.(sub_names{3}) + sub_flipYZ.(sub_names{3});
    subbands.(sub_names{3}) = subbands.(sub_names{3}) + sub_flipZ.(sub_names{5});
    subbands.(sub_names{3}) = subbands.(sub_names{3}) + sub_flipYZ.(sub_names{2});
    subbands.(sub_names{3}) = subbands.(sub_names{3}) + sub_flipY.(sub_names{3});
    subbands.(sub_names{3}) = subbands.(sub_names{3}) + sub_noflip.(sub_names{2});
    subbands.(sub_names{3}) = subbands.(sub_names{3}) + sub_flipX.(sub_names{2});
    subbands.(sub_names{3}) = subbands.(sub_names{3}) + sub_flipXZ.(sub_names{5});
    subbands.(sub_names{3}) = subbands.(sub_names{3}) + sub_flipZ.(sub_names{2});
    subbands.(sub_names{3}) = subbands.(sub_names{3}) + sub_flipXYZ.(sub_names{2});
    subbands.(sub_names{3}) = subbands.(sub_names{3}) + sub_flipXY.(sub_names{5});
    subbands.(sub_names{3}) = subbands.(sub_names{3}) + sub_flipY.(sub_names{2});
    subbands.(sub_names{3}) = subbands.(sub_names{3})/24;
elseif strcmp(pooling,'max')
    temp = [subbands.(sub_names{3})(:)';...
                    sub_flipX.(sub_names{3})(:)';...
                    sub_flipXZ.(sub_names{3})(:)';...
                    sub_flipZ.(sub_names{3})(:)';...
                    sub_noflip.(sub_names{5})(:)';...
                    sub_flipYZ.(sub_names{5})(:)';...
                    sub_flipY.(sub_names{5})(:)';...
                    sub_flipXY.(sub_names{3})(:)';...
                    sub_flipX.(sub_names{5})(:)';...
                    sub_flipXY.(sub_names{2})(:)';...
                    sub_flipXYZ.(sub_names{3})(:)';...
                    sub_flipXZ.(sub_names{2})(:)';...
                    sub_flipXYZ.(sub_names{5})(:)';...
                    sub_flipYZ.(sub_names{3})(:)';...
                    sub_flipZ.(sub_names{5})(:)';...
                    sub_flipYZ.(sub_names{2})(:)';...
                    sub_flipY.(sub_names{3})(:)';...
                    sub_noflip.(sub_names{2})(:)';...
                    sub_flipX.(sub_names{2})(:)';...
                    sub_flipXZ.(sub_names{5})(:)';...
                    sub_flipZ.(sub_names{2})(:)';...
                    sub_flipXYZ.(sub_names{2})(:)';...
                    sub_flipXY.(sub_names{5})(:)';...
                    sub_flipY.(sub_names{2})(:)'];
    subbands.(sub_names{3}) = nanmax(temp);
    subbands.(sub_names{3}) = reshape(subbands.(sub_names{3})',szVol(1),szVol(2),szVol(3));
end

% SUBBAND LHH
if strcmp(pooling,'average')
    subbands.(sub_names{4}) = subbands.(sub_names{4}) + sub_flipX.(sub_names{7});
    subbands.(sub_names{4}) = subbands.(sub_names{4}) + sub_flipXZ.(sub_names{4});
    subbands.(sub_names{4}) = subbands.(sub_names{4}) + sub_flipZ.(sub_names{7});
    subbands.(sub_names{4}) = subbands.(sub_names{4}) + sub_noflip.(sub_names{7});
    subbands.(sub_names{4}) = subbands.(sub_names{4}) + sub_flipYZ.(sub_names{7});
    subbands.(sub_names{4}) = subbands.(sub_names{4}) + sub_flipY.(sub_names{6});
    subbands.(sub_names{4}) = subbands.(sub_names{4}) + sub_flipXY.(sub_names{4});
    subbands.(sub_names{4}) = subbands.(sub_names{4}) + sub_flipX.(sub_names{6});
    subbands.(sub_names{4}) = subbands.(sub_names{4}) + sub_flipXY.(sub_names{6});
    subbands.(sub_names{4}) = subbands.(sub_names{4}) + sub_flipXYZ.(sub_names{7});
    subbands.(sub_names{4}) = subbands.(sub_names{4}) + sub_flipXZ.(sub_names{6});
    subbands.(sub_names{4}) = subbands.(sub_names{4}) + sub_flipXYZ.(sub_names{6});
    subbands.(sub_names{4}) = subbands.(sub_names{4}) + sub_flipYZ.(sub_names{4});
    subbands.(sub_names{4}) = subbands.(sub_names{4}) + sub_flipZ.(sub_names{6});
    subbands.(sub_names{4}) = subbands.(sub_names{4}) + sub_flipYZ.(sub_names{6});
    subbands.(sub_names{4}) = subbands.(sub_names{4}) + sub_flipY.(sub_names{7});
    subbands.(sub_names{4}) = subbands.(sub_names{4}) + sub_noflip.(sub_names{6});
    subbands.(sub_names{4}) = subbands.(sub_names{4}) + sub_flipX.(sub_names{4});
    subbands.(sub_names{4}) = subbands.(sub_names{4}) + sub_flipXZ.(sub_names{7});
    subbands.(sub_names{4}) = subbands.(sub_names{4}) + sub_flipZ.(sub_names{4});
    subbands.(sub_names{4}) = subbands.(sub_names{4}) + sub_flipXYZ.(sub_names{4});
    subbands.(sub_names{4}) = subbands.(sub_names{4}) + sub_flipXY.(sub_names{7});
    subbands.(sub_names{4}) = subbands.(sub_names{4}) + sub_flipY.(sub_names{4});
    subbands.(sub_names{4}) = subbands.(sub_names{4})/24;
elseif strcmp(pooling,'max')
    temp = [subbands.(sub_names{4})(:)';...
                    sub_flipX.(sub_names{7})(:)';...
                    sub_flipXZ.(sub_names{4})(:)';...
                    sub_flipZ.(sub_names{7})(:)';...
                    sub_noflip.(sub_names{7})(:)';...
                    sub_flipYZ.(sub_names{7})(:)';...
                    sub_flipY.(sub_names{6})(:)';...
                    sub_flipXY.(sub_names{4})(:)';...
                    sub_flipX.(sub_names{6})(:)';...
                    sub_flipXY.(sub_names{6})(:)';...
                    sub_flipXYZ.(sub_names{7})(:)';...
                    sub_flipXZ.(sub_names{6})(:)';...
                    sub_flipXYZ.(sub_names{6})(:)';...
                    sub_flipYZ.(sub_names{4})(:)';...
                    sub_flipZ.(sub_names{6})(:)';...
                    sub_flipYZ.(sub_names{6})(:)';...
                    sub_flipY.(sub_names{7})(:)';...
                    sub_noflip.(sub_names{6})(:)';...
                    sub_flipX.(sub_names{4})(:)';...
                    sub_flipXZ.(sub_names{7})(:)';...
                    sub_flipZ.(sub_names{4})(:)';...
                    sub_flipXYZ.(sub_names{4})(:)';...
                    sub_flipXY.(sub_names{7})(:)';...
                    sub_flipY.(sub_names{4})(:)'];
    subbands.(sub_names{4}) = nanmax(temp);
    subbands.(sub_names{4}) = reshape(subbands.(sub_names{4})',szVol(1),szVol(2),szVol(3));    
end

% SUBBAND HLL
if strcmp(pooling,'average')
    subbands.(sub_names{5}) = subbands.(sub_names{5}) + sub_flipX.(sub_names{2});
    subbands.(sub_names{5}) = subbands.(sub_names{5}) + sub_flipXZ.(sub_names{5});
    subbands.(sub_names{5}) = subbands.(sub_names{5}) + sub_flipZ.(sub_names{2});
    subbands.(sub_names{5}) = subbands.(sub_names{5}) + sub_noflip.(sub_names{2});
    subbands.(sub_names{5}) = subbands.(sub_names{5}) + sub_flipYZ.(sub_names{2});
    subbands.(sub_names{5}) = subbands.(sub_names{5}) + sub_flipY.(sub_names{3});
    subbands.(sub_names{5}) = subbands.(sub_names{5}) + sub_flipXY.(sub_names{5});
    subbands.(sub_names{5}) = subbands.(sub_names{5}) + sub_flipX.(sub_names{3});
    subbands.(sub_names{5}) = subbands.(sub_names{5}) + sub_flipXY.(sub_names{3});
    subbands.(sub_names{5}) = subbands.(sub_names{5}) + sub_flipXYZ.(sub_names{2});
    subbands.(sub_names{5}) = subbands.(sub_names{5}) + sub_flipXZ.(sub_names{3});
    subbands.(sub_names{5}) = subbands.(sub_names{5}) + sub_flipXYZ.(sub_names{3});
    subbands.(sub_names{5}) = subbands.(sub_names{5}) + sub_flipYZ.(sub_names{5});
    subbands.(sub_names{5}) = subbands.(sub_names{5}) + sub_flipZ.(sub_names{3});
    subbands.(sub_names{5}) = subbands.(sub_names{5}) + sub_flipYZ.(sub_names{3});
    subbands.(sub_names{5}) = subbands.(sub_names{5}) + sub_flipY.(sub_names{2});
    subbands.(sub_names{5}) = subbands.(sub_names{5}) + sub_noflip.(sub_names{3});
    subbands.(sub_names{5}) = subbands.(sub_names{5}) + sub_flipX.(sub_names{5});
    subbands.(sub_names{5}) = subbands.(sub_names{5}) + sub_flipXZ.(sub_names{2});
    subbands.(sub_names{5}) = subbands.(sub_names{5}) + sub_flipZ.(sub_names{5});
    subbands.(sub_names{5}) = subbands.(sub_names{5}) + sub_flipXYZ.(sub_names{5});
    subbands.(sub_names{5}) = subbands.(sub_names{5}) + sub_flipXY.(sub_names{2});
    subbands.(sub_names{5}) = subbands.(sub_names{5}) + sub_flipY.(sub_names{5});
    subbands.(sub_names{5}) = subbands.(sub_names{5})/24;
elseif strcmp(pooling,'max')
    temp = [subbands.(sub_names{5})(:)';...
                    sub_flipX.(sub_names{2})(:)';...
                    sub_flipXZ.(sub_names{5})(:)';...
                    sub_flipZ.(sub_names{2})(:)';...
                    sub_noflip.(sub_names{2})(:)';...
                    sub_flipYZ.(sub_names{2})(:)';...
                    sub_flipY.(sub_names{3})(:)';...
                    sub_flipXY.(sub_names{5})(:)';...
                    sub_flipX.(sub_names{3})(:)';...
                    sub_flipXY.(sub_names{3})(:)';...
                    sub_flipXYZ.(sub_names{2})(:)';...
                    sub_flipXZ.(sub_names{3})(:)';...
                    sub_flipXYZ.(sub_names{3})(:)';...
                    sub_flipYZ.(sub_names{5})(:)';...
                    sub_flipZ.(sub_names{3})(:)';...
                    sub_flipYZ.(sub_names{3})(:)';...
                    sub_flipY.(sub_names{2})(:)';...
                    sub_noflip.(sub_names{3})(:)';...
                    sub_flipX.(sub_names{5})(:)';...
                    sub_flipXZ.(sub_names{2})(:)';...
                    sub_flipZ.(sub_names{5})(:)';...
                    sub_flipXYZ.(sub_names{5})(:)';...
                    sub_flipXY.(sub_names{2})(:)';...
                    sub_flipY.(sub_names{5})(:)'];
    subbands.(sub_names{5}) = nanmax(temp);
    subbands.(sub_names{5}) = reshape(subbands.(sub_names{5})',szVol(1),szVol(2),szVol(3));
end

% SUBBAND HLH
if strcmp(pooling,'average')
    subbands.(sub_names{6}) = subbands.(sub_names{6}) + sub_flipX.(sub_names{6});
    subbands.(sub_names{6}) = subbands.(sub_names{6}) + sub_flipXZ.(sub_names{6});
    subbands.(sub_names{6}) = subbands.(sub_names{6}) + sub_flipZ.(sub_names{6});
    subbands.(sub_names{6}) = subbands.(sub_names{6}) + sub_noflip.(sub_names{4});
    subbands.(sub_names{6}) = subbands.(sub_names{6}) + sub_flipYZ.(sub_names{4});
    subbands.(sub_names{6}) = subbands.(sub_names{6}) + sub_flipY.(sub_names{4});
    subbands.(sub_names{6}) = subbands.(sub_names{6}) + sub_flipXY.(sub_names{6});
    subbands.(sub_names{6}) = subbands.(sub_names{6}) + sub_flipX.(sub_names{4});
    subbands.(sub_names{6}) = subbands.(sub_names{6}) + sub_flipXY.(sub_names{7});
    subbands.(sub_names{6}) = subbands.(sub_names{6}) + sub_flipXYZ.(sub_names{6});
    subbands.(sub_names{6}) = subbands.(sub_names{6}) + sub_flipXZ.(sub_names{7});
    subbands.(sub_names{6}) = subbands.(sub_names{6}) + sub_flipXYZ.(sub_names{4});
    subbands.(sub_names{6}) = subbands.(sub_names{6}) + sub_flipYZ.(sub_names{6});
    subbands.(sub_names{6}) = subbands.(sub_names{6}) + sub_flipZ.(sub_names{4});
    subbands.(sub_names{6}) = subbands.(sub_names{6}) + sub_flipYZ.(sub_names{7});
    subbands.(sub_names{6}) = subbands.(sub_names{6}) + sub_flipY.(sub_names{6});
    subbands.(sub_names{6}) = subbands.(sub_names{6}) + sub_noflip.(sub_names{7});
    subbands.(sub_names{6}) = subbands.(sub_names{6}) + sub_flipX.(sub_names{7});
    subbands.(sub_names{6}) = subbands.(sub_names{6}) + sub_flipXZ.(sub_names{4});
    subbands.(sub_names{6}) = subbands.(sub_names{6}) + sub_flipZ.(sub_names{7});
    subbands.(sub_names{6}) = subbands.(sub_names{6}) + sub_flipXYZ.(sub_names{7});
    subbands.(sub_names{6}) = subbands.(sub_names{6}) + sub_flipXY.(sub_names{4});
    subbands.(sub_names{6}) = subbands.(sub_names{6}) + sub_flipY.(sub_names{7});
    subbands.(sub_names{6}) = subbands.(sub_names{6})/24;
elseif strcmp(pooling,'max')
    temp = [subbands.(sub_names{6})(:)';...
                    sub_flipX.(sub_names{6})(:)';...
                    sub_flipXZ.(sub_names{6})(:)';...
                    sub_flipZ.(sub_names{6})(:)';...
                    sub_noflip.(sub_names{4})(:)';...
                    sub_flipYZ.(sub_names{4})(:)';...
                    sub_flipY.(sub_names{4})(:)';...
                    sub_flipXY.(sub_names{6})(:)';...
                    sub_flipX.(sub_names{4})(:)';...
                    sub_flipXY.(sub_names{7})(:)';...
                    sub_flipXYZ.(sub_names{6})(:)';...
                    sub_flipXZ.(sub_names{7})(:)';...
                    sub_flipXYZ.(sub_names{4})(:)';...
                    sub_flipYZ.(sub_names{6})(:)';...
                    sub_flipZ.(sub_names{4})(:)';...
                    sub_flipYZ.(sub_names{7})(:)';...
                    sub_flipY.(sub_names{6})(:)';...
                    sub_noflip.(sub_names{7})(:)';...
                    sub_flipX.(sub_names{7})(:)';...
                    sub_flipXZ.(sub_names{4})(:)';...
                    sub_flipZ.(sub_names{7})(:)';...
                    sub_flipXYZ.(sub_names{7})(:)';...
                    sub_flipXY.(sub_names{4})(:)';...
                    sub_flipY.(sub_names{7})(:)'];
    subbands.(sub_names{6}) = nanmax(temp);
    subbands.(sub_names{6}) = reshape(subbands.(sub_names{6})',szVol(1),szVol(2),szVol(3));
end

% SUBBAND HHL
if strcmp(pooling,'average')
    subbands.(sub_names{7}) = subbands.(sub_names{7}) + sub_flipX.(sub_names{4});
    subbands.(sub_names{7}) = subbands.(sub_names{7}) + sub_flipXZ.(sub_names{7});
    subbands.(sub_names{7}) = subbands.(sub_names{7}) + sub_flipZ.(sub_names{4});
    subbands.(sub_names{7}) = subbands.(sub_names{7}) + sub_noflip.(sub_names{6});
    subbands.(sub_names{7}) = subbands.(sub_names{7}) + sub_flipYZ.(sub_names{6});
    subbands.(sub_names{7}) = subbands.(sub_names{7}) + sub_flipY.(sub_names{7});
    subbands.(sub_names{7}) = subbands.(sub_names{7}) + sub_flipXY.(sub_names{7});
    subbands.(sub_names{7}) = subbands.(sub_names{7}) + sub_flipX.(sub_names{7});
    subbands.(sub_names{7}) = subbands.(sub_names{7}) + sub_flipXY.(sub_names{4});
    subbands.(sub_names{7}) = subbands.(sub_names{7}) + sub_flipXYZ.(sub_names{4});
    subbands.(sub_names{7}) = subbands.(sub_names{7}) + sub_flipXZ.(sub_names{4});
    subbands.(sub_names{7}) = subbands.(sub_names{7}) + sub_flipXYZ.(sub_names{7});
    subbands.(sub_names{7}) = subbands.(sub_names{7}) + sub_flipYZ.(sub_names{7});
    subbands.(sub_names{7}) = subbands.(sub_names{7}) + sub_flipZ.(sub_names{7});
    subbands.(sub_names{7}) = subbands.(sub_names{7}) + sub_flipYZ.(sub_names{4});
    subbands.(sub_names{7}) = subbands.(sub_names{7}) + sub_flipY.(sub_names{4});
    subbands.(sub_names{7}) = subbands.(sub_names{7}) + sub_noflip.(sub_names{4});
    subbands.(sub_names{7}) = subbands.(sub_names{7}) + sub_flipX.(sub_names{6});
    subbands.(sub_names{7}) = subbands.(sub_names{7}) + sub_flipXZ.(sub_names{6});
    subbands.(sub_names{7}) = subbands.(sub_names{7}) + sub_flipZ.(sub_names{6});
    subbands.(sub_names{7}) = subbands.(sub_names{7}) + sub_flipXYZ.(sub_names{6});
    subbands.(sub_names{7}) = subbands.(sub_names{7}) + sub_flipXY.(sub_names{6});
    subbands.(sub_names{7}) = subbands.(sub_names{7}) + sub_flipY.(sub_names{6});
    subbands.(sub_names{7}) = subbands.(sub_names{7})/24;
elseif strcmp(pooling,'max')
    temp = [subbands.(sub_names{7})(:)';...
                    sub_flipX.(sub_names{4})(:)';...
                    sub_flipXZ.(sub_names{7})(:)';...
                    sub_flipZ.(sub_names{4})(:)';...
                    sub_noflip.(sub_names{6})(:)';...
                    sub_flipYZ.(sub_names{6})(:)';...
                    sub_flipY.(sub_names{7})(:)';...
                    sub_flipXY.(sub_names{7})(:)';...
                    sub_flipX.(sub_names{7})(:)';...
                    sub_flipXY.(sub_names{4})(:)';...
                    sub_flipXYZ.(sub_names{4})(:)';...
                    sub_flipXZ.(sub_names{4})(:)';...
                    sub_flipXYZ.(sub_names{7})(:)';...
                    sub_flipYZ.(sub_names{7})(:)';...
                    sub_flipZ.(sub_names{7})(:)';...
                    sub_flipYZ.(sub_names{4})(:)';...
                    sub_flipY.(sub_names{4})(:)';...
                    sub_noflip.(sub_names{4})(:)';...
                    sub_flipX.(sub_names{6})(:)';...
                    sub_flipXZ.(sub_names{6})(:)';...
                    sub_flipZ.(sub_names{6})(:)';...
                    sub_flipXYZ.(sub_names{6})(:)';...
                    sub_flipXY.(sub_names{6})(:)';...
                    sub_flipY.(sub_names{6})(:)'];
    subbands.(sub_names{7}) = nanmax(temp);
    subbands.(sub_names{7}) = reshape(subbands.(sub_names{7})',szVol(1),szVol(2),szVol(3));
end

% SUBBAND HHH
if strcmp(pooling,'average')
    subbands.(sub_names{8}) = subbands.(sub_names{8}) + sub_flipX.(sub_names{8});
    subbands.(sub_names{8}) = subbands.(sub_names{8}) + sub_flipXZ.(sub_names{8});
    subbands.(sub_names{8}) = subbands.(sub_names{8}) + sub_flipZ.(sub_names{8});
    subbands.(sub_names{8}) = subbands.(sub_names{8}) + sub_flipYZ.(sub_names{8});
    subbands.(sub_names{8}) = subbands.(sub_names{8}) + sub_flipY.(sub_names{8});
    subbands.(sub_names{8}) = subbands.(sub_names{8}) + sub_flipXY.(sub_names{8});
    subbands.(sub_names{8}) = subbands.(sub_names{8}) + sub_flipXYZ.(sub_names{8});
    subbands.(sub_names{8}) = subbands.(sub_names{8})/8;
elseif strcmp(pooling,'max')
    temp = [subbands.(sub_names{8})(:)';...
                    sub_flipX.(sub_names{8})(:)';...
                    sub_flipXZ.(sub_names{8})(:)';...
                    sub_flipZ.(sub_names{8})(:)';...
                    sub_flipYZ.(sub_names{8})(:)';...
                    sub_flipY.(sub_names{8})(:)';...
                    sub_flipXY.(sub_names{8})(:)';...
                    sub_flipXYZ.(sub_names{8})(:)'];
    subbands.(sub_names{8}) = nanmax(temp);
    subbands.(sub_names{8}) = reshape(subbands.(sub_names{8})',szVol(1),szVol(2),szVol(3));        
end

end



function subbands = unflip(subbands,dims)
sub_names = fieldnames(subbands); nSub = numel(sub_names);
nDim = numel(dims);
for s = 1:nSub
    sub_name = sub_names{s};
    for d = 1:nDim
        dim = dims(d);
        subbands.(sub_name) = flip(subbands.(sub_name),dim);
    end
end

end