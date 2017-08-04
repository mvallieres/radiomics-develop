function [subbands] = getWaveletSubbands(vol)

% MAKING SURE THE VOLUME HAS EVEN SIZE (necessary for swt2)
% Adding a layer identical to the last one of the volume. This should not
% create problems for sufficiently large bounding boxes. Is box10 ok?

% Getting all sub-bands

% IF VOLUME WAS UNEVEN, REMOVE VOXELS

end