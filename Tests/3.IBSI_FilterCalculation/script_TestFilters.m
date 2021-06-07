

% IMPORTANT.
% --> Orientation between Alex nifti phantom (IBSI 1) and nifitread is not
% the same. We need to transpose every image first.
% STEPS:
% 1. niftiread
% 2. Transpose every image
% 3. X: left to right, Y: top to bottom, Z: slice dimension
% 4. Do your stuff
% 5. Transpose again before saving.

% This is implemented, let's confirm with the team.
clear, clc
tic


%%%%%%%%%%%%%%%%%%%%%%%%
%      OPTIONS         %
%%%%%%%%%%%%%%%%%%%%%%%%
pathWORK = pwd;
pathMaps = fullfile(pathWORK,'ResponseMaps');
pathPhantoms = fullfile(pathWORK,'Phantoms');
%addpath(genpath(fullfile(pathWORK,'Functions_Filters')));

% ---> Options for mean filters. This defines the experiments to perform.
% ---> Make sure the number of entries correspond
filtersInfo.Mean.kernelSize = [15,15,15,15];
filtersInfo.Mean.phantom = {'checkerboard.nii.gz','checkerboard.nii.gz','checkerboard.nii.gz','checkerboard.nii.gz'};
filtersInfo.Mean.boundary = {0,'replicate','circular','symmetric'};
filtersInfo.Mean.infos = {'Phase1_1a1','Phase1_1a2','Phase1_1a3','Phase1_1a4'};

% ---> Options for LoG filters. This defines the experiments to perform.
% ---> Make sure the number of entries correspond
filtersInfo.LoG.kernelSize = [13,21]; % In intrinsic image units. 2*4sigma + 1 (spatial to mm). 
filtersInfo.LoG.sigma = [3,5]; % In mm
filtersInfo.LoG.phantom = {'response.nii.gz','checkerboard.nii.gz'}; % Phantoms to analyze
filtersInfo.LoG.boundary = {0,'symmetric'}; % Either a numeric number (constant padding), 'replicate' (nearest value padding), 'circular' (periodization padding) or 'symmetric' (mirror padding)
filtersInfo.LoG.infos = {'Phase1_2a','Phase1_2b'};

% ---> Options for Laws filters. This defines the experiments to perform.
% ---> Make sure the number of entries in phantom, boundary, kernels,
%      pooling, energy and deltas correspond.
E5 = 1/sqrt(10)*[-1,-2,0,2,1];
L5 = 1/sqrt(70)*[1,4,6,4,1];
S5 = 1/sqrt(6)*[-1,0,2,0,-1];
W5 = 1/sqrt(10)*[-1,2,0,-2,1];
R5 = 1/sqrt(70)*[1,-4,6,-4,1];
single_rot_E5L5S5 = {E5,L5,S5};
single_rot_E5W5R5 = {E5,W5,R5};
rot_invariant_E5L5S5 = set_rotInv_3D(single_rot_E5L5S5);
rot_invariant_E5W5R5 = set_rotInv_3D(single_rot_E5W5R5);
filtersInfo.Laws.phantom = {'response.nii.gz','response.nii.gz','response.nii.gz','checkerboard.nii.gz','checkerboard.nii.gz','checkerboard.nii.gz'}; 
filtersInfo.Laws.boundary = {0,0,0,'symmetric','symmetric','symmetric'};
filtersInfo.Laws.kernels = {single_rot_E5L5S5,rot_invariant_E5L5S5,rot_invariant_E5L5S5,single_rot_E5W5R5,rot_invariant_E5W5R5,rot_invariant_E5W5R5}; % We will always use rotation invariance. Add another cell for other types of kernels.
filtersInfo.Laws.pooling = {'average','max','max','average','max','max'};
filtersInfo.Laws.energy = [false,false,true,false,false,true];
filtersInfo.Laws.deltas = [NaN,NaN,7,NaN,NaN,7]; % In units of voxels
filtersInfo.Laws.infos = {'Phase1_3a1','Phase1_3a2','Phase1_3a3','Phase1_3b1','Phase1_3b2','Phase1_3b3'};
 
% % ---> Options for Gabor filters. This defines the experiments to perform.
% % ---> Make sure the number of entries correspond.
filtersInfo.Gabor.phantom = {'response.nii.gz','response.nii.gz','sphere.nii.gz','sphere.nii.gz'};
filtersInfo.Gabor.boundary = {0,0,'symmetric','symmetric'}; % Use SFB to calculate an efficient padding size?
filtersInfo.Gabor.sigma = [10,10,20,20]; % In units of mm
filtersInfo.Gabor.lambda = [4,4,8,8]; % In units of mm
filtersInfo.Gabor.gamma = [1/2,1/2,5/2,5/2];
filtersInfo.Gabor.theta = {-60,-(45:45:360),-225,-((45/2):(45/2):360)}; % Minus signs: difference between Matlab and IBSI definition
filtersInfo.Gabor.pooling = {'average','average','average','average'};
filtersInfo.Gabor.threeD = [false,true,false,true]; % If false, filters are only applied using the in-plane orientation, no 3-way pass average through orthogonal planes.
filtersInfo.Gabor.infos = {'Phase1_4a1','Phase1_4a2','Phase1_4b1','Phase1_4b2'};

% ---> Options for Wavelet filters. This defines the experiments to perform.
% ---> Make sure the number of entries correspond
filtersInfo.Wavelet.phantom = {'response.nii.gz','response.nii.gz','sphere.nii.gz','sphere.nii.gz','checkerboard.nii.gz','checkerboard.nii.gz'};
filtersInfo.Wavelet.boundary = {0,0,'circular','circular','symmetric','symmetric'};
filtersInfo.Wavelet.basisFunction = {'db2','db2','coif1','coif1','haar','haar'};
filtersInfo.Wavelet.direction = {'LHL','LHL','HHL','HHL','LLL','HHH'};
filtersInfo.Wavelet.level = [1,1,1,1,2,2];
filtersInfo.Wavelet.rotInv = [false,true,false,true,true,true];
filtersInfo.Wavelet.pooling = {'','average','','average','average','average'};
filtersInfo.Wavelet.infos = {'Phase1_5a1','Phase1_5a2','Phase1_6a1','Phase1_6a2','Phase1_7a1','Phase1_7a2'};

% % ---> Options for Riesz filters. This defines the experiments to perform.
% % ---> Make sure the number of entries correspond.
% ----------------------------------------------------






%%%%%%%%%%%%%%%%%%%%%%%%
%     MAIN CODE        %
%%%%%%%%%%%%%%%%%%%%%%%%
% --> It is assumed that any imagingVolume has isoptropic voxel sizes.

filters = fieldnames(filtersInfo); nFilters = numel(filters);
for f = 1:nFilters
    filter = filters{f};
    filterInfo = filtersInfo.(filter);
    phantoms = filterInfo.phantom; nPhantoms = numel(phantoms);
    boundaries = filterInfo.boundary; nBoundaries = numel(boundaries);
    
    for p = 1:nPhantoms
        phantom = phantoms{p};
        vol = single(niftiread(fullfile(pathPhantoms,phantom)));  % CHANGING THE VOLUME TO FLOAT 32 --> Must change the header accordingly before writing the new nifiti file.
        vol = transpose3D(vol);
        info = niftiinfo(fullfile(pathPhantoms,phantom)); voxDim = info.PixelDimensions;
        info = adaptNiiHeader(info,'single'); % CHANGING THE VOLUME TO FLOAT 32 --> Must change the header accordingly before writing the new nifiti file.
        phantom = getPhantomName(phantom);
        boundary = boundaries{p};

        % Mean filter specifics
        if strcmp(filter,'Mean')
            kernelSize = filterInfo.kernelSize(p);
            volFilter = single(imfilter(vol,ones(kernelSize,kernelSize,kernelSize),boundary)/(kernelSize^3));
            nameSave = filterInfo.infos{p};
            niftiwrite(transpose3D(volFilter),fullfile(pathMaps,nameSave),info,'Compressed',true);
        end

        % LoG filter specifics
        if strcmp(filter,'LoG')
            kernelSize = filterInfo.kernelSize(p);
            sigma = filterInfo.sigma(p);
            nameSave = filterInfo.infos{p};
            sigma = sigma/voxDim(1); % In intrinsic image units, also assuming that the imagingVolume has isotropic voxel size
            kernel = getLoG_3Dkernel(kernelSize,sigma);
            getAndSaveNiiMap(transpose3D(vol),kernel,boundary,info,fullfile(pathMaps,nameSave));
        end


        % Laws filters
        if strcmp(filter,'Laws')
            kernel = filterInfo.kernels{p}; 
            pooling = filterInfo.pooling{p};
            energy = filterInfo.energy(p);
            delta = filterInfo.deltas(p); 
            volFilter = single(filter_separable(vol,kernel,boundary,pooling,energy,delta));
            nameSave = filterInfo.infos{p};
            niftiwrite(transpose3D(volFilter),fullfile(pathMaps,nameSave),info,'Compressed',true);
        end
        
        % Gabor filters
        if strcmp(filter,'Gabor')
            
            % Initialization
            sigma = filterInfo.sigma(p)/voxDim(1); % In intrinsic image units
            lambda = filterInfo.lambda(p)/voxDim(1); % In intrinsic image units
            gamma = filterInfo.gamma(p); % VERIFY THAT THIS IS DIRECTLY THE RATIO OF SEMI-MAJOR TO SEMI-MINOR AXES
            theta = filterInfo.theta{p}; 
            pooling = filterInfo.pooling{p}; 
            threeD = filterInfo.threeD(p);
            nameSave = filterInfo.infos{p};
            
            % Finding other relevant parameters 
            F_b = log2((sigma/lambda*pi + sqrt(log(2)/2)) / (sigma/lambda*pi - sqrt(log(2)/2)));
            padSize = 2*ceil(sigma*7)+1; % Pad size in all directions need to be set to the size of the filter (spatial support) --> see gabor.m
            padSize = [padSize,padSize,padSize]; % In 3D.
            
            % Computation
            vol = padVol(vol,padSize,boundary);
            vol = apply_gaborfilt(vol,lambda,theta,F_b,gamma,pooling,threeD); % New function below
            vol = unPadVol(vol,padSize);
            
            % Saving the response map
            niftiwrite(transpose3D(vol),fullfile(pathMaps,nameSave),info,'Compressed',true);
        end

        % Wavelet filters
        if strcmp(filter,'Wavelet')
                
            % Initialization
            waveletName = filterInfo.basisFunction{p};
            nameSave = filterInfo.infos{p};
            direction = filterInfo.direction{p};
            level = filtersInfo.Wavelet.level(p);
            rotInv = filtersInfo.Wavelet.rotInv(p);
            pooling = filtersInfo.Wavelet.pooling{p};

            % Computation
            padSize = size(vol)/4; % This should be sufficient
            vol = padVol(vol,padSize,boundary);
            [L_filt,H_filt,~,~] =  wfilters(waveletName);
            if ~rotInv
                [subbands] = getWaveletSubbandsLevel(vol,level,waveletName,L_filt,H_filt); % Computation is much faster when the filters are provided, why?
            else
                [subbands] = getWaveletSubbandsLevel_rotInv(vol,level,waveletName,L_filt,H_filt,pooling);
            end
            field_subbands = fieldnames(subbands);
            nSubbands = numel(field_subbands);
            for s = 1:nSubbands
                field_subband = field_subbands{s};
                subbands.(field_subband) = unPadVol(subbands.(field_subband),padSize);
            end
            wavNameSave = replaceCharacter(waveletName,'.','dot'); % Using 'getWaveletSubbandsLevel' may change the name of the wavelet (if there is a dot). 
            volFilter = single(subbands.([direction,'_',wavNameSave]));
            niftiwrite(transpose3D(volFilter),fullfile(pathMaps,nameSave),info,'Compressed',true);
        end
    end
    
end
toc
% ----------------------------------------------------






%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    IN-LINE FUNCTIONS    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --> This is temporary. Proper functions and/or classes will be created later.

%%%%%%%%%% FILTERING FUNCTIONS %%%%%%%%%%

% 1. getLoG_3Dkernel
function kernel = getLoG_3Dkernel(kernelSize,sigma)
    % --> Here, both "kernelSize" and "sigma" are in intrinsic image units.
    % --> Kernel size and sigma  are expected to be single numbers specifying the size of
    %     the kernel/std in all 3 dimensions (isotropic voxel size is used in radiomics)

    if mod(kernelSize,1) 
        error('The kernel size must be an integer')
    end

    % Creation of meshgrid
    siz = (kernelSize-1)/2;
    [X,Y,Z] = meshgrid(-siz:siz,-siz:siz,-siz:siz);
    normSquare = X.*X + Y.*Y + Z.*Z;

    % Calculating the Gaussian
    gaussian = exp(-normSquare/(2*sigma^2));
    gaussian(gaussian<eps*max(gaussian(:))) = 0;


    % --> LOG CALCULATION
    % OPTION 1: Calculate the Laplacian (according to current IBSI equation)
    %kernel = (-1/(pi*sigma^2)) * (1 - normSquare/(2*sigma^2)) .* gaussian; 
    
    % OPTION 2: USING (transformed) MATLAB 2018b EQUATION
    %sumGaussian = sum(gaussian(:));
    %if sumGaussian ~= 0
    %   gaussian = gaussian/sumGaussian;
    %end
    %kernel = (-1/(sigma^2)) * (3 - normSquare/(sigma^2)) .* gaussian;
    
    % --> OPTION 2 AND 3 GIVES THE SAME RESULTS AS MATLAB 2018b with fspecial3, i.e. same
    % as Phil.
    
    % OPTION 3: USING UPDATED EQUATION WITH UPDATED NORMALIZATION FACTOR
    Nd = 3;
    normFactor = (1/(sqrt(2*pi) * sigma))^Nd;
    kernel = normFactor * (-1/(sigma^2)) * (Nd - normSquare/(sigma^2)) .* gaussian;

    % Final step: Make the filter sum to 0 (this may not be necessary for a Laplacian though)
    kernel  = kernel - sum(kernel(:))/(kernelSize^3);

end

function vol_filt = apply_gaborfilt(vol,lambda,theta,F_b,gamma,pooling,threeD)

% INITIALIZATIONS
if threeD
    % Rotations that would be required to go in "axial", "sagittal" and "coronal" planes
    perms = {[1,2,3],[1,3,2],[3,2,1]};
    rots = [0,90,180];
else % Staying in the "axial" plane
    perms = {[1,2,3]};
    rots = [0]; 
end
nPerms = numel(perms);
nTheta = numel(theta);

% COMPUTATIONS
vol_filt = zeros(size(vol));
for p = 1:nPerms
    perm = perms{p};
    rot = rots(p);
    
    % Rotating the volume
    vol_rot = imrotate(permute(vol,perm),rot);
    szRot = size(vol_rot); num_in_slice = szRot(1) * szRot(2);
    
    % Applying the Gabor filter for all slices
    nSlices = szRot(3);
    for s = 1:nSlices
        
        % Initialization
        im = vol_rot(:,:,s);
        stack = zeros(nTheta,num_in_slice);
        
        % Computation of gabor filter
        for t = 1:nTheta
            orientation = theta(t);
            im_filt = imgaborfilt(im,lambda,orientation,'SpatialFrequencyBandwidth',F_b,'SpatialAspectRatio',gamma);
            stack(t,:) = im_filt(:)';
        end
        
        % Pooling
        if strcmp(pooling,'average')
            im_filt = nanmean(stack,1);
        elseif strcmp(pooling,'max')
            stack = [stack;ones(1,size(stack,2))*-Inf]; %  In case there is only one row in stack, MATLAB specific.
            im_filt = nanmax(stack);
        end
        im_filt = reshape(im_filt',szRot(1),szRot(2));
        
        % Slice done
        vol_rot(:,:,s) = im_filt;
    end
    
    % Unrotating the volume
    vol_rot = ipermute(imrotate(vol_rot,-rot),perm);
    
    % Summing up the 2D response of the rotate volume (to average in the
    % three orthogonal planes, just outside the loop)
    vol_filt = vol_filt + vol_rot;
    
end
vol_filt = vol_filt/nPerms; % Averaging over the three orthogonal planes

end

% 2. filter_separable
function vol_filt = filter_separable(vol,kernel,boundary,pooling,energy,delta)
    
    % Pad size
    padSize = numel(kernel{1,1})*2; % This should be sufficient.
    padSize = [padSize,padSize,padSize];
    
    % Getting all filtered volumes from all defined rotations.
    nRot = size(kernel,1);
    cell_rot = cell(1,nRot);
    for r = 1:nRot
        temp = vol;
            
        % Padding the volume
        cell_rot{r} = padVol(temp,padSize,boundary);
        
        % X direction: left to right
        cell_rot{r} = imfilter(cell_rot{r},kernel{r,1},'conv');
        
        % Y direction: top to bottom
        cell_rot{r} = imfilter(cell_rot{r},kernel{r,2}','conv'); % Transposing the kernel to apply in the y-direction 
        
        % Z direction: slice dimension
        cell_rot{r} = permute(cell_rot{r},[1,3,2]); % Permuting around the y-axis
        cell_rot{r} = imfilter(cell_rot{r},kernel{r,3},'conv');
        cell_rot{r} = permute(cell_rot{r},[1,3,2]); % Permuting back around the y-axis
        
        % Unpadding the volume
        cell_rot{r} = unPadVol(cell_rot{r},padSize);
        
    end
    szVol = size(cell_rot{1});
    
    % Pooling results
    if nRot == 1 && strcmp(pooling,'average')
        vol_filt = cell_rot{1};
    else
        nVox = numel(cell_rot{1});
        temp = zeros(nRot,nVox);
        for r = 1:nRot
            temp(r,:) = cell_rot{r}(:)';
        end
        if strcmp(pooling,'average')
            vol_filt = nanmean(temp);
        elseif strcmp(pooling,'max')
            vol_filt = nanmax(temp);
        end
        vol_filt = reshape(vol_filt',szVol(1),szVol(2),szVol(3));
    end
    
    % Getting energy image
    if energy
        szKernel = delta*2 + 1; % delta must not be NaN!
        kernel_energy = ones(szKernel,szKernel,szKernel); % Kernel used to get the energy image.
        %normFactor = imfilter(ones(size(cell_rot{r})),kernel_energy,0); % This lines produce the effective neighborhood of each voxel, similarly to NGTDM (assuming zero padding is used for the actual filter). Keep it there for reference, but will probably never be used. 
        normFactor = sum(kernel_energy(:)); % Same as Phil, and as recommended by IBSI for simplicity.
        vol_filt = imfilter(abs(vol_filt),kernel_energy,boundary) ./ normFactor; 
    end
    
end

%3. rot_invariant
function rot_invariant = set_rotInv_3D(single_rot)
    % single_rot: 1X3 cell of vectors
    
    rot_invariant = single_rot; % g_0_0
    rot_invariant = [rot_invariant;[{flip(single_rot{3})},single_rot(2),single_rot(1)]]; % g_0_pi/2_0
    rot_invariant = [rot_invariant;[{flip(single_rot{1})},single_rot(2),single_rot(3)]]; % g_0_pi_0
    rot_invariant = [rot_invariant;[single_rot(3),single_rot(2),{flip(single_rot{1})}]]; % g_0_3pi/2_0
    rot_invariant = [rot_invariant;[single_rot(2),single_rot(3),single_rot(1)]]; % g_pi/2_0_pi/2
    rot_invariant = [rot_invariant;[single_rot(2),{flip(single_rot{3})},{flip(single_rot{1})}]]; % g_pi/2_0_3pi/2
    rot_invariant = [rot_invariant;[single_rot(2),{flip(single_rot{1})},single_rot(3)]]; % g_pi/2_0_0
    rot_invariant = [rot_invariant;[{flip(single_rot{1})},{flip(single_rot{2})},single_rot(3)]]; % g_pi_0_0
    rot_invariant = [rot_invariant;[{flip(single_rot{2})},single_rot(1),single_rot(3)]]; % g_3pi/2_0_0
    rot_invariant = [rot_invariant;[{flip(single_rot{3})},{flip(single_rot{1})},single_rot(2)]]; % g_0_pi/2_3pi/2
    rot_invariant = [rot_invariant;[{flip(single_rot{3})},{flip(single_rot{2})},{flip(single_rot{1})}]]; % g_0_pi/2_pi
    rot_invariant = [rot_invariant;[{flip(single_rot{3})},single_rot(1),{flip(single_rot{2})}]]; % g_0_pi/2_pi/2
    rot_invariant = [rot_invariant;[{flip(single_rot{2})},{flip(single_rot{1})},{flip(single_rot{3})}]]; % g_pi/2_pi_0
    rot_invariant = [rot_invariant;[single_rot(1),{flip(single_rot{2})},{flip(single_rot{3})}]]; % g_pi_pi_0
    rot_invariant = [rot_invariant;[single_rot(2),single_rot(1),{flip(single_rot{3})}]]; % g_3pi/2_pi_0
    rot_invariant = [rot_invariant;[single_rot(3),{flip(single_rot{1})},{flip(single_rot{2})}]]; % g_0_3pi/2_pi/2
    rot_invariant = [rot_invariant;[single_rot(3),{flip(single_rot{2})},single_rot(1)]]; % g_0_3pi/2_pi
    rot_invariant = [rot_invariant;[single_rot(3),single_rot(1),single_rot(2)]]; % g_0_3pi/2_3pi/2
    rot_invariant = [rot_invariant;[{flip(single_rot{1})},single_rot(3),single_rot(2)]]; % g_pi_0_pi/2
    rot_invariant = [rot_invariant;[{flip(single_rot{2})},single_rot(3),{flip(single_rot{1})}]]; % g_3pi/2_0_pi/2
    rot_invariant = [rot_invariant;[single_rot(1),single_rot(3),{flip(single_rot{2})}]]; % g_0_0_pi/2
    rot_invariant = [rot_invariant;[{flip(single_rot{1})},{flip(single_rot{3})},{flip(single_rot{2})}]]; % g_pi_0_3pi/2
    rot_invariant = [rot_invariant;[{flip(single_rot{2})},{flip(single_rot{3})},single_rot(1)]]; % g_3pi/2_0_3pi/2
    rot_invariant = [rot_invariant;[single_rot(1),{flip(single_rot{3})},single_rot(2)]]; % g_0_0_3pi/2
end

%4.  getAndSaveNiiMap
function getAndSaveNiiMap(vol,kernel,boundary,info,nameSave)
volFilter = imfilter(vol,kernel,boundary,'conv');
niftiwrite(volFilter,nameSave,info,'Compressed',true);
end

%5. padVol
function vol = padVol(vol,padSize,boundary)
vol = padarray(vol,padSize,boundary);
end

%6. unPadVol
function vol = unPadVol(vol,padSize)
% Assummes option 'both' was used with padarray
padI = padSize(1);
padJ = padSize(2);
padK = padSize(3);
vol = vol((padI+1):(end-padI),(padJ+1):(end-padJ),(padK+1):(end-padK));
end


%%%%%%%%%% UTILITIES FUNCTIONS %%%%%%%%%%
function name = getPhantomName(nameWithExtension)
    name = nameWithExtension(1:end-7);
end
function info = adaptNiiHeader(info,format)
    switch format
        case 'single'
            info.raw.datatype = 16;
            info.raw.bitpix = 32;
            info.Datatype = 'single';
            info.BitsPerPixel = 32;
            info.SpaceUnits = 'Millimeter';
    end
end
function vol_out = transpose3D(vol_in)
    sizeVol = size(vol_in);
    vol_out = zeros(sizeVol(2),sizeVol(1),sizeVol(3));
    for s = 1:sizeVol(3)
        vol_out(:,:,s) = vol_in(:,:,s)';
    end
    vol_out = single(vol_out);
end

% -------------------------------------------------------------------------



% % PIECE OF CODE COMMENTED TO QUICKY UNDERSTAND HOW GABOR WORK IN MATLAB.
% % --> Copy and paste in Matlab to see what happens.
% g = gabor(4,-45,'SpatialFrequencyBandwidth',[2,1,0.5],'SpatialAspectRatio',1.5);
% figure;
% subplot(1,3,1)
% for p = 1:length(g)
% subplot(1,3,p);
% imshow(real(g(p).SpatialKernel),[]);
% lambda = g(p).Wavelength;
% sfb  = g(p).SpatialFrequencyBandwidth;
% title(sprintf('Re[h(x,y)], \\lambda = %d, \\sfb = %0.2f',lambda,sfb));
% end