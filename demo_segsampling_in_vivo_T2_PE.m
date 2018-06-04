%% Reconstruction of In Vivo T2 weighted data using segregated masks
% This is a demo that shows an example implementation of segregated sampling algorithm
% over in vivo T2 weighted data with three differnt TE values (three contrast). 
% Profile Encoding reconstruction is applied to undersampled k-space data where the coils are combined prior to the reconstruction. 
%
% This demo is based on Senel et. al, "Statistically Segregated k-Space Sampling for Accelerating Multiple-Acquisition MRI".
% Segregated sampling is the technique that aims to minimize k-space overlap across 
% patterns for separate acquisitions while maintainnig randomness and similiar sampling
% density in individual patterns.

clear
close all

% For better visuals
set(groot, 'defaultTextInterpreter','latex')
set(groot, 'defaultAxesTickLabelInterpreter','latex')
set(groot, 'defaultLegendInterpreter','latex')

% load data
load('T2_weighted_in_vivo.mat')

im_size = [size(ims,1), size(ims,2)];
TEs = [145,257,320];

% Parameters for PDF and mask generation
lambda = 0;                     % segregation parameter for segregated sampling
iter = 1000;                    % number of trials for mask generation  
ACC = 4;                        % Undersampling factor
no_of_masks = 3;                % 3 masks for the 3 contrast data
center_size = 1/(4*ACC);        % Ratio of fully sampled k-space radis to maximum frequency
tol = im_size(1)*im_size(2)/100;% tolerance for deviation from desired sampling acceleration
 
% random and segregated masks for ACC = [3,4,6] are
% available in the masks folder. If a different acceleration factor is used
% then parameters "power"  should be adjusted.

% power of the polynomial used to generate PDF
switch ACC
    case 3
        power = 3;
    case 4 
        power = 4;
    case 6 
        power = 5;
end

% Generate initial PDF for sampling masks
pdf = genPDF(im_size, power, 1/ACC, 2, center_size, 0);

% Combine coils for PE reconstruction
delta = 1e-11;
for TE = 1:size(ims,4)
    for coil_no = 1:size(ims,3)
        coil_data = fft2c(ims(:,:,coil_no,TE));

        N = size(coil_data);
        phmask = zpad(hamming(16)*hamming(16)',N(1),N(2)); 
        phmask = phmask/max(phmask(:));    
        Ci(:,:,coil_no) = ifft2c(coil_data.*phmask);
    end
    imcomb = sum(ims(:,:,:,TE).*conj(Ci)./repmat(sum(abs(Ci).^2 + delta,3), 1, 1, coil_no),3);
    combined_images(:,:,TE) = cut_and_normalize(imcomb,1,1);
end

clear ims
ims = combined_images;
kDATA = fft2c(combined_images);
 
% Generate initial PDF for sampling masks
pdf = genPDF(im_size, power, 1/ACC, 2, center_size, 0);
 

% Reconstruction Parameters 		
kSize = [11,11];  % SPIRiT kernel size
CalibSize = [96, 96]; % Calibration region size
nIterCG = 15;     % number of CG iterations for the PI part
nIterSplit = 30;  % number of splitting iterations for CS part
CalibTyk = 0.01;  % Tykhonov regularization in the calibration
ReconTyk = 1e-6;  % Tykhovon regularization in the reconstruction (SPIRiT only)
wavWeight = 5e-3; % L1-Wavelet threshold


% For normalization after reconstructions
upperLim = 1;
cutPer = 0.98;

% Initialize 
RMSEs = zeros(2,3);
SSims = zeros(2,3); 
PSNRs = zeros(2,3); 
recon_images = zeros(im_size(1), im_size(2), 2, 3);

ref_im = cut_and_normalize(abs(ims), upperLim, cutPer);
for mask_type_no = 1:2
    if mask_type_no == 1  
        load(sprintf('T2_weighted_in_vivo_acc_%d_random_masks.mat',ACC))
%         masks = genMasks(pdf, iter, tol, no_of_masks);
        fprintf('Random masks aggregate coverage: %.0f%%\n', coverage(masks)*100)
    elseif mask_type_no == 2
        load(sprintf('T2_weighted_in_vivo_acc_%d_segregated_masks.mat',ACC))
%         masks = genSegregatedMasks(pdf, iter, tol, no_of_masks, lambda);
        fprintf('Segregated masks aggregate coverage: %.0f%%\n', coverage(masks)*100)
    end 

    pe = size(kDATA,2); fe = size(kDATA,1); coils = size(kDATA,3); % get sizes

    DATA = kDATA.*masks; % multiply with sampling matrix

    [~, dcomp] = getCalibSize(masks(:,:,1));  % get size of calibration area from mask
    DATAcomp = DATA.*repmat(dcomp,[1,1,coils]);
    scale_fctr = norm(DATAcomp(:))/sqrt(coils)/20;
    DATA = DATA/scale_fctr;

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%        Perform Calibration                       %%%%%%%%%%
    disp('performing calibration for SPIRiT')
    kCalib = crop(DATA,[CalibSize,coils]);

    kernel = calibSPIRiT(kCalib, kSize, coils, CalibTyk);
    GOP = SPIRiT(kernel, 'fft',[fe,pe]);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%                  Reconstruction                        %%%%%%%%%

    disp('performing CG reconstruction')
    tic;
    [res_cg, RESVEC] = cgL1SPIRiT(DATA,GOP,nIterCG,ReconTyk, wavWeight, nIterSplit);
    toc 

    im_cgspirit = abs(ifft2c(res_cg));

    % Normalize the output images
    im_cgspirit = cut_and_normalize(im_cgspirit, upperLim, cutPer);

    for T2_no = 1:3
        RMSEs(mask_type_no, T2_no) = round(my_rmse(ref_im(:,:,T2_no), im_cgspirit(:,:,T2_no)),4);
        SSims(mask_type_no, T2_no) = round(ssim(ref_im(:,:,T2_no), im_cgspirit(:,:,T2_no))*1e2,2);
        PSNRs(mask_type_no, T2_no) = psnr(ref_im(:,:,T2_no), im_cgspirit(:,:,T2_no));

        recon_images(:,:,mask_type_no,T2_no) = im_cgspirit(:,:,T2_no);
    end
end

% Display the resulting reconstructions
disp_recons_contrast(ref_im, recon_images, RMSEs, SSims, PSNRs, TEs)