%% Reconstruction of In Vivo bSSFP Images Using Segregated Masks
% This is a demo that shows an example implementation of segregated sampling algorithm
% over in vivo bSSFP data with multiple phase-cycled acquisitions. Profile Encoding 
% reconstruction is applied to undersampled k-space data where the coils are combined prior to the reconstruction. 
%
% This demo is based on Senel et. al, "Statistically Segregated k-Space Sampling for Accelerating Multiple-Acquisition MRI".
% Segregated sampling is the technique that aims to minimize k-space overlap across 
% patterns for separate acquisitions while maintainnig randomness and similiar sampling
% density in individual patterns. 

% Reconstructions may take long time depending on the machine and number of
% phase cycles used ("number_of_masks") where the maximum is 8. For accs =
% 2:2:8 it may take 15-20 to complete all the reconstructions when
% pregenerated masks are used. 

clear
close all

% For better visuals
set(groot, 'defaultTextInterpreter','latex')
set(groot, 'defaultAxesTickLabelInterpreter','latex')
set(groot, 'defaultLegendInterpreter','latex')

% load data
load('data/bSSFP_in_vivo.mat')

original_image = ims;
coil_count = size(original_image,3);
im_size = [size(ims,1), size(ims,2)];

% Parameters for mask generation
lambda = 0; % segregation parameter for segregated sampling
iter = 1000; % number of trials for mask generation  
tol = im_size(1)*im_size(2)/100; % tolerance for deviation from desired sampling acceleration 
                                 % in number of sampling points (i.e. 1% deviaiton is allowed)
    
% Acceleration factors. random and segregated masks for ACC = [2,4,6,8] are
% available in the masks folder. If a different acceleration factor is used
% then parameters "power" and number_of_masks" should be adjusted.
accs = 4:2:8;

% Reconstruction Parameters 		
kSize = [11,11];  % SPIRiT kernel size
CalibSize = [96, 96]; % Calibration region size
nIterCG = 40;     % number of CG iterations for the PI part
nIterSplit = 30;  % number of splitting iterations for CS part
CalibTyk = 0.01;  % Tykhonov regularization in the calibration
ReconTyk = 1e-6;  % Tykhovon regularization in the reconstruction (SPIRiT only)
wavWeight = 1e-3; % L1-Wavelet threshold
     

% For normalization after reconstructions
upperLim = 1; % Maximum intensity values after normalization (in magnitude)
cutPer = 0.98; % Threshold percentage



% Combine coils for PE reconstruction   
delta = 1e-3;        
for phase_cycle_no = 1:size(original_image,4)
    phase_cycle_image = original_image(:,:,:,phase_cycle_no);
    for coil_no = 1:coil_count
        coil_data = fft2c(phase_cycle_image(:,:,coil_no));

        N = size(coil_data);
        phmask = zpad(hamming(16)*hamming(16)',N(1),N(2)); 
        phmask = phmask/max(phmask(:));    
        Ci(:,:,coil_no) = ifft2c(coil_data.*phmask);
    end

    imcomb = sum(phase_cycle_image.*conj(Ci)./repmat(sum(abs(Ci).^2 + delta,3), 1, 1, coil_no),3);
    combined_images(:,:,phase_cycle_no) = cut_and_normalize(imcomb,1,1);

    combined_data = fft2c(combined_images(:,:,phase_cycle_no));
    image(:,:,phase_cycle_no) = ifft2c(combined_data);
end
ims = image;

% Initialize
RMSEs = zeros(2,length(accs));
SSims = zeros(2,length(accs)); 
PSNRs = zeros(2,length(accs)); 
recon_images = zeros(im_size(1), im_size(2), 2, length(accs));

acc_no = 0;
for ACC = accs
    acc_no = acc_no + 1;
    fprintf('Reconstruction for R = %d\n',ACC)
    
    clear image
    if ACC == 2  
        for k = 1:ACC
            image(:,:,k) = cut_and_normalize(ims(:,:,4*k), 1, 1); 
            kDATA(:,:,k) = fft2c(image(:,:,k));
        end
        power = 2;
    elseif ACC == 4
        for k = 1:ACC
            image(:,:,k) = cut_and_normalize(ims(:,:,2*k), 1, 1);  
            kDATA(:,:,k) = fft2c(image(:,:,k));
        end
        power = 4;
    elseif ACC == 6 
        for k = 1:ACC
            image(:,:,k) = cut_and_normalize(ims(:,:,k), 1, 1);  
            kDATA(:,:,k) = fft2c(image(:,:,k));
        end
        power = 5;
    elseif ACC == 8
        for k = 1:ACC
            image(:,:,k) = cut_and_normalize(ims(:,:,k), 1, 1);  
            kDATA(:,:,k) = fft2c(image(:,:,k));
        end
        power = 6;
    end         
    ref_im = double(cut_and_normalize(sos(image),upperLim,cutPer));
 
    center_size = 1/(4*ACC); % Ratio of fully sampled k-space radis to maximum frequency
    no_of_masks = ACC;       % Acceleration and acquisitions are matched
    pdf = genPDF(im_size, power, 1/ACC, 2, center_size, 0);
    

    for mask_type_no = 1:2

        if mask_type_no == 1   
            load(sprintf('Masks/bSSFP_in_vivo_acc_%d_random_masks.mat',ACC))         
%             masks = genMasks(pdf, iter, tol, no_of_masks);
            fprintf('Random masks aggregate coverage: %.0f%%\n', coverage(masks)*100)
        elseif mask_type_no == 2
            load(sprintf('Masks/bSSFP_in_vivo_acc_%d_segregated_masks.mat',ACC))
%             masks = genSegregatedMasks(pdf, iter, tol, no_of_masks, lambda);
            fprintf('Segregated masks aggregate coverage: %.0f%%\n', coverage(masks)*100)
        end  
        
        pe = size(kDATA,2); fe = size(kDATA,1); coils = size(kDATA,3); % get sizes

        DATA = kDATA.*masks; % multiply with sampling matrix

        [~, dcomp] = getCalibSize(masks(:,:,1));  % get size of calibration area from mask
        DATAcomp = DATA.*repmat(dcomp,[1,1,coils]);
        scale_fctr = norm(DATAcomp(:))/sqrt(coils)/20;
        DATA = DATA/scale_fctr;


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%		 Perform Calibration                       %%%%%%%%%%
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

        im_cgspirit = ifft2c(res_cg); 
        im_cgspirit_sqr = sos(im_cgspirit);

        % Normalize the output images
        im_cgspirit_sqr = cut_and_normalize(im_cgspirit_sqr, upperLim, cutPer);

        RMSEs(mask_type_no, acc_no) = round(my_rmse(ref_im, im_cgspirit_sqr),4);
        SSims(mask_type_no, acc_no) = round(ssim(ref_im, im_cgspirit_sqr)*1e2,2);
        PSNRs(mask_type_no, acc_no) = psnr(ref_im, im_cgspirit_sqr);

        recon_images(:,:,mask_type_no,acc_no) = im_cgspirit_sqr;
    end
end

% Display the resulting reconstructions
disp_recons_acc(ref_im, recon_images, RMSEs, SSims, PSNRs, accs)