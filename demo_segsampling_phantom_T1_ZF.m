%% Reconstruction of T1 weighted numerical phantom data using segregated masks
% This is a demo that shows an example implementation of segregated sampling algorithm
% over T1 weighted numerical phantom data for a single TE value. 
% Zero filled (ZF) Fourier reconstruction is applied to single channel undersampled k-space
% where data is compansated for variable sampling density. 
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
load('data/T1_Weighted_simulation.mat')
image = double(cut_and_normalize(ims,1,1)); % Normalize complex image to 0-1 range in magnitude
kspace = fft2c(ims);                        % Obtain k-space

im_size = size(image);

% Parameters for mask generation
lambda = 0; % segregation parameter for segregated sampling
iter = 1000; % number of trials for mask generation  
tol = im_size(1)*im_size(2)/100; % tolerance for deviation from desired sampling acceleration 
                                 % in number of sampling points (i.e. 1% deviaiton is allowed)
                                 
% Acceleration factors. random and segregated masks for ACC = [2,4,6,8] are
% available in the masks folder. If a different acceleration factor is used
% then parameters "power" and number_of_masks" should be adjusted.                     
accs = 2:2:8; % If different acceleration factors are used variables "power" and "no_of_masks" should be adjusted

% Initialize
RMSEs = zeros(2,length(accs));
SSims = zeros(2,length(accs)); 
PSNRs = zeros(2,length(accs)); 
recon_images = zeros(im_size(1), im_size(2), 2, length(accs));

acc_no = 0;
for ACC = accs
    fprintf('Reconstruction for R = %d\n',ACC)
    acc_no = acc_no + 1;
    switch ACC
        case 2
            power = 2;
        case 4
            power = 4;
        case 6 
            power = 5;
        case 8
            power = 6;
    end
    
    no_of_masks = ACC;                  % same data is acquired multiple (ACC) times with differnt masks
    center_size = 1/(4*ACC);            % Ratio of full sampled center radius to maximum spatial frequency 

    pdf = genPDF(im_size, power, 1/ACC, 2, center_size, 0);
    
    % Load the pre-generaetd masks
    load(sprintf('T1_weighted_acc_%d_random_masks.mat',ACC))
    load(sprintf('T1_weighted_acc_%d_segregated_masks.mat',ACC))
    
%     % Generate sampling masks
%     rand_masks = genMasks(pdf, iter, tol, no_of_masks);
%     seg_masks = genSegregatedMasks(pdf, iter, tol, no_of_masks, 1);
      
    fprintf('Random masks aggregate coverage: %.0f%%\n', coverage(rand_masks)*100)
    fprintf('Segregated masks aggregate coverage: %.0f%%\n', coverage(seg_masks)*100)
    
    rand_Sum_Recon = zeros(im_size);
    seg_Sum_Recon = zeros(im_size);

    rand_ZF_Recons = zeros(im_size(1),im_size(2),no_of_masks);
    seg_ZF_Recons = zeros(im_size(1),im_size(2),no_of_masks);
    
    % Obtain ZF reconstruction with density compensation and combine them by taking their average
    for k = 1:no_of_masks
        rand_ZF_Recons(:,:,k) = ifft2c(kspace.*rand_masks(:,:,k)./pdf);
        rand_Sum_Recon = rand_Sum_Recon + rand_ZF_Recons(:,:,k);
        
        seg_ZF_Recons(:,:,k) = ifft2c(kspace.*seg_masks(:,:,k)./pdf);
        seg_Sum_Recon = seg_Sum_Recon + seg_ZF_Recons(:,:,k);
    end

    rand_Sum_Recon = abs(rand_Sum_Recon)./no_of_masks;   
    rand_Sum_Recon = cut_and_normalize(rand_Sum_Recon, 1, 1); % Adjust the scale
        
    seg_Sum_Recon = abs(seg_Sum_Recon)./no_of_masks;   
    seg_Sum_Recon = cut_and_normalize(seg_Sum_Recon, 1, 1);
    
    recon_images(:,:,1,acc_no) = abs(rand_Sum_Recon);
    recon_images(:,:,2,acc_no) = abs(seg_Sum_Recon);
    
    RMSEs(1,acc_no) = round(sqrt(my_rmse(abs(image), abs(rand_Sum_Recon))),4);
    SSims(1,acc_no) = round(ssim(abs(image), abs(rand_Sum_Recon))*1e2,2);
    PSNRs(1,acc_no) = psnr(abs(image), abs(rand_Sum_Recon));
    
    RMSEs(2,acc_no) = round(sqrt(my_rmse(abs(image), abs(seg_Sum_Recon))),4);
    SSims(2,acc_no) = round(ssim(abs(image), abs(seg_Sum_Recon))*1e2,2);
    PSNRs(2,acc_no) = psnr(abs(image), abs(seg_Sum_Recon));
end

% Display the resulting reconstructions
disp_recons_acc(abs(image), recon_images, RMSEs, SSims, PSNRs, accs)