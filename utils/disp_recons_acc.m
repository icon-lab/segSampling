function  disp_recons_acc(ref_im, recon_ims, RMSEs, SSims, PSNRs, accs)
    im_size = size(ref_im);
    mask_type_count = size(recon_ims, 3);
    acc_count = length(accs);
    separation = 20;
    
    % Show reference image separately
    figure, imshow(ref_im,[]), title('Fully Sampled')
    
    image = zeros(im_size(1)*mask_type_count + separation*(mask_type_count-1), im_size(2)*acc_count + separation*(acc_count-1));

    title_str = cell(2,1);
    title_str{1} = '';
    title_str{2} = 'Rand.:  ';  
    title_str{3} = '  Seg.:  ';  
    
    
    for acc_no = 1:acc_count
        for mask_type_no = 1:mask_type_count
            image( (mask_type_no-1)*(im_size(1)+separation)+1 : (mask_type_no-1)*(im_size(1)+separation)+im_size(1), (acc_no-1)*(im_size(2)+separation)+1 : (acc_no-1)*(im_size(2)+separation)+im_size(2)) = recon_ims(:,:,mask_type_no, acc_no);
        end
        title_str{1} = [title_str{1} , sprintf('~~~~~~~~~~~~~~~~~~~~~~~~~~R = %d~~~~~~~~~~~~~~~~~~~~~~~~~~',accs(acc_no))];
        title_str{2} = [title_str{2} , sprintf('RMSE: %.4f SSim: %.2f PSNR: %.2f~~~~~~~',RMSEs(1,acc_no), SSims(1,acc_no), PSNRs(1,acc_no))];
        title_str{3} = [title_str{3} , sprintf('RMSE: %.4f SSim: %.2f PSNR: %.2f~~~~~~~',RMSEs(2,acc_no), SSims(2,acc_no), PSNRs(2,acc_no))];
    end
   
    
    figure, imshow(image,[])
    title(title_str)
    
    truesize(gcf,[700, 2500])
    ylabel('Segregated~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Random','fontweight','bold','fontsize',14),
    pos = get(gcf,'Position'); pos = [100,100,1500, 900]; set(gcf,'Position',pos);
end

