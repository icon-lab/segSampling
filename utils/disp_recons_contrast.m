function  disp_recons_contrast(ref_ims, recon_ims, RMSEs, SSims, PSNRs, TE_values)
    im_size = [size(ref_ims,1), size(ref_ims,2)];
    mask_type_count = size(recon_ims, 3);
    contrast_count = size(recon_ims, 4);
    separation = 20;
    
    image = zeros(im_size(1)*(mask_type_count+1) + separation*(mask_type_count), im_size(2)*contrast_count + separation*(contrast_count-1));
    
    % Show reference images on top
    for TE_no = 1:contrast_count
        image( 1:im_size(1), (TE_no-1)*(im_size(2)+separation)+1 : (TE_no-1)*(im_size(2)+separation)+im_size(2)) = ref_ims(:,:,TE_no);    
    end
    
    title_str = cell(2,1);
    title_str{1} = '';
    title_str{2} = 'Rand.:  ';  
    title_str{3} = '  Seg.:  ';  
    
    % Show reconstructed images below
    for TE_no = 1:contrast_count
        for mask_type_no = 1:mask_type_count
            image( (mask_type_no)*(im_size(1)+separation)+1 : (mask_type_no)*(im_size(1)+separation)+im_size(1), (TE_no-1)*(im_size(2)+separation)+1 : (TE_no-1)*(im_size(2)+separation)+im_size(2)) = recon_ims(:,:,mask_type_no, TE_no);
        end
        title_str{1} = [title_str{1} , sprintf('~~~~~~~~~~~~~~~~~~~~~~~~~~TE = %d ms~~~~~~~~~~~~~~~~~~~~~~~~~~',TE_values(TE_no))];
        title_str{2} = [title_str{2} , sprintf('RMSE: %.4f SSim: %.2f PSNR: %.2f~~~~~~~',RMSEs(1,TE_no), SSims(1,TE_no), PSNRs(1,TE_no))];
        title_str{3} = [title_str{3} , sprintf('RMSE: %.4f SSim: %.2f PSNR: %.2f~~~~~~~',RMSEs(2,TE_no), SSims(2,TE_no), PSNRs(2,TE_no))];
    end
   
    
    figure, imshow(image,[])
    title(title_str)
    
    truesize(gcf,[800,800])
    ylabel('Segregated~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Random~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Fully Sampled','fontweight','bold','fontsize',14),
    pos = get(gcf,'Position'); pos = [200,0,1400, 1000]; set(gcf,'Position',pos);
end

