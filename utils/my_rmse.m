function output = my_rmse(ref_im, im)
    output = sqrt(sum(sum(abs(ref_im - im).^2))/numel(im));
end

