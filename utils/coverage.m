function [coverage] = coverage( masks )
    total = size(masks,1)*size(masks,2);
    sum = zeros(size(masks,1), size(masks,2));
    
    for i = 1:size(masks, 3)
        sum = sum + masks(:,:,i);
    end
    coverage = nnz(sum)/total;
end

