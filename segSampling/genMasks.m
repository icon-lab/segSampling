function [ masks, Intrs ] = genMasks( PDF, iter, tol, numOfMasks )
    masks = false(size(PDF,1), size(PDF,2), numOfMasks); 
    Intrs = zeros(1,numOfMasks); 
    for i = 1:numOfMasks
        [masks(:,:,i), stat, N] = samplingPattern(PDF, iter, tol);
        Intrs(i) = stat(N);
    end
end

