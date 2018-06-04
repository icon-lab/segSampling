function [ masks, PDFs, Intrs ] = genSegregatedMasks(PDF, iter, tol, numOfMasks, lambda, iter2, generatedMasks, generatedPDFs, generatedIntrs)
    % Using given PDF generates segregated sampling patterns in a recursive manner.
    %
    % INPUTS  
    % PDF, iter, tol - passed to mySamplingPattern fucntion.
    % numOfMasks - number of segregated masks to be generated
    % lambda - segregation parameter. 0 is maximum segregation, 1 is no
    % segregation
    % other inputs are used for recursion.
    %
    % OUTPUTS
    % masks - generated sampling masks
    % PDFs - contains generated PDFs where the first PDF is the input PDF
    % Intrs - measured interference for the generated mask obtained from 
    % segSamplingPattern function
    %
    % Note: This might not be the optimal implementation of this algorithm.

    if nargin < 6
        maskInit = false(size(PDF,1), size(PDF,2), numOfMasks);
        PDFInit = NaN(size(PDF,1), size(PDF,2), numOfMasks);
        IntrInit = NaN(1,numOfMasks);
        [masks, PDFs, Intrs] = genSegregatedMasks(PDF, iter, tol, numOfMasks, lambda, 1, maskInit, PDFInit, IntrInit);
        
    elseif iter2 <= numOfMasks
        masks = generatedMasks;
        PDFs = generatedPDFs;
        Intrs = generatedIntrs;
        
        %Create mask with given PDF and save the PDF
        PDFs(:,:,iter2) = PDF;
        [masks(:,:,iter2), stat, N] = segSamplingPattern(PDF, iter, tol, PDFs(:,:,1));
        Intrs(iter2) = stat(N);
        
        sumMask  = zeros(size(masks, 1), size(masks, 2));
        for i = 1:iter2
            sumMask = sumMask + masks(:,:,i);
        end
        sumMask = logical(sumMask);
               
        newPDF = zeros(size(PDF,1), size(PDF,2));

        % Calculate new PDF using generated patterns
        for i= 1:size(newPDF,1) % i is rows
            for j = 1:size(newPDF,2) % j is columns
                prob = PDFs(i, j, 1); %Probability of the selected point in the initial PDF
                
                if prob < 1
                    num = prob;
                    den = 1;
                    check = 0;
                    sumNum = num;
                    
                    for k = 1:iter2+1 %To check if all points are sampled or all points will be sampled in this mask or there are  a lot more points to be sampled                   
                        lim = num/den;                  
                        if lim > 1 || lim < 0
                            check = check + 1;
                        end                       
                        if check == 1
                            den1 = den;
                        end                       
                        den = den-num;
                        num = prob-prob*lambda*(sumNum);
                        sumNum = sumNum + num;
                    end
                    
                    if (check < 2) %If not all points with this probability are sampled                    
                        if (lim <= 1) && (sumMask(i, j) == 1)%If point is sampled and there are many unsampled points
                            prob  = prob*lambda;
                        elseif (lim <= 1) && (sumMask(i, j) == 0)%If point is not sampled and there are many unsampled points
                            prob = lim;
                        elseif (lim > 1) && (sumMask(i, j) == 1)%If point is sampled and there are only few unsampled points
                            prob = (prob-den1)/(1-den1);
                        elseif (lim > 1) && (sumMask(i, j) == 0)%If point is not sampled and there are only few unsampled points
                            prob = 1;
                        else
                            disp('There is something wrong');
                        end
                    end
                end
          
                newPDF(i, j) = prob;
            end
        end
               
        %Recall fuction with new PDF and all masks created
        [masks,PDFs, Intrs] = genSegregatedMasks(newPDF, iter, tol, numOfMasks, lambda, iter2+1, masks, PDFs, Intrs);
    else
        masks = generatedMasks;
        PDFs = generatedPDFs;
        Intrs = generatedIntrs;
    end
end
