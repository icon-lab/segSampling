function [ out_im ] = cut_and_normalize( image, upperLim, cutPer)

% This function takes the absolute value of the image, cuts the pixels
% above cutPer percentage and normalizes to the interval [0, upperLim]
    complex = 1;
    
    if complex
        abs_image = abs(image);   
        out_im = zeros(size(image));
        x_size = size(image, 1);
        y_size = size(image, 2);

        len = x_size*y_size;

        if size(image,3) == 1

            x = reshape(abs_image,len, 1);
            x = sort(x);

            cutValue = x(round(len*cutPer));

            image(abs_image>cutValue) = image(abs_image>cutValue)./abs(image(abs_image>cutValue)).*cutValue; %maximum absolute value becomes cut value

            image = image/cutValue*upperLim; % To make maximum value 1

            out_im = image;

        else
            for k = 1:size(image,3)
                abs_im = abs_image(:,:,k);

                x = reshape(abs_im,len, 1);
                x = sort(x);
                cutValue = x(round(len*cutPer));
                temp_im = image(:,:,k);
                temp_im(abs_im>cutValue) = temp_im(abs_im>cutValue)./abs(temp_im(abs_im>cutValue)).*cutValue;
                temp_im = temp_im/cutValue*upperLim;
                out_im(:,:,k) = temp_im;       
            end
        end  
        
        
        
    else
        image = abs(image);   
        out_im = zeros(size(image));
        x_size = size(image, 1);
        y_size = size(image, 2);

        len = x_size*y_size;

        if size(image,3) == 1

            x = reshape(image,len, 1);
            x = sort(x);

            cutValue = x(round(len*cutPer));

            image(image>cutValue) = image(image>cutValue)./abs(image(image>cutValue)).*cutValue;

            image = image/cutValue*upperLim;

            out_im = image;

        else
            for k = 1:size(image,3)
                im = image(:,:,k);

                x = reshape(im,len, 1);
                x = sort(x);
                cutValue = x(round(len*cutPer));
                temp_im = image(:,:,k);
                temp_im(im>cutValue) = temp_im(im>cutValue)./abs(temp_im(im>cutValue)).*cutValue;
                temp_im = temp_im/cutValue*upperLim;
                out_im(:,:,k) = temp_im;       
            end
        end  
        
        
    end
end

