%%%%%%%%%%%%%  Function nBitPlane %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:  
%      Change the gray-level quantization of the original 512x512 image by reducing the number
%       of bits per pixel from 8 to 7, 6, 5, 4, 3, 2 and 1 bits/pixel. 
%
% Input Variables:
%      I       512x512 input 2D gray-scale image to be resolved
%      x|i       x coordinate of a pixel for the input|output image
%      y|j       y coordinate of a pixel for the input|output image
%      row      Number of rows of input image
%      col      Number of coloumns of input image
%      factor   Denotes the number of bits in the plan
% 
% Returned Results:
%     ImageSpaced     new image cantaining the filtered results
%
% Processing Flow:
%      1.  Set a new image full of ZEROS
%      2.  For each valid pixel,
%             compute the mean of the 3x3 neighborhood about the
%             pixel and assign this value to the mean image
%
%  Notes:
%      This function takes an 8-bit image as input, the number of bit-plane
%       that is to be reduced to, the number of rows and columns of the input image.  
%
%  The following functions are called:
%      none
%
%       Author:      Pavan Gurudath
%       Date:        09/08/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

function ImageSliced = nBitPlane(I,factor,row,col)
    ImageSliced = uint8(zeros(512,512)); 
    i=1;j=1;
    for x=1:1:row
        for y=1:1:col
            ImageSliced(i,j) = uint8(floor(double(I(x,y))./factor)).*factor;
            j=j+1;
        end
        i=i+1;
        j=1;
    end
    ImageSliced = uint8(ImageSliced);
%     figure; imshow(ImageSliced);
end
