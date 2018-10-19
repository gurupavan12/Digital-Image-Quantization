clc;
clear all;
% close all;


%% Read the given Image
[I,cmap] = imread('walkbridge.tif');
Image512= I(:,:,1);  % Obtain the first layer since that contains the information in the .tif image

[rows,col]=size(Image512); % Obtain the size of the original image 

Image32 = zeros(rows/16,col/16); % Creating a zero matrix of size 32x32.

%% Downsampling the 512x512 to 32x32 image
% A 512x512 image is created with the pixels at every multiple of 16 in the rows and
% columns. The values are filled accordingly from the 32x32 image. 

i=1;j=1;
for x=1:16:rows         
    for y=1:16:col      
        Image32(i,j) = Image512(x,y);
        j=j+1;
    end
    i=i+1;
    j=1;
end

%Padding the 32x32 image with a layer of zero to ensure bilinear 
%       interpolation occurs accurately to the last set of 16pixels.
Image32(:,33)=0; 
Image32(33,:)=0;


%% Creating the 512x512 image
% A 512x512 image is created with the pixels at every multiple of 16 in the rows and
% columns. The values are filled accordingly from the 32x32 image. 

Image32= uint8((Image32));
ImageUpsampled512 =uint8(zeros(513,513));

i=1;
j=1;

[m,n]=size(Image32);

for x=1:16:513      % Jumping rows in multiples of 16
    for y=1:16:513  % Jumping coloumn in multiples of 16
        ImageUpsampled512(x,y) = Image32(i,j);
        j=j+1;
    end
    i=i+1;
    j=1;
end
% figure; imshow(ImageUpsampled512);   %Image with just the multiple of 16
%                                       pixel values amongst the rows and column.

ImageUpsampled512 = double(ImageUpsampled512);
Image32 = double(Image32);
% Interpolate the image in the x-direction
for i=1:16:513
    for j=1:1:513
        if rem(j,16) == 1
        else
            ImageUpsampled512(i,j) = (ceil(j/16)*16+1-j)/16.*Image32(ceil(i/16),ceil(j/16)) + (j-((ceil(j/16)-1)*16+1))/16*Image32(ceil(i/16),ceil(j/16)+1);
%             Implementation of bilinear interpolation in the x-direction
        end
    end
end

%  imtool((ImageUpsampled512));
a=ImageUpsampled512; %Store the x-directed interpolated image

% Interpolate the x-directed image in the y-direction.
for j=1:1:513
    for i=1:1:513
        if rem(i,16) == 1
        else
             ImageUpsampled512(i,j) = (ceil(i/16)*16+1-i)/16.*a((ceil(i/16)-1)*16+1,j) + (i-(ceil(i/16)-1)*16+1)/16*a(ceil(i/16)*16+1,j);
%             Implementation of bilinear interpolation in the y-direction

        end
    end
end

% figure; imshow(uint8(ImageUpsampled512));

%Removing the extra padded row and coloumn thereby converting the 513x513
%to a 512x512 image
for i=1:1:512
    for j=1:1:512
        if j<513
        finalImageUpsampled512(i,j) = ImageUpsampled512(i,j);
        end
    end
end

finalImageUpsampled512= uint8(finalImageUpsampled512);
figure; imshow(uint8(finalImageUpsampled512));
imwrite(finalImageUpsampled512,'BilinearInterpolated.tif');
