%%%%%%%%%%%%% main.m file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:  
%      Familiarize concepts of digital image processing 
%
%
% Code flow: 
%      1.  Load input image of size 512x512.
%      2.  Downsample the spatial resolution of the input image to 256x256, 128x128, and 32x32
%          pixels. 
%      3.  Interpolate the obtained images in step 2 back to 512x512. 
%      3.  Create an interpolated 512x512 image from the downsampled 32x32 image by using
%          bilinear interpolation.
%      4.  Change the gray-level quantization of the original 512x512 image by reducing the number
%          of bits per pixel from 8 to 7, 6, 5, 4, 3, 2 and 1 bits/pixel.
%      5.  Create a 512x512 image that changes the spatial resolution to 256x256 pixels and gray-scale
%          resolution to 6 bits/pixel.
% 
%  The following functions are called:
%      BilinearInterpolation    To perform bilinear interpolation
%      nBitPlane                To acheive gray level quantization 
%          
%
%  Author:      Pavan Gurudath
%  Date:        09/08/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;      % Clear out all memory 

clc;
clear all;
close all;


%% Read the given Image
[I,cmap] = imread('walkbridge.tif');
Image512= I(:,:,1);  % Obtain the first layer since it contains the necessary information in the .tif image

I=Image512;
Itemp=I;            %Storing the input image such that it can be used later
[rows,col]=size(Image512); % Obtain the size of the original image 

%% Objective 1: To downsample the spatial resolution of a 512x512 image to different sizes.

[rows,col]=size(Image512); % Obtain the size of the original image 


%% Create a new image that is reduced from 512x512 to 256x256

Image256 = zeros(rows/2,col/2); % Create a 256x256 image with zeroes
% i and j are coordinates of Image256 
% x and y are coordinates of Image512
i=1;j=1;              
for x=1:2:rows         
    for y=1:2:col
        Image256(i,j) = Image512(x,y);
        j=j+1;
    end
    i=i+1;
    j=1;
end

Image256 = uint8(Image256);   %Converting the image into an unsigned 8bit integer format.

figure; imshow(Image256);
imwrite(Image256,'Image256.tif'); %Writing the image 

%% Create a new image that is reduced from 512x512 to 128x128
Image128 = zeros(rows/4,col/4); % Create a 128x128 image with zeroes
% i and j are coordinates of Image128 
% x and y are coordinates of Image512
i=1;j=1;
for x=1:4:rows
    for y=1:4:col
        Image128(i,j) = Image512(x,y);
        j=j+1;
    end
    i=i+1;
    j=1;
end

Image128 = uint8(Image128);   %Converting the image into an unsigned 8bit integer format.

figure; imshow(Image128);
imwrite(Image128,'Image128.tif'); %Writing the image 


%% Create a new image that is reduced from 512x512 to 32x32
Image32 = zeros(rows/16,col/16); % Create a 32x32 image with zeroes
% i and j are coordinates of Image32 
% x and y are coordinates of Image512
i=1;j=1;
for x=1:16:rows
    for y=1:16:col
        Image32(i,j) = Image512(x,y);
        j=j+1;
    end
    i=i+1;
    j=1;
end

Image32 = uint8(Image32);   %Converting the image into an unsigned 8bit integer format.

figure; imshow(Image32);
imwrite(Image32,'Image32.tif'); %Writing the image 






%% Objective 1.2: To upsample the downsampled image using nearest-neighbor interpolation i.e. from a 256x256 to 512x512.

Image256resample = uint8(Image256); %Making sure the image is in unsigned 8bit integer format.

[M, N] = size(Image256resample);    %Number of rows and columns of the image from which it is being upsampled.

ImageUpsampled256 = zeros(512,512); %Create a 512x512 zero image
i=1;
j=1;
for x =  1 : M
    for y =  1 : N
        for i1 = 0 : 1
            for j1 = 0 : 1
                ImageUpsampled256(i+i1,j+j1) = Image256resample(x,y);
            end
        end
        j=j+2;
    end
    i=i+2;
    j=1;
end
ImageUpsampled256 = uint8(ImageUpsampled256);   %Converting the image into an unsigned 8bit integer format.

figure, imshow(ImageUpsampled256);
imwrite(ImageUpsampled256,'ImageUpsampled256.tif');




%% Upsample the 128x128 image to 512x512

Image128resample = uint8(Image128);%Make sure the image is in unsigned 8bit integer format.

[M, N] = size(Image128resample);%Number of rows and columns of the image from which it is being upsampled.

ImageUpsampled128 = zeros(512,512); %Create a zero 512x512 image

i=1;
j=1;
for x =  1 : M
    for y =  1 : N
        for i1 = 0 : 3
            for j1 = 0 : 3
                ImageUpsampled128(i+i1,j+j1) = Image128resample(x,y);
            end
        end
        j=j+4;
    end
    i=i+4;
    j=1;
end

ImageUpsampled128 = uint8(ImageUpsampled128);   %Converting the image into an unsigned 8bit integer format.

figure, imshow(ImageUpsampled128);
imwrite(ImageUpsampled128,'ImageUpsampled128.tif');

%% Upsample the 32x32 image to 512x512
Image32resample = uint8(Image32);%Make sure the image is in unsigned 8bit integer format.

[M, N] = size(Image32resample);%Number of rows and columns of the image from which it is being upsampled.

ImageUpsampled32 = zeros(512,512); %Create a zero 512x512 image 
i=1;
j=1;
for x =  1 : M
    for y =  1 : N
        for i1 = 0 : 15
            for j1 = 0 : 15
                ImageUpsampled32(i+i1,j+j1) = Image32resample(x,y);
            end
        end
        j=j+16;
    end
    i=i+16;
    j=1;
end
ImageUpsampled32 = uint8(ImageUpsampled32);%Converting the image into an unsigned 8bit integer format.

figure(), imshow(uint8(ImageUpsampled32));
imwrite(ImageUpsampled32,'ImageUpsampled32.tif');

%% Objective 2: To upsample the downsampled 32x32 image using bilinear interpolation to a 512x512 image

BilinearInterpolation;
%% Objective 3: To reduce the gray-level quantization of the image by reducing the number of bits per pixel.

I=imread('walkbridge.tif');
I=I(:,:,1);
[row col] = size(I);
%% 7bit plane
Inew7=nBitPlane(I,2,row,col);   %Calling nBitPlane function to do the necessary processing
figure; imshow(Inew7);
imwrite(Inew7, '7BitPlane.tif');

%% 6bit plane
Inew6=nBitPlane(I,4,row,col);   %Calling nBitPlane function to do the necessary processing
figure; imshow(Inew6);
imwrite(Inew6, '6BitPlane.tif');
 
%% 5bit plane
Inew5=nBitPlane(I,8,row,col);   %Calling nBitPlane function to do the necessary processing
figure; imshow(Inew5);
imwrite(Inew5, '5BitPlane.tif');
 %% 4bit plane
Inew4=nBitPlane(I,16,row,col);   %Calling nBitPlane function to do the necessary processing
figure; imshow(Inew4);
imwrite(Inew4, '4BitPlane.tif');
%% 3bit plane
Inew3=nBitPlane(I,32,row,col);   %Calling nBitPlane function to do the necessary processing
figure; imshow(Inew3);
imwrite(Inew3, '3BitPlane.tif');
 
%% 2bit plane
Inew2=nBitPlane(I,64,row,col);   %Calling nBitPlane function to do the necessary processing
figure; imshow(Inew2);
imwrite(Inew2, '2BitPlane.tif');

%% 1bit plane
Inew1=nBitPlane(I,128,row,col);   %Calling nBitPlane function to do the necessary processing
figure; imshow(Inew1);
imwrite(Inew1, '1BitPlane.tif');


%% Objective 4: To infer the changes compared with the original image when both spatial and gray-level resolution are changed.

ImageSpatial256 = zeros(rows/2,col/2); % Create a 256x256 image with zeroes
% i and j are coordinates of ImageSpatial256 
% x and y are coordinates of Image512
i=1;j=1;              
for x=1:2:rows         
    for y=1:2:col
        ImageSpatial256(i,j) = Image512(x,y);
        j=j+1;
    end
    i=i+1;
    j=1;
end

ImageSpatial256 = uint8(ImageSpatial256);   %Converting the image into an unsigned 8bit integer format.
imwrite(ImageSpatial256, 'SpatialResolution only.tif');

[rowgray, colgray] = size(ImageSpatial256);
igray=1;
jgray=1;
Inewgray = uint8(zeros(256,256));
for xgray=1:1:rowgray
    for ygray=1:1:colgray
        Inewgray(igray,jgray) = uint8(floor(double(ImageSpatial256(xgray,ygray))./2).*2);
        jgray=jgray+1;
    end
    igray=igray+1;
    jgray=1;
end
figure;imshow(Inewgray);
imwrite(Inewgray,'Spatial and gray resolution.tif');

