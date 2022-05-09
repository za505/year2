%%denoise.m
%%Author: Zarina Akbary
%%Date: 7 April 2021
%%Purpose: to denoise fluorescent images

%%User Input
basename='04062021_Exp2_colony1';
dirname=['/Users/zarina/Downloads/NYU/Year2_2021_Spring/04062021_analysis/' basename '/' basename '_AF/' basename '_aligned'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load the denoising convolution neural network
net = denoisingNetwork('DnCNN');


%change directory
cd(dirname)
directory=dir('*.tif');
T=length(directory);

% %%%This is a test
% imagename=directory(1).name;
% im=imread(imagename);
% %imshow(im)
% 
% %now denoise the image
% im2=denoiseImage(im, net);
% 
% %display them side by side
% montage({im,im2})
% title('Original Image (Left) and Denoised Image (Right)')

%now, let's try it for all the images
cd('../')

mkdir([basename '_denoised'])
cd(['./' basename '_denoised'])

    
    for t=1:T
        
        cd(dirname)
        
        imagename=directory(t).name;
        im=imread(imagename);
        b=sprintf(['%4.4d'],t);
        savename=[basename '_b' b '.tif'];
        
        im2=denoiseImage(im, net);
        
        cd('../')
        cd([basename '_denoised'])
        imwrite(im2,savename);
    end


