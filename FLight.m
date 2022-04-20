function imPrep = FLight(inputImage)
% This function adjust lighting in RGB colonoscopy images using 
% subtracting weighted image of local mean from 20x20 px neighborhood
% -------------------------------------------------------------------------
% Input: 
% inputImage - input RGB image obtained during colonoscopy (after cropping
% in evaluation function and removal of specular highlights)
%
% Output:
% imPrep - output RGB image with adjusted lighting
% -------------------------------------------------------------------------
% Authors: Terezie Dobrovolná, Ondřej Nantl, Jan Šíma
% =========================================================================
[m,n,o] = size(inputImage);
mm = zeros(m,n,o);
N = 20; % size of the local neighborhood
meanMask = 1/(N^2).*ones(N,N);
% calculating local mean in 20x20 window
for j = 1:o
    mm(:,:,j) = 0.3.*conv2(inputImage(:,:,j),meanMask,'same'); % slight change in constant compared to Sanchez2018
end
% subtracting local mean image
imPrep = inputImage - mm;
% transforming into different color systems
% imPrepLab = rgb2lab(imPrep);
% imPrepGray = rgb2gray(imPrep);
% imPrepHSV = rgb2hsv(imPrep);
end