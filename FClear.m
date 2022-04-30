function imNoHigh = FClear(inputImage,bEdgeMask)
% This function removes specular highlights in RGB colonoscopy images using 
% range filtering, Otsu thresholding of the parametric map and inpainting
% of the original inputImage
% -------------------------------------------------------------------------
% Input: 
% inputImage - input RGB image obtained during colonoscopy (after cropping
% in evaluation function)
%
% bEdgeMask - mask identifying the remaining part of black edges present in
% the input image (it is necessary to avoid inpainting black edges)
%
% Output:
% imNoHigh - output RGB image with removed specular highlights
% -------------------------------------------------------------------------
% Authors: Terezie Dobrovolná, Ondřej Nantl, Jan Šíma
% =========================================================================
%finding regions with very big dynamic range
pm = rangefilt(rgb2gray(inputImage),true(7)); 
% estimating threshold for extracting only specular highlights locations
T = graythresh(pm); 
% extracting only specular highlights locations
reflMask = imbinarize(imfill(pm,'holes'),T);
% inpaiting the specular highlights
imNoHigh = inpaintCoherent(inputImage,logical((~bEdgeMask).*reflMask),'SmoothingFactor',5,'Radius',5);
end