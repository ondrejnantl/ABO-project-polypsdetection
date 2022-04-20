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
% clear bEdgeMask bEdgeMask2 bEdgeMask3 imCropped imCroppedRow GTCropped GTCroppedRow
pm = rangefilt(rgb2gray(inputImage),true(7)); %finding regions with very big dynamic range
T = graythresh(pm); % estimating threshold for extracting only specular highlights locations
reflMask = imbinarize(imfill(pm,'holes'),T); % extracting only specular highlights locations
imNoHigh = inpaintCoherent(inputImage,logical((~bEdgeMask).*reflMask),'SmoothingFactor',5,'Radius',5); % inpaiting the specular highlights
end