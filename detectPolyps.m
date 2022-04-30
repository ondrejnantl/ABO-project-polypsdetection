function [binaryMap] = detectPolyps(inputImage,bEdgeMask,method)
% This function detects and segments polyps in RGB colonoscopy images using 
% hysteresis thresholding and region growing technique
% 
% Image is before detection preprocessed by eliminating specular highlights
% and correction of variant lighting
% -------------------------------------------------------------------------
% Input: 
% inputImage - input RGB image obtained during colonoscopy (after cropping
% in evaluation function)
%
% bEdgeMask - binary mask identifying the remaining part of black edges 
% present in the input image (it is necessary to avoid inpainting black 
% edges) - this input is also necessary
%
% method - string representing one of our designed methods, 'HTRGRd' for
% Hough transform followed by region growing of red channel pr 'HysThRGRd'
% for hysteresis thresholding followed by region growing of red channel
% 
% Output:
% binaryMap - binary mask of segmented polyp as a matrix
% -------------------------------------------------------------------------
% Authors: Terezie Dobrovolná, Ondřej Nantl, Jan Šíma
% =========================================================================
%% Method with hysteresis thresholding and Region Growing of red channel
if strcmp(method,'HysThRGRd')
    inputImage=FClear(inputImage,bEdgeMask);
    imPrep = FLight(inputImage);
    [x,y]  = FHysThres(imPrep);
    binaryMap = FRegionGrow(imPrep,x,y);
%% Method with Hough Transform and Region Growing of red channel
elseif strcmp(method,'HTRGRd')
    inputImage=FClear(inputImage,bEdgeMask);
    imPrep = FLight(inputImage);
    [x,y]  = FHouTrans(imPrep);
    binaryMap = FRegionGrow(imPrep,x,y);
else
    error('Incorrect input for method, see help')
end
end

