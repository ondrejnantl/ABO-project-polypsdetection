function [final,areaMean,areaSize] = FRegionGrow(imPrep,x,y)
% This function performs segmentation of polyp in preprocessed RGB 
% colonoscopy image using region growing technique with seed coordinates
% obtained using function FHysThres
% -------------------------------------------------------------------------
% Input: 
% imPrep - input preprocessed RGB image obtained during colonoscopy (after 
% cropping in evaluation function, removal of specular highlights and 
% adjustment of lighting)
%
% x - x coordinate of seed for region growing
% 
% y - y coordinate of seed for region growing
% 
% colorID - index of color component for performing region growing, 1 - red,
% 2 - green,3 - blue
% 
% Output:
% final - final mask of segmented polyp after region growing and final 
% adjustments
% 
% areaMean - mean value of area which is assigned as potential polyp region
% 
% areaSize - size of area which is assigned as potential polyp region
% -------------------------------------------------------------------------
% Authors: Terezie Dobrovolná, Ondřej Nantl, Jan Šíma
% =========================================================================
[m,n,~] = size(imPrep);
% segmenting using red component of an input image
segIm = zeros(size(imPrep(:,:,1)));
% threshold estimated from standard deviation
Trg = 0.8*std(imPrep(:,:,1),[],'all');
% change threshold and segment again until the object is at least of 0.00005
% of size on the input image
while sum(segIm == 1)< 0.00005*m*n
segIm = grayconnected(imPrep(:,:,1),y,x,Trg);
Trg = 1.25*Trg;
end

% final adjustments
final = imfill(segIm,'holes');

%% features for data analysis
imPrepCol = imPrep(:,:,1);
areaPixels = imPrepCol(final);
%pixel size of segmented area
areaSize = sum(final(:));
% mean of brightness of segmented area
areaMean = mean(areaPixels(:));
end
