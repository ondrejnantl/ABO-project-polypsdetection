function [binaryMap] = detectPolyps(inputImage,bEdgeMask)
% UNTITLED Summary of this function goes here
% Detailed explanation goes here
% 
% Authors: Ondřej Nantl, Terezie Dobrovolná, Jan Šíma
% =========================================================================
%% elimination of specular highlights and correction of variant lighting
% elimination of specular highlights
pm = rangefilt(rgb2gray(inputImage),true(7));
T = graythresh(pm);
reflMask = imbinarize(imfill(pm,'holes'),T);
imCropped = inpaintCoherent(inputImage,logical((~bEdgeMask).*reflMask),'SmoothingFactor',5,'Radius',5);

% correction of variant lighting
[m,n,o] = size(imCropped);
mm = zeros(m,n,o);
N = 19;
meanMask = 1/(N^2).*ones(N,N);
for j = 1:o
    mm(:,:,j) = 0.15.*conv2(imCropped(:,:,j),meanMask,'same'); % slight change in constant compared to Sanchez2018
end
imPrep = imCropped - mm;

%% Hough transform for circles
imEdge = edge(rgb2gray(imPrep),'canny',[.03 .1],sqrt(2)); % constants set according to Sanchez2018
rs = 5:40; % range of diameters
HS = zeros(size(imPrep,1),size(imPrep,2),length(rs));
r_ind = 1;
[X,Y] = find(imEdge == 1);
for r = rs
    tmp_c = gen_circle(r);
    for i = 1:length(X)
        c1 = X(i);
        c2 = Y(i);
        if c1 > r && c1< (size(imPrep,1) - r)
            if c2 > r && c2< (size(imPrep,2) - r)
                HS((c1-r):(c1+r),(c2-r):(c2+r),r_ind) = HS((c1-r):(c1+r),(c2-r):(c2+r),r_ind)+tmp_c;
            end
        end
    end
    r_ind = r_ind + 1;
end
% finding the center of the most probable circle in edge representation
[linInd] = find(HS == max(HS,[],'all'));
[y,x,r] = ind2sub(size(HS),linInd); 

if length(x)>1 || length(y)>1
    x = floor(mean(x));
    y = floor(mean(y));
end
% 
%% region growing
for i = 1:o
segIm(:,:,i) = grayconnected(imPrep(:,:,i),y,x,0.02); % position is defined by Hough t.
end
%% geometric contours
% finding the smallest object in 3 results of region growing
sumRegion = reshape(sum(sum(segIm)),[3 1 1]);
[~,smallObjChannel] = min(sumRegion);
% level sets
binaryMap = activecontour(rgb2gray(imCropped),imdilate(segIm(:,:,smallObjChannel),[1 1 1; 1 1 1; 1 1 1]));
end

