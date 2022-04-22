function [binaryMap] = detectPolyps(inputImage,bEdgeMask)
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
% bEdgeMask - mask identifying the remaining part of black edges present in
% the input image (it is necessary to avoid inpainting black edges)
%
% Output:
% binaryMap - binary mask of segmented polyp as a matrix
% -------------------------------------------------------------------------
% Authors: Terezie Dobrovolná, Ondřej Nantl, Jan Šíma
% =========================================================================
%% elimination of specular highlights and correction of variant lighting

% % elimination of specular highlights
% pm = rangefilt(rgb2gray(inputImage),true(7)); %finding regions with very big dynamic range
% T = graythresh(pm); % estimating threshold for extracting only specular highlights locations
% reflMask = imbinarize(imfill(pm,'holes'),T); % extracting only specular highlights locations
% imCropped = inpaintCoherent(inputImage,logical((~bEdgeMask).*reflMask),'SmoothingFactor',5,'Radius',5); % inpaiting the specular highlights
% 
% % correction of variant lighting
% [m,n,o] = size(imCropped);
% mm = zeros(m,n,o);
% N = 20;
% meanMask = 1/(N^2).*ones(N,N);
% % calculating local mean in 20x20 window
% for j = 1:o
%     mm(:,:,j) = 0.3.*conv2(imCropped(:,:,j),meanMask,'same'); % slight change in constant compared to Sanchez2018
% end
% % subtracting mean image
% imPrep = imCropped - mm;
% % transforming into different color systems
% imPrepLab = rgb2lab(imPrep);
% imPrepGray = rgb2gray(imPrep);
% imPrepHSV = rgb2hsv(imPrep);
%% hysteresis thresholding
% % using red component of an image
% T = multithresh(imPrep(:,:,1),2);
% BW = hysthresh(imPrep(:,:,1),T(2),T(1));
% % using grayscale image
% T = multithresh(imPrepGray,2);
% BW = hysthresh(imPrepGray,T(2),T(1));
% % using lab space
% % T = multithresh(imPrepLab(:,:,2),2);
% % BW = hysthresh(imPrepLab(:,:,2),T(2),T(1));
% % BW = hysthresh(imPrepLab(:,:,2),0.75*max(imPrepLab(:,:,2),[],'all'),0.65*max(imPrepLab(:,:,2),[],'all'));
% 
% % BW = imerode(BW,strel('disk',2));
% props = regionprops(BW,'Area','Centroid');%,'Circularity','ConvexHull','ConvexImage','FilledImage','MajorAxisLength','MinorAxisLength'
% [~,idx] = sort([props.Area],'descend');
% biggest = idx(1);
% seedRow = round(props(biggest).Centroid(2));
% seedCol = round(props(biggest).Centroid(1));
% % if (props(biggest).Area>0.6*m*n) && length(props)>1
% %     sbiggest = idx(2);
% %     seedRow = round(props(sbiggest).Centroid(2));
% %     seedCol = round(props(sbiggest).Centroid(1));
% % end
% 
% % using red component of an image
% segIm = zeros(size(imPrep(:,:,1)));
% Trg = 0.6*std(imPrep(:,:,1),[],'all');
% while sum(segIm == 1)< 0.00005*m*n
% % segIm = grayconnected(imPrep(:,:,1),seedRow,seedCol,Trg);
% segIm = regiongrowing(imPrep(:,:,1),seedRow,seedCol,Trg);
% Trg = 1.25*Trg;
% end

% % using grayscale image
% segIm = zeros(size(imPrepGray));
% Trg = 0.61*std(imPrepGray,[],'all');
% while sum(segIm == 1)< 0.00005*m*n
% segIm = grayconnected(imPrepGray,seedRow,seedCol,Trg);
% Trg = 1.25*Trg;
% end
% segIm = grayconnected(imPrepLab(:,:,2),seedRow,seedCol,0.5*std(imPrepGray,[],'all'));

% binaryMap = imfill(segIm,'holes');


%% Hough transform for circles

% stdPic = stdfilt(imPrep(:,:,1),true(5));
% Ts = graythresh(stdPic);
% Tv = graythresh(imPrepHSV(:,:,3)); % elimination of edges in dark background
% imEdge = (stdPic>Ts & imPrepHSV(:,:,3)>Tv);
% % imEdge = edge(imPrepGray,'canny',[.03 .1],sqrt(2)); % constants set according to Sanchez2018
% rs = 5:2:100; % range of radia
% HS = zeros(size(imPrep,1),size(imPrep,2),length(rs));
% r_ind = 1;
% [X,Y] = find(imEdge == 1); % finding edges
% % filling the Hough space
% for r = rs
%     tmp_c = gen_circle(r);
%     for i = 1:length(X)
%         c1 = X(i);
%         c2 = Y(i);
%         if c1 > r && c1< (size(imPrep,1) - r)
%             if c2 > r && c2< (size(imPrep,2) - r)
%                 HS((c1-r):(c1+r),(c2-r):(c2+r),r_ind) = HS((c1-r):(c1+r),(c2-r):(c2+r),r_ind)+tmp_c;
%             end
%         end
%     end
%     r_ind = r_ind + 1;
% end
% % finding the center of the most probable circle in edge representation
% [linInd] = find(HS == max(HS,[],'all'),1,'first');
% [y,x,r] = ind2sub(size(HS),linInd); 
% 
% if length(x)>1 || length(y)>1
%     x = floor(mean(x));
%     y = floor(mean(y));
% end
% % 
%% region growing
% segIm = zeros(size(imPrep(:,:,1)));
% Trg = 1.2*std(imPrep(:,:,1),[],'all');
% while sum(segIm == 1)< 0.00005*m*n
% segIm = grayconnected(imPrep(:,:,1),y,x,Trg);
% % segIm = regiongrowing(imPrep(:,:,1),seedRow,seedCol,Trg);
% Trg = 1.25*Trg;
% end
% binaryMap = imfill(segIm,'holes');

% for i = 1:o
% segIm(:,:,i) = grayconnected(imPrep(:,:,i),y,x,0.1*std(imPrepGray,[],'all')); % position is defined by Hough t.
% end
% sumRegion = reshape(sum(sum(segIm)),[3 1 1]);
% [~,smallObjChannel] = min(sumRegion);
% binaryMap = imfill(segIm(:,:,smallObjChannel),'holes');
%% geometric contours
% % % finding the smallest object in 3 results of region growing
% % sumRegion = reshape(sum(sum(segIm)),[3 1 1]);
% % [~,smallObjChannel] = min(sumRegion);
% % % level sets
% % binaryMap = activecontour(rgb2gray(imCropped),imdilate(segIm(:,:,smallObjChannel),[1 1 1; 1 1 1; 1 1 1]));

%% Method with hysteresis thresholding and Region Growing

% inputImage=FClear(inputImage,bEdgeMask);
% imPrep = FLight(inputImage);
% [x,y]  = FHysThres(imPrep);
% binaryMap = FRegionGrow(imPrep,x,y);
%% Method with Hough Transform and Region Growing

inputImage=FClear(inputImage,bEdgeMask);
% imPrep = inputImage;
imPrep = FLight(inputImage);
[x,y]  = FHouTrans(imPrep);
binaryMap = FRegionGrow(imPrep,x,y);


end

