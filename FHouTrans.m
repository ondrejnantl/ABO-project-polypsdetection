function [x,y,r,HTMax] = FHouTrans(imPrep)
% This function performs estimation of seed coordinates for region growing
% from preprocessed RGB colonoscopy image using Hough transform for circles
% -------------------------------------------------------------------------
% Input: 
% imPrep - input preprocessed RGB image obtained during colonoscopy (after 
% cropping in evaluation function, removal of specular highlights and 
% adjustment of lighting)
%
% Output:
% x - x coordinate of center of the most likely present circle in the image 
% (the seed)
% 
% y - y coordinate of center of the most likely present circle in the image
% (the seed)
% 
% r - diameter of the most likely present circle in the image
% -------------------------------------------------------------------------
% Authors: Terezie Dobrovolná, Ondřej Nantl, Jan Šíma
% =========================================================================
imPrepHSV = rgb2hsv(imPrep);
% variant using local standard deviation thresholding
stdPic = stdfilt(imPrep(:,:,1),true(5));
Ts = graythresh(stdPic); % Otsu method
Tv = graythresh(imPrepHSV(:,:,3)); % elimination of edges in dark background
imEdge = (stdPic>Ts & imPrepHSV(:,:,3)>Tv);

rs = 10:2:100; % range of radia
HS = zeros(size(imPrep,1),size(imPrep,2),length(rs));
r_ind = 1;
[X,Y] = find(imEdge == 1); % finding edges
% filling the Hough space
for r = rs
    tmp_c = gen_circle(r); % generate circle of apropriate radius
    for i = 1:length(X)
        c1 = X(i);
        c2 = Y(i);
        if c1 > r && c1< (size(imPrep,1) - r)
            if c2 > r && c2< (size(imPrep,2) - r)
                % adding circle into Hough space
                HS((c1-r):(c1+r),(c2-r):(c2+r),r_ind) = HS((c1-r):(c1+r),(c2-r):(c2+r),r_ind)+tmp_c;
            end
        end
    end
    r_ind = r_ind + 1;
end
% finding the center of the most probable circle in edge representation
% as maxima of Hough space
HTMax = max(HS,[],'all');
[linInd] = find(HS == HTMax,1,'last');
[y,x,r] = ind2sub(size(HS),linInd); 

end