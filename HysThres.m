function [x,y]  = HysThres(imPrep)
imPre = rgb2gray(im2double(imPrep));
thresholds = multithresh(imPre,2);
T1 = thresholds(2);
T2 = thresholds(1);
h=hyst_thresh(imPre,T1,T2);

centers = regionprops(h);
loc = find(max([centers.Area]));
[x,y] = centers(loc).Centroids;
end