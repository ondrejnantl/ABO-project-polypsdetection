function [x,y]  = FHysThres(imPrep)
imPre = rgb2gray(im2double(imPrep));
thresholds = multithresh(imPre,2);
T1 = thresholds(2);
T2 = thresholds(1);
h=hyst_thresh(imPre,T1,T2);

centers = regionprops(h);
loc = find(max([centers.Area]));
positions = centers(loc).Centroid;
x = positions(1);
y = positions(2);
% figure
% imshow(imPrep.*~h,[])
% hold on
% plot(x,y,'r+','LineWidth',20)

end