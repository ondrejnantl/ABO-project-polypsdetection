function [x,y,r] = FHouTtrans(imPrep)
%% Houghova transforHouTtransace pro kruh - mohla by fungovat
imEdge = edge(rgb2gray(imPrep),'canny',[.03 .1],sqrt(2)); % varianta s rgb
% imEdge = edge(imPrepLab(:,:,3),'canny'); % varianta s Lab

rs = 5:50;
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
% imshow5(HS)
[linInd] = find(HS == max(HS,[],'all'));
[y,x,r] = ind2sub(size(HS),linInd);

imshow(imPrep);hold on; for i = 1:length(x);h = images.roi.Circle(gca,'Center',[x(i) y(i)],'Radius',r(i));end

% if length(x)>1 || length(y)>1
%     x = floor(mean(x));
%     y = floor(mean(y));
% end
% roiMask = createMask(h);
end