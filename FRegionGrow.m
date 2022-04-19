function [final] = FRegionGrow(imPrep,x,y)
%% region growing
% seedRow = randi(size(imPrep,1)); % je nutne vymyslet jak zjistit pozici seminka
% seedCol = randi(size(imPrep,2));
% [~,seedRow] = max(sum(imPrep(:,:,1),2)); % je nutne vymyslet jak zjistit pozici seminka
% [~,seedCol] = max(sum(imPrep(:,:,1),1));
[~,~,o] = size(imPrep);
% figure;
for i = 1:o
segIm(:,:,i) = grayconnected(imPrep(:,:,i),round(x),round(y),0.02); 
% segIm(:,:,i) = grayconnected(imPrep(:,:,i),y,x,0.02); % pozici definuje Houghova t.
% segIm(:,:,i) = regiongrowing(imPrep(:,:,i),y,x,0.02); % jina region growing funkce, pozici definuje Houghova t.
% subplot(1,3,i)
% imshow(segIm(:,:,i));hold on; plot(x,y,'rx')
end
sumRegion = reshape(sum(sum(segIm)),[3 1 1]);
[~,smallObjChannel] = min(sumRegion);
final = imfill(segIm(:,:,smallObjChannel),'holes');
% figure
% imshowpair(maskCropped,final)
end