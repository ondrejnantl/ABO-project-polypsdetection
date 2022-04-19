function imCropped = FClear(inputImage,bEdgeMask)
% smazani odlesku - alternativa 2
% clear bEdgeMask bEdgeMask2 bEdgeMask3 imCropped imCroppedRow GTCropped GTCroppedRow
pm = rangefilt(rgb2gray(inputImage),true(5));
T = graythresh(pm);
reflMask = imbinarize(imfill(pm,'holes'),T);
imCropped = inpaintCoherent(inputImage,logical((~bEdgeMask).*reflMask),'SmoothingFactor',5,'Radius',5);
end