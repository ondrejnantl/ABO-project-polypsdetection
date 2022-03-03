function [resultDataMatrix,diceCoef,IoU] = polypsEval(datasetPath)
% This function can be used for evaluating the performance of polyp
% detection/segmentation algorithm
% 
% 
% =========================================================================
% Gaining the names of images  - original and ground truth
imDS = imageDatastore([datasetPath '\Original'],'ReadFcn', @read2gray);
groundTruthDS = ls([datasetPath '\Ground Truth']);

% Obtaining information about dimensions of input images and their count
imForDim = size(imread(imDS.Files{1}));
numImages = size(imDS.Files,1);

resultDataMatrix = double(zeros(imForDim(1),imForDim(2),numImages));
diceCoef = zeros(numImages,1);
IoU = zeros(numImages,1);

for imIter = 1:numImages
    % loading one image and its ground truth mask 
    image = readimage(imDS,imIter);
    GT = readimage(groundTruthDS,imIter);
    % analysis of the image using our algorithm
    % this is important     
    resultDataMatrix(:,:,imIter) = detectPolyps(image);
    % evaluation of our algorithm using Dice and Jaccard coefficients
    diceCoef(imIter) = dice(resultDataMatrix(:,:,imIter),GT);
    IoU(imIter) = jaccard(resultDataMatrix(:,:,imIter),GT);
end

end

function image = read2gray(path)
    image = rgb2gray(im2double(imread(path)));
end
