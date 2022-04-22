%% Script for training the classficator
clear all,clc,close all
%% loading names of sorted images, preparing training and test datasets
pathCVC_Orig = 'D:\andyn\OneDrive - Vysoké učení technické v Brně\materialy_4r_moje\MPA-ABO\projekt\CVC-ClinicDB\Original\';
pathCVC_Mask = 'D:\andyn\OneDrive - Vysoké učení technické v Brně\materialy_4r_moje\MPA-ABO\projekt\CVC-ClinicDB\Ground Truth\';

filenames = table2array(readtable('2Skupiny.xlsx'));

permIDs1 = randperm(size(filenames,1));
permIDs2 = randperm(size(filenames,1));

trainIDs = [permIDs1(1:0.6*size(filenames,1)) permIDs2(1:0.6*size(filenames,1))];
trainCats = [ones(1,0.6*size(filenames,1)) repmat(2,1,0.6*size(filenames,1))];
testIDs = [permIDs1((0.6*size(filenames,1)+1):end) permIDs2((0.6*size(filenames,1)+1):end)];
testCats = [ones(1,0.4*size(filenames,1)) repmat(2,1,0.4*size(filenames,1))];
allIDs = [trainIDs testIDs];
allCats = [trainCats testCats];
%% obtaining parameters
featureSpace = zeros(length(allIDs),8);
catNames = {'Oval','Unclassifiable'};
featureNames = {'HysThMean','HysThArea','RGMean1','RGArea1','HTRadius','HTMax','RGMean2','RGArea2'};
diceCoef = zeros(length(allIDs),2);
IoU = zeros(length(allIDs),2);
for imIter = 1:length(allIDs)
    image = im2double(imread([pathCVC_Orig, num2str(allIDs(imIter)), '.tif']));
    GT = im2double(imread([pathCVC_Mask, num2str(allIDs(imIter)) '.tif']));
    GT(GT<1) = 0;
    % cropping of the black frame
    clear bEdgeMask bEdgeMask2 bEdgeMask3 imCropped imCroppedRow GTCropped GTCroppedRow
    imHSV = rgb2hsv(image); % transfer into HSV color space
    bEdgeMask = (imHSV(:,:,3) <= 0.2); % obtaining mask of black edge
    newRowCount = 0;
    % cropping the rows which are only dark 
    for i = 1:size(bEdgeMask,1)
        if any(bEdgeMask(i,:) ~= 1)
            newRowCount = newRowCount + 1;
            bEdgeMask2(newRowCount,:) = bEdgeMask(i,:);
            imCroppedRow(newRowCount,:,:) = image(i,:,:);
            GTCroppedRow(newRowCount,:) = GT(i,:);
        end
    end
    newColCount = 0;
    % cropping the colums which are only dark 
    for j = 1:size(bEdgeMask2,2)
        if any(bEdgeMask2(:,j) ~= 1)
            newColCount = newColCount + 1;
            bEdgeMask3(:,newColCount) = bEdgeMask2(:,j);
            imCropped(:,newColCount,:) = imCroppedRow(:,j,:);
            GTCropped(:,newColCount) = GTCroppedRow(:,j);
        end
    end
    % Method with hysteresis thresholding and region growing
    imCropped=FClear(imCropped,bEdgeMask3);
    imPrep = FLight(imCropped);
    [x,y,featureSpace(imIter,1),featureSpace(imIter,2)]  = FHysThres(imPrep);
    [binaryMap,featureSpace(imIter,3),featureSpace(imIter,4)] = FRegionGrow(imPrep,x,y);
    % Method with Hough transform and region growing
    [x,y,featureSpace(imIter,5),featureSpace(imIter,6)]  = FHouTrans(imPrep);
    [binaryMap2,featureSpace(imIter,7),featureSpace(imIter,8)] = FRegionGrow(imPrep,x,y);
    
    diceCoef(imIter,1) = dice(binaryMap,logical(GTCropped));
    IoU(imIter,1) = jaccard(binaryMap,logical(GTCropped));
    diceCoef(imIter,2) = dice(binaryMap2,logical(GTCropped));
    IoU(imIter,2) = jaccard(binaryMap2,logical(GTCropped));
end
trainFeatureSpace = featureSpace(1:(0.6*size(featureSpace,1)),:);
testFeatureSpace = featureSpace((0.6*size(featureSpace,1)+1):end,:);

% creating boxplots - for all features according to their category in 
% manual sorting
figure
for i = 1:size(trainFeatureSpace,2)
    subplot(4,2,i)
    boxplot(trainFeatureSpace(:,i),trainCats,'Labels',catNames);
    title(featureNames{i})
end

%% classificator
% training
Mdl = TreeBagger(100,trainFeatureSpace,trainCats,'Method','classification');
% predicting and evaluating on test data
predCats = str2num(cell2mat(predict(Mdl,testFeatureSpace)));
confusionchart(testCats',predCats) % 1- oval, 2 - unclassifiable
