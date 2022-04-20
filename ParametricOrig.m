%% Parametric Field for original images
% @JanSima,@OndrejNantl,@TerezieDobrovolna
clear all; clc;
%% loading
% Change the pathway to the dataset!
pathCVC_Orig ='D:\HONZA\Honza VUT\Ing\SEMESTR2\ABO\Projekt\polypy\CVC-ClinicDB\CVC-ClinicDB\Original\';
% pathCVC_Mask = 'D:\HONZA\Honza VUT\Ing\SEMESTR2\ABO\Projekt\polypy\CVC-ClinicDB\CVC-ClinicDB\Ground Truth\';
a = dir([pathCVC_Orig '*.tif']);
n = numel(a);
ParametricFieldOrig = [];
% Gray - grayscale image
% R,G,B - color channels
% Mean, SD, E - mean, standard deviation and entropy
Labels = {'MeanGray','MeanR','MeanG','MeanB','SDGray','SDR','SDG',...
    'SDB','EGray','ER','EG','EB'};
%% calculation of selected parameters for every image in dataset
for idx = 1:n
    im = rgb2gray(im2double(imread([pathCVC_Orig, num2str(idx) '.tif'])));
    imColor = im2double(imread([pathCVC_Orig, num2str(idx) '.tif']));
%     mask = im2double(imread([pathCVC_Mask, num2str(idx) '.tif']));

    % removal of black edge
    clear bEdgeMask bEdgeMask2 bEdgeMask3 imCropped imCroppedRow %maskCropped maskCroppedRow
    imHSV = rgb2hsv(imColor); % transfer into HSV color space
    bEdgeMask = (imHSV(:,:,3) <= 0.2); % obtaining mask of black edge
    newRowCount = 0;
    % cropping the rows which are only dark 
    for i = 1:size(bEdgeMask,1)
        if any(bEdgeMask(i,:) ~= 1)
            newRowCount = newRowCount + 1;
            bEdgeMask2(newRowCount,:) = bEdgeMask(i,:);
            imCroppedRow(newRowCount,:,:) = image(i,:,:);
%             GTCroppedRow(newRowCount,:) = GT(i,:);
        end
    end
    newColCount = 0;
    % cropping the colums which are only dark 
    for j = 1:size(bEdgeMask2,2)
        if any(bEdgeMask2(:,j) ~= 1)
            newColCount = newColCount + 1;
            bEdgeMask3(:,newColCount) = bEdgeMask2(:,j);
            imCropped(:,newColCount,:) = imCroppedRow(:,j,:);
%             GTCropped(:,newColCount) = GTCroppedRow(:,j);
        end
    end
    % transfer into grayscale
    imGray = rgb2gray(im2double(imCropped));

    % Calculating mean, SD and entropy in RGB and grayscale image +
    % difference
    % Mean imGray
    ParametricFieldOrig(idx,1) = mean(mean(imGray));
    % a) R
    Color = imCropped(:,:,1);
    ParametricFieldOrig(idx,2) = mean(mean(Color));
    % b) G
    Color = imCropped(:,:,2);
    ParametricFieldOrig(idx,3) = mean(mean(Color));
    % b) B
    Color = imCropped(:,:,3);
    ParametricFieldOrig(idx,4) = mean(mean(Color));

    % SD imGray
    ParametricFieldOrig(idx,5) = std(std(imGray));
    % SD imColor 
    % a) R
    Color = imCropped(:,:,1);
    ParametricFieldOrig(idx,6) = std(std(Color));
    % b) G
    Color = imCropped(:,:,2);
    ParametricFieldOrig(idx,7) = std(std(Color));
    % b) B
    Color = imCropped(:,:,3);
    ParametricFieldOrig(idx,8) = std(std(Color));
    
    % Entropy
    ParametricFieldOrig(idx,9) = entropy(imGray);
    % R
    Color = imCropped(:,:,1);
    ParametricFieldOrig(idx,10) = entropy(Color);
    % G
    Color = imCropped(:,:,2);
    ParametricFieldOrig(idx,11) = entropy(Color);
    % B
    Color = imCropped(:,:,3);
    ParametricFieldOrig(idx,12) = entropy(Color);
end

%% saving the results
save("Parametric_Field_Original","ParametricFieldOrig","Labels")