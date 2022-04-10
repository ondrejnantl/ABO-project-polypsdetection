%% ABO - Projekt c.10 Polypy
% @JanSima,@OndrejNantl
clear all; clc;
%% nacteni
% Zmen si cestu k souboru!
pathCVC_Orig ='D:\HONZA\Honza VUT\Ing\SEMESTR2\ABO\Projekt\polypy\CVC-ClinicDB\CVC-ClinicDB\Original\';
pathCVC_Mask = 'D:\HONZA\Honza VUT\Ing\SEMESTR2\ABO\Projekt\polypy\CVC-ClinicDB\CVC-ClinicDB\Ground Truth\';
a = dir([pathCVC_Orig '*.tif']);
n = numel(a);
ParametricField = [];
% Gray - šedotónovy obraz
% RGB - barevné kanály
% BG - background
% - a / - rozdíl a poměr mezi polypem a pozadím
% Mean, SD, E - průměr, směrodatná odchylka a entopie
Labels = {'MeanGray','MeanGrayBG','MeanGray-','MeanGray/','MeanR',...
    'MeanRBG','MeanR-','MeanR/','MeanG','MeanGBG','MeanG-','MeanG/',...
    'MeanB','MeanBBG','MeanB-','MeanB/','SDGray','SDGrayBG','SDGray-',...
    'SDGray/','SDR','SDRBG','SDR-','SDR/','SDG','SDGBG','SDG-','SDG/',...
    'SDB','SDBBG','SDB-','SDB/','EGray','EGrayBG','EGray-','EGray/',...
    'ER','ERBG','ER-','ER/','EG','EGBG','EG-','EG/','EB','EBBG','EB-','EB/'};
%%
for idx = 1:n
    im = rgb2gray(im2double(imread([pathCVC_Orig, num2str(idx) '.tif'])));
    imColor = im2double(imread([pathCVC_Orig, num2str(idx) '.tif']));
    mask = im2double(imread([pathCVC_Mask, num2str(idx) '.tif']));

    % odstraneni ramecku
    clear bEdgeMask bEdgeMask2 bEdgeMask3 imCropped imCroppedRow maskCropped maskCroppedRow
    imHSV = rgb2hsv(imColor); % prevod do HSV
    bEdgeMask = (imHSV(:,:,3) <= 0.2); % konstanta podle Sanchez2018
    newRowCount = 0;
    for i = 1:size(bEdgeMask,1)
        if any(bEdgeMask(i,:) ~= 1)
            newRowCount = newRowCount + 1;
            bEdgeMask2(newRowCount,:) = bEdgeMask(i,:);
            imCroppedRow(newRowCount,:,:) = imColor(i,:,:);
            maskCroppedRow(newRowCount,:) = mask(i,:);
        end
    end
    newColCount = 0;
    for j = 1:size(bEdgeMask2,2)
        if any(bEdgeMask2(:,j) ~= 1)
            newColCount = newColCount + 1;
            bEdgeMask3(:,newColCount) = bEdgeMask2(:,j);
            imCropped(:,newColCount,:) = imCroppedRow(:,j,:);
            maskCropped(:,newColCount) = maskCroppedRow(:,j);
        end
    end
    % převod na gray
    imGray = rgb2gray(im2double(imCropped));

    % Výpočet mean, median a SD v RGB a v šedotónovém obraze + rozdíl
    % Mean imGray
    ParametricField(idx,1) = mean(mean(imGray(maskCropped==1)));
    ParametricField(idx,2) = mean(mean(imGray(maskCropped==0)));
    ParametricField(idx,3) = ParametricField(idx,2)-ParametricField(idx,1);
    ParametricField(idx,4) = ParametricField(idx,2)/ParametricField(idx,1);
    % Mean imColor 
    % a) R
    Color = imCropped(:,:,1);
    ParametricField(idx,5) = mean(mean(Color(maskCropped==1)));
    ParametricField(idx,6) = mean(mean(Color(maskCropped==0)));
    ParametricField(idx,7) = ParametricField(idx,6)-ParametricField(idx,5);
    ParametricField(idx,8) = ParametricField(idx,6)/ParametricField(idx,5);
    % b) G
    Color = imCropped(:,:,2);
    ParametricField(idx,9) = mean(mean(Color(maskCropped==1)));
    ParametricField(idx,10) = mean(mean(Color(maskCropped==0)));
    ParametricField(idx,11) = ParametricField(idx,10)-ParametricField(idx,9);
    ParametricField(idx,12) = ParametricField(idx,10)/ParametricField(idx,9);
    % b) B
    Color = imCropped(:,:,3);
    ParametricField(idx,13) = mean(mean(Color(maskCropped==1)));
    ParametricField(idx,14) = mean(mean(Color(maskCropped==0)));
    ParametricField(idx,15) = ParametricField(idx,14)-ParametricField(idx,13);
    ParametricField(idx,16) = ParametricField(idx,14)/ParametricField(idx,13);

    % SD imGray
    ParametricField(idx,17) = std(imGray(maskCropped==1));
    ParametricField(idx,18) = std(imGray(maskCropped==0));
    ParametricField(idx,19) = ParametricField(idx,18)-ParametricField(idx,17);
    ParametricField(idx,20) = ParametricField(idx,18)/ParametricField(idx,17);
    % SD imColor 
    % a) R
    Color = imCropped(:,:,1);
    ParametricField(idx,21) = std(Color(maskCropped==1));
    ParametricField(idx,22) = std(Color(maskCropped==0));
    ParametricField(idx,23) = ParametricField(idx,22)-ParametricField(idx,21);
    ParametricField(idx,24) = ParametricField(idx,22)/ParametricField(idx,21);
    % b) G
    Color = imCropped(:,:,2);
    ParametricField(idx,25) = std(Color(maskCropped==1));
    ParametricField(idx,26) = std(Color(maskCropped==0));
    ParametricField(idx,27) = ParametricField(idx,26)-ParametricField(idx,25);
    ParametricField(idx,28) = ParametricField(idx,26)/ParametricField(idx,25);
    % b) B
    Color = imCropped(:,:,3);
    ParametricField(idx,29) = std(Color(maskCropped==1));
    ParametricField(idx,30) = std(Color(maskCropped==0));
    ParametricField(idx,31) = ParametricField(idx,30)-ParametricField(idx,29);
    ParametricField(idx,32) = ParametricField(idx,30)/ParametricField(idx,29);
    
    % Entropy
    ParametricField(idx,33) = entropy(imGray(maskCropped==1));
    ParametricField(idx,34) = entropy(imGray(maskCropped==0));
    ParametricField(idx,35) = ParametricField(idx,34)-ParametricField(idx,33);
    ParametricField(idx,36) = ParametricField(idx,34)/ParametricField(idx,33);
    % R
    Color = imCropped(:,:,1);
    ParametricField(idx,37) = entropy(Color(maskCropped==1));
    ParametricField(idx,38) = entropy(Color(maskCropped==0));
    ParametricField(idx,39) = ParametricField(idx,38)-ParametricField(idx,37);
    ParametricField(idx,40) = ParametricField(idx,38)/ParametricField(idx,37);
    % G
    Color = imCropped(:,:,2);
    ParametricField(idx,41) = entropy(Color(maskCropped==1));
    ParametricField(idx,42) = entropy(Color(maskCropped==0));
    ParametricField(idx,43) = ParametricField(idx,42)-ParametricField(idx,41);
    ParametricField(idx,44) = ParametricField(idx,42)/ParametricField(idx,41);
    % B
    Color = imCropped(:,:,3);
    ParametricField(idx,45) = entropy(Color(maskCropped==1));
    ParametricField(idx,46) = entropy(Color(maskCropped==0));
    ParametricField(idx,47) = ParametricField(idx,46)-ParametricField(idx,45);
    ParametricField(idx,48) = ParametricField(idx,46)/ParametricField(idx,45);
end

%%
save("Parametric_Field","ParametricField","Labels")