%% ABO - Projekt č.10 Polypy
% @JanSima
clear all; clc;
%% načtení
% Změň si cestu k souboru!

pathCVC_Orig = 'D:\HONZA\Honza VUT\Ing\SEMESTR2\ABO\Projekt\polypy\CVC-ClinicDB\CVC-ClinicDB\Original\';
pathCVC_Mask = 'D:\HONZA\Honza VUT\Ing\SEMESTR2\ABO\Projekt\polypy\CVC-ClinicDB\CVC-ClinicDB\Ground Truth\';
for idx = 99
    im = rgb2gray(im2double(imread([pathCVC_Orig, num2str(idx) '.tif'])));
    mask = im2double(imread([pathCVC_Mask, num2str(idx) '.tif']));

    %     % Vykreslení polypu
    %     figure
    %     imshow(im.*mask,[])

    % Smazání odlesků
    mira_potlaceni = 1.5;
    im_odlesk = im>(mean(im(:))+mira_potlaceni*std(im(:))); % Práh pro odstranění odlesků
    image = regionfill(im,im_odlesk);

    %     % Vykreslení smazání odlesků
    %     figure
    %     subplot 121
    %     imshow(im_odlesk)
    %     title('Odelsky')
    %     subplot 122
    %     imshow(im_bezodlesk)
    %     title('Bez odlesků')

    % Vykreslení upraveného obrazu + histogram
    figure
    subplot 121
    imshow(image,[])
    title('Original')
    subplot 122
    imhist(image)
    title('Histogram')

    % Výpočet gradientu (1. derivace)
    grad = imgradient(image);

    figure
    subplot 121
    imshow(grad,[])
    title('Gradient')
    subplot 122
    imhist(grad)
    title('Histogram')

%     % Funkce Edge
%     figure
%     E = edge(image,'sobel');
%     imshow(E)

    % Pokus o filtraci 
    filt = [-1 -1 -1;0 0 0;-1 -1 -1];
    edgeIm = filter2(filt,image);
    figure
    subplot 121
    imshow(edgeIm,[])
    subplot 122
    imshow(edge(edgeIm,'zerocross'),[])
        
    K = 9.2;
    h1 = (1/(K-8)).*[-1,-1,-1;-1,K,-1;-1,-1,-1];
    sh1 = conv2(image,h1,'same');
    figure
    subplot 121
    imshow(sh1,[])
    subplot 122
    imshow(edge(sh1,'approxcanny'),[])
end

%% detekce v zóně, kde je náznak elipsovitého útvaru
