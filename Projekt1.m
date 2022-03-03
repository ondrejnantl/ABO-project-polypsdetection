%% ABO - Projekt č.10 Polypy
% @JanSima
clear all; clc;
%% načtení
% Změň si cestu k souboru!

pathCVC_Orig = 'D:\HONZA\Honza VUT\Ing\SEMESTR2\ABO\Projekt\polypy\CVC-ClinicDB\CVC-ClinicDB\Original\';
pathCVC_Mask = 'D:\HONZA\Honza VUT\Ing\SEMESTR2\ABO\Projekt\polypy\CVC-ClinicDB\CVC-ClinicDB\Ground Truth\';
for idx = 333
    im = rgb2gray(im2double(imread([pathCVC_Orig, num2str(idx) '.tif'])));
    mask = im2double(imread([pathCVC_Mask, num2str(idx) '.tif']));

    %     % Vykreslení polypu
    %     figure
    %     imshow(im.*mask,[])

    % Smazání odlesků
    mira_potlaceni = 1.5;
    im_odlesk = im>(mean(im(:))+mira_potlaceni*std(im(:))); % Práh pro odstranění odlesků
    image = regionfill(im,im_odlesk);
    image2 = imfill(ordfilt2(im,1,ones(5,5)),'holes');

    %     image = image1-image2;

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

    % Lawovi filtry 
    load("Laws.mat")
    L = 9;
    [m,n] = size(image);
    pm = zeros(m,n,L);

    K = 35;
    mask = ones(K,K)./(K^2);
    MASK = fft2(mask,m,n);
    figure
    for i=1:L
        out = abs(conv2(image, rot90(law(:,:,i),2), "same"));
        OUT = fft2(out) .* MASK;
        out = real(ifft2(OUT));
        pm(:,:,i) = out;
        subplot(3,3,i)
        imshow(out,[])
    end

    % Homomorfická filtrace na srovnání jasu v obraze
    figure
    subplot(121)
    imshow(image);
    homofil(image,10,m,n,2);

    %     % Funkce Edge
    %     figure
    %     E = edge(image,'sobel');
    %     imshow(E)

end

%% detekce v zóně, kde je náznak elipsovitého útvaru
%% 3d vlnka a použít korelaci