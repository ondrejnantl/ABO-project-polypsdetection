%% ABO - Projekt c.10 Polypy
% @JanSima,@OndrejNantl
clear all; clc;
%% nacteni
% Zmen si cestu k souboru!
pathCVC_Orig = 'D:\andyn\OneDrive - Vysoké učení technické v Brně\materialy_4r_moje\MPA-ABO\projekt\CVC-ClinicDB\Original\';
pathCVC_Mask = 'D:\andyn\OneDrive - Vysoké učení technické v Brně\materialy_4r_moje\MPA-ABO\projekt\CVC-ClinicDB\Ground Truth\';
for idx = 199
    im = rgb2gray(im2double(imread([pathCVC_Orig, num2str(idx) '.tif'])));
    imColor = im2double(imread([pathCVC_Orig, num2str(idx) '.tif']));
    mask = im2double(imread([pathCVC_Mask, num2str(idx) '.tif']));
end

%% odstraneni ramecku, smazani odlesku a korekce osvetleni

% odstraneni ramecku
clear bEdgeMask bEdgeMask2 bEdgeMask3 imCropped imCroppedRow maskCropped maskCroppedRow gaussFilt
imHSV = rgb2hsv(imColor); % prevod do HSV
bEdgeMask = (imHSV(:,:,3) <= 0.2); % konstanta podle Sanchez2018
% imCropped = zeros(size(imHSV));
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

% %!! nove - rozmazani do okraju
% N = 39;
% meanMask = 1/(N^2).*ones(N,N);
% blurMask = imdilate(bEdgeMask3,[strel('line',2,45) strel('line',2,135) strel('line',1,0) strel('line',1,90)]);
% 
% meanFilt = imfilter(imCropped,meanMask);
% imCropped = imCropped + 1.5.*blurMask.*meanFilt;
% imshow(imCropped)

% smazani odlesku - alternativa 2
pm = rangefilt(rgb2gray(imCropped),true(5));
T = graythresh(pm);
reflMask = imbinarize(imfill(pm,'holes'),T);
imCropped = inpaintCoherent(imCropped,logical((~bEdgeMask3).*reflMask),'SmoothingFactor',5,'Radius',5);
% imshow(imCropped)
% imCropped = inpaintExemplar(imCropped,bEdgeMask3,'PatchSize',[30 30]);
% imCropped = imgaussfilt(imCropped,0.8);
% figure
% imshow(imCropped)

% uprava osvetleni
[m,n,o] = size(imCropped);
mm = zeros(m,n,o);
N = 20;
meanMask = 1/(N^2).*ones(N,N);
for j = 1:o
    mm(:,:,j) = 0.3.*conv2(imCropped(:,:,j),meanMask,'same'); % vaha 0.3 podle Sanchez2018 - lehce zmenena
end
imPrep = imCropped - mm;
imPrepGray = rgb2gray(imPrep);
imPrepHSV = rgb2hsv(imPrep);
imPrepLab = rgb2lab(imPrep);
figure
imshowpair(imPrep,maskCropped,'montage')

%% uprava preprocesovaneho obrazu - zatim nespoustet
figure;
for j = 1:3
    edgedImage(:,:,j) = edge(imPrep(:,:,j),'canny',[.03 .1],sqrt(2)); % konstanty nastaveny podle Sanchez2018
    subplot(1,3,j)
    imshow(edgedImage(:,:,j))
end
% edgedImage = edgedImage(:,:,1).*edgedImage(:,:,2).*edgedImage(:,:,3);

%% Houghova transformace pro kruh - mohla by fungovat
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
%% region growing
seedRow = 140; % je nutne vymyslet jak zjistit pozici seminka
seedCol = 69;
% [~,seedRow] = max(sum(imPrep(:,:,1),2)); % je nutne vymyslet jak zjistit pozici seminka
% [~,seedCol] = max(sum(imPrep(:,:,1),1));
figure;
for i = 1:o
% segIm(:,:,i) = grayconnected(imPrep(:,:,i),seedRow,seedCol,0.02); % pozici definujeme my
segIm(:,:,i) = grayconnected(imPrep(:,:,i),y,x,0.02); % pozici definuje Houghova t.
% segIm(:,:,i) = regiongrowing(imPrep(:,:,i),y,x,0.02); % jina region growing funkce, pozici definuje Houghova t.
subplot(1,3,i)
imshow(segIm(:,:,i));hold on; plot(x,y,'rx')
end
sumRegion = reshape(sum(sum(segIm)),[3 1 1]);
[~,smallObjChannel] = min(sumRegion);
final = imfill(segIm(:,:,smallObjChannel),'holes');
figure
imshowpair(maskCropped,final)

%% hysterezni prahovani
% BW = hysthresh(rgb2gray(imPrep),0.6,0.4);
BW = hysthresh(imPrepLab(:,:,2),0.75*max(imPrepLab(:,:,2),[],'all'),0.65*max(imPrepLab(:,:,2),[],'all'));
BW = imerode(BW,strel('disk',3));
imshow(BW)
props = regionprops(BW,'Area','Centroid','Circularity','ConvexHull','ConvexImage','FilledImage','MajorAxisLength','MinorAxisLength');
[~,idx] = sort([props.Area]);
biggest = idx(1) ;
seedRow = round(props(biggest).Centroid(2));
seedCol = round(props(biggest).Centroid(1));
if (props(biggest).Area>0.4*m*n) && length(props)>1
    sbiggest = idx(2);
    seedRow = round(props(sbiggest).Centroid(2));
    seedCol = round(props(sbiggest).Centroid(1));
end
% segIm = grayconnected(imPrepGray,seedRow,seedCol,0.1*std(imPrepGray,[],'all'));
segIm = grayconnected(imPrepLab(:,:,2),seedRow,seedCol,0.5*std(imPrepLab(:,:,2),[],'all'));
segIm = imfill(segIm,'holes');
imshowpair(maskCropped,segIm);hold on; plot(seedCol,seedRow,'rx')
%% mikrostrukturni analyza, clustering a nasledna analyza s geometrickou konturou
load('Laws.mat')
imPrepV = imPrepLab(:,:,2);
% creating parametric maps
L = size(law,3);
pm = zeros(m,n,L);

for i = 1:L
    pm(:,:,i) = abs(conv2(imPrepGray,rot90(law(:,:,i),2),'same')); % rotation of mask must be done
end
 
% k-means 
% clus = kmeans(reshape(pm,m*n,L),3);
% clusim = reshape(clus,m,n); 
% c-means
[~,clus] = fcm(reshape(imPrepV,m*n,1),3,[2 100 1e-4 0]);
% [~,clus] = fcm(reshape(pm,m*n,L),3,[2 100 1e-4 0]);
clusim = reshape(clus',m,n,3);
imshow(clusim)
% creating threshold
thresArray = zeros(o,3);
for i = 1:o
    clusim(:,:,i) = medfilt2(clusim(:,:,i),[3 3]);
    % for k-means
%     thresArray(i,1) = min(imPrepGray(clusim(:,:)==i));
%     thresArray(i,2) = max(imPrepGray(clusim(:,:)==i));
    % pro c-means
    thresArray(i,1) = min(imPrepGray(clusim(:,:,i)>0.75));
    thresArray(i,2) = median(imPrepGray(clusim(:,:,i)>0.75));
    thresArray(i,3) = max(imPrepGray(clusim(:,:,i)>0.75));
end
thresArray = sort(thresArray);
threshold1 = thresArray(2,1);
threshold2 = thresArray(2,2);
% threshold variant 2
imPrepGray = double(grayslice(imPrepGray,255));
imValues = unique(imPrepGray);
fuzzyMat = zeros(length(imValues),3);
for valIter = 1:length(imValues)
    fuzzyMat(valIter,1) = mean(clus(1,imPrepGray(:) == imValues(valIter)));
    fuzzyMat(valIter,2) = mean(clus(2,imPrepGray(:) == imValues(valIter)));
    fuzzyMat(valIter,3) = mean(clus(3,imPrepGray(:) == imValues(valIter)));
end
repMed = max(min(repmat(fuzzyMat(:,2),1,size(fuzzyMat,1)),repmat(fuzzyMat(:,1)',size(fuzzyMat,1),1)));
repHigh = max(min(repmat(fuzzyMat(:,3),1,size(fuzzyMat,1)),repmat(fuzzyMat(:,2)',size(fuzzyMat,1),1)));
% creating rough segmentation
% medHighBW = imPrepGray>threshold2;
medHighBW = (threshold1<imPrepGray) & (imPrepGray<threshold2);
medHighBW = imerode(medHighBW,strel('disk',7));
medHighBW = imerode(medHighBW,strel('disk',7));
imshow(medHighBW)
% obtaining region properties
% props = regionprops(medHighBW,'Area','Centroid','Circularity','ConvexHull','ConvexImage','MajorAxisLength','MinorAxisLength');
% propsCell = struct2cell(props)';
% areas = cell2mat(propsCell(:,1));
% [~,idxVec] = sort(areas,'descend');
% idx = idxVec(1);
% BW = roipoly(medHighBW,props(idx).ConvexHull(:,1),props(idx).ConvexHull(:,2));
% binaryMap = activecontour(imPrepGray,BW);
imshowpair(maskCropped,medHighBW)
%% k-means z gradientu
% imForGrad = medfilt2(rgb2gray(imPrep),[5 5]);
grad = imgradient(rgb2gray(imPrep));
grad = double(grayslice(grad,7));
cat = kmeans(reshape(grad,[],1),3);
cat = reshape(cat,size(imPrep,1),size(imPrep,2));
imshow(cat,[])
% imshow(imclose(cat,[0,1,0;1,1,1;0,1,0]),[])
imEdge = (cat == 3 & grad>1);
%% geometricke kontury - nepouzivat, nefunguje podle Vicara
% stanoveni pocatecnich hranic pro vystup HT - zatim nevyuzite
% X = bwboundaries(roiMask == 1, 8); 

% stanoveni pocatecnich hranic pro vystup region growing
% sumRegion = reshape(sum(sum(segIm)),[3 1 1]);
% [~,smallObjChannel] = min(sumRegion);
% final = activecontour(rgb2gray(imCropped),segIm(:,:,smallObjChannel));%imdilate(segIm(:,:,smallObjChannel),[1 1 1; 1 1 1; 1 1 1])
% imshowpair(maskCropped,final)
%% parametricke kontury - nepouzivat, nefunguje podle Vicara
% stanoveni pocatecnich hranic pro vystup HT - zatim nevyuzite
% X = bwboundaries(roiMask == 1, 8); 

% stanoveni pocatecnich hranic pro vystup region growing
% sumRegion = reshape(sum(sum(segIm)),[3 1 1]);
% [~,smallObjChannel] = min(sumRegion);
% X = bwboundaries(imdilate(segIm(:,:,smallObjChannel),[1 1 1; 1 1 1; 1 1 1]), 8);
% % switch x and y coordinates
% X = X{1};
% X = X(:,[2,1]);
% % figure;
% % imshow(imPrep,[]); hold on
% % plot(X(:,1), X(:,2), 'r')
% 
% % Parameters
% W = 2.0;    
% alpha = 1;
% beta = 1;
% step = 1;
% num_ite = 100;
% 
% % Optimization
% G = imgaussfilt(rgb2gray(imPrep),5);
% [aGoG,pGoG] = imgradient(G);
% [Gx,Gy] = imgradientxy(aGoG.^2);
% % [Gx,Gy] = vpocitat_gvf(aGoG.^2,5,0.75); % eliminate the problem of contour in the background
% 
% % figure;
% % imshow(imPrep,[])
% % hold on
% % quiver(Gx,Gy)
% % hold off
% 
% % Gx = 20.*Gx;
% % Gy = 20.*Gy;
% 
% for ite = 1 : num_ite
%     h = sqrt(diff([X(:,1);X(1,1)]).^2 + diff([X(:,2);X(1,2)]).^2);
%     h(h==0) = 1; % due to prevent problems with 2 same points in contour
%     % Internal Force:
%     a = beta./((h.^4));
%     b = -(4*beta./(h.^4)) - (alpha./(h.^2));
%     c = (6*beta./(h.^4)) + (2*alpha./(h.^2));
%   
%     %Create A matrix
%     numP = size(X,1);
%     A = zeros(numP);
%     A = diag(c);
%     A = A + diag(b(1:numP-1),1);
%     A = A + diag(a(3:numP),-2);
%     A = A + diag(b(2:numP),-1);
%     A = A + diag(a(1:numP-2),2);
%     A(1,end)=b(1);
%     A(1,end-1)=a(1);
%     A(2,end)=a(2);
%     A(end,1)=b(end);
%     A(end-1,1)=a(end-1);
%     A(end,2)=a(end);
%     
%     % External forces
%     Fext = [];
%     for i = 1:size(X,1)
%         Fext(i,1) = Gx(floor(X(i,2)),floor(X(i,1)));
%         Fext(i,2) = Gy(floor(X(i,2)),floor(X(i,1)));
%     end
%     
%     % Balloon force
%     N = comp_normal(X);
%     w_bal = 0.1;
% %     
%     % Pohyb kontury
% %     X = (eye(numP)+step.*A)^(-1) * (X + step.*W.*Fext); % basic optimalization
% %     X = (eye(numP)+step.*A)^(-1) * (X); % the contour will collapse => we use balloon force
%     X = (eye(numP)+step.*A)^(-1) * (X + step.*W.*Fext + step.*w_bal.*N);
%     
%     imshow(imPrep,[]);
%     hold on
%     plot(X(:,1), X(:,2), 'r')
%     hold off
%     drawnow;
% 
%     % interpolation of distant contour points
% %     if ite==num_ite
% %         D =[0; cumsum(sum(abs(diff(X)),2))];
% %         X = interp1(D,X,D(1):5:D(end)); % ...to close the gaps
% %     end
% %     
%     if any(X(:,1)>(size(imPrep,2)-2)) || any(X(:,1)<2) || any(X(:,2)>(size(imPrep,1)-2)) || any(X(:,2)<2)
%         break
%     end
% 
% end
%% Detekce hran v obraze a vyplneni objektu - nefunguje
% % % tvorba masky ramecku
% % figure
% % imshow(im)
% % roi = drawpolygon;
% % borderMask = createMask(roi);
% edgedImage = edge(image,'canny',1.25*std(image(:))); % ,std(image(:))
% % [edgeRows,edgeCols] = find(edgedImage == 1);
% % labelImage = bwlabel(edgedImage);
% 
% se90 = strel('line',15,90);
% se0 = strel('line',15,0);
% se45 = strel('line',15,45);
% se135 = strel('line',15,135);
% dst = strel('disk',10,4);
% BWsdil = imdilate(edgedImage,[dst]);
% imshow(BWsdil)
% title('Dilated Gradient Mask')
% pause(3)
% % BWdfill = imfill(BWsdil,'holes');
% % BWnobord = imclearborder(BWdfill);
% % seD = strel('diamond',1);
% % BWfinal = imerode(BWnobord,seD);
% % BWfinal = imerode(BWfinal,seD);
% BWfinal = imerode(BWsdil,[se0 se90 se45 se135]);
% imshow(BWfinal)
% title('Segmented Image');


%% detekce v zóně, kde je náznak elipsovitého útvaru
% %% wavelet transform - spis na vyzkouseni
% [ca,chd,cvd,cdd] = swt2([rgb2gray(imPrep); zeros(3,size(imPrep,2))],2,'haar');
% 
% A1 = wcodemat(ca(:,:,1),255);
% H1 = wcodemat(chd(:,:,1),255);
% V1 = wcodemat(cvd(:,:,1),255);
% D1 = wcodemat(cdd(:,:,1),255);
% 
% A2 = wcodemat(ca(:,:,2),255);
% H2 = wcodemat(chd(:,:,2),255);
% V2 = wcodemat(cvd(:,:,2),255);
% D2 = wcodemat(cdd(:,:,2),255);
% 
% subplot(2,2,1)
% imagesc(A1)
% title('Approximation Coef. of Level 1')
% 
% subplot(2,2,2)
% imagesc(H1)
% title('Horizontal Detail Coef. of Level 1')
% 
% subplot(2,2,3)
% imagesc(V1)
% title('Vertical Detail Coef. of Level 1')
% 
% subplot(2,2,4)
% imagesc(D1)
% title('Diagonal Detail Coef. of Level 1')
% 
% figure
% 
% subplot(2,2,1)
% imagesc(A2)
% title('Approximation Coef. of Level 2')
% 
% subplot(2,2,2)
% imagesc(H2)
% title('Horizontal Detail Coef. of Level 2')
% 
% subplot(2,2,3)
% imagesc(V2)
% title('Vertical Detail Coef. of Level 2')
% 
% subplot(2,2,4)
% imagesc(D2)
% title('Diagonal Detail Coef. of Level 2')

%% nepouzite veci
    %     % Vykresleni polypu
    %     figure
    %     imshow(im.*mask,[])

    % Smazani odlesku
%     mira_potlaceni = 1.75;
%     im_odlesk = im>(mean(im(:))+mira_potlaceni*std(im(:))); % Práh pro odstranění odlesků
%     image = regionfill(im,im_odlesk);
    
    % Smazani odlesku - alternativa
%     image = ordfilt2(im,1,ones(5,5));

    % Smazani odlesku - alternativa 2
%       pm = rangefilt(im,true(7));
%       T = graythresh(pm);
%       reflmask = imbinarize(imfill(pm,'holes'),T);
%       image = regionfill(im,reflmask);
% 
%         % Vykresleni smazani odlesku
%         figure
%         subplot 131
%         imshow(im)
%         title('S odlesky')
%         subplot 132
%         imshow(reflmask)
%         title('Odlesky')
%         subplot 133
%         imshow(image)
%         title('Bez odlesků')
% % 
%     % Vykresleni upraveného obrazu + histogram
%     figure
%     subplot 121
%     imshow(image,[])
%     title('Original bez odlesku')
%     subplot 122
%     imhist(image)
%     title('Histogram')

    % Vypocet gradientu (1. derivace)
%     grad = imgradient(image,'prewitt');

%     figure
%     subplot 121
%     imshow(grad,[])
%     title('Gradient')
%     subplot 122
%     imhist(grad)
%     title('Histogram')

%     % Funkce Edge
%     figure
%     E = edge(image,'sobel');
%     imshow(E)

%     % Pokus o filtraci 
%     filt = [-1 -1 -1;0 0 0;-1 -1 -1];
%     edgeIm = filter2(filt,image);
%     figure
%     subplot 121
%     imshow(edgeIm,[])
%     subplot 122
%     imshow(edge(edgeIm,'zerocross'),[])
%         
%     K = 9.2;
%     h1 = (1/(K-8)).*[-1,-1,-1;-1,K,-1;-1,-1,-1];
%     sh1 = conv2(image,h1,'same');
%     figure
%     subplot 121
%     imshow(sh1,[])
%     subplot 122
%     imshow(edge(sh1,'approxcanny'),[])
    
%     % parametricka mapa pomoci std
%     pm = rangefilt(im,true(5));
%     figure
%     imshow(pm,[])
% grad = grayslice(grad,11);
% w = watershed(grad);
% 
% imshow(label2rgb(w,'jet','k'),[])

%% Lab

imLab = rgb2lab(imPrep);
figure
subplot 121
imshow(imLab(:,:,2),[])
subplot 122
imhist(zscore(imLab(:,:,2)),256)

T = graythresh(imLab(:,:,2));
figure
imshow(imLab(:,:,3)>0.6*max(imLab(:,:,3),[],'all') | imLab(:,:,2)>0.6*max(imLab(:,:,2),[],'all'),[])


%%
% Správný postup? - ořezání rámečku
%               - smazání odlesků
%               - korekce osvětlení
%               - Hlednání výchozího bodu pro segmentaci
%                    a) HT kružnic - občasné odchylky
%                    b) Lab barevný prostor - prahování v ab - jen někde
%               - RegionGrowing - prozatím nejlepší metoda s nízkým prahem
%               - Kontura - parametrická - problém s gradientem (obrysy jsou slabě kontrastní)
% Analýza na základě rozdělení do skupin podle pohledu
%Prahovaní s histerezí
% klasifikatory - randomforest