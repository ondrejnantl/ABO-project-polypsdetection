%% ABO - Projekt c.10 Polypy
% @JanSima,@OndrejNantl
clear all; clc;
%% nacteni
% Zmen si cestu k souboru!

pathCVC_Orig = 'D:\andyn\OneDrive - Vysoké učení technické v Brně\materialy_4r_moje\MPA-ABO\projekt\CVC-ClinicDB\Original\';
pathCVC_Mask = 'D:\andyn\OneDrive - Vysoké učení technické v Brně\materialy_4r_moje\MPA-ABO\projekt\CVC-ClinicDB\Ground Truth\';
for idx = 20
    im = rgb2gray(im2double(imread([pathCVC_Orig, num2str(idx) '.tif'])));
    imColor = im2double(imread([pathCVC_Orig, num2str(idx) '.tif']));
    mask = im2double(imread([pathCVC_Mask, num2str(idx) '.tif']));
end

%% odstraneni ramecku, smazani odlesku a korekce osvetleni

% odstraneni ramecku
clear bEdgeMask bEdgeMask2 bEdgeMask3 imCropped imCroppedRow maskCropped maskCroppedRow
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

% smazani odlesku - alternativa 2
pm = rangefilt(rgb2gray(imCropped),true(7));
T = graythresh(pm);
reflmask = imbinarize(imfill(pm,'holes'),T);
imCropped = inpaintCoherent(imCropped,logical((~bEdgeMask3).*reflmask),'SmoothingFactor',5,'Radius',5);
% imshow(imCropped)
% imCropped = inpaintExemplar(imCropped,bEdgeMask3,'PatchSize',[30 30]);
% imCropped = imgaussfilt(imCropped,0.8);
% figure
% imshow(imCropped)

% uprava osvetleni
[m,n,o] = size(imCropped);
mm = zeros(m,n,o);
N = 19;
meanMask = 1/(N^2).*ones(N,N);
for j = 1:o
    mm(:,:,j) = 0.3.*conv2(imCropped(:,:,j),meanMask,'same'); % vaha 0.3 podle Sanchez2018
end
imPrep = imCropped - mm;
figure
imshow(imPrep,[])
%% uprava preprocesovaneho obrazu - zatim nespoustet
% figure;
% for j = 1:3
%     edgedImage(:,:,j) = edge(imPrep(:,:,j),'canny',[.03 .1],sqrt(2)); % konstanty nastaveny podle Sanchez2018
%     subplot(1,3,j)
%     imshow(edgedImage(:,:,j))
% end
% edgedImage = edgedImage(:,:,1).*edgedImage(:,:,2).*edgedImage(:,:,3);

%% region growing
seedRow = 201; % je nutne vymyslet jak zjistit pozici seminka
seedCol = 181;
figure;
for i = 1:o
segIm(:,:,i) = grayconnected(imPrep(:,:,i),seedRow,seedCol,0.06);
subplot(1,3,i)
imshow(segIm(:,:,i))
end
%% Houghova transformace pro kruh - mohla by fungovat - zatim zakomentovana
% imEdge = edge(rgb2gray(imPrep),'canny',[.03 .1],sqrt(2));
% rs = 5:40;
% HS = zeros(size(imPrep,1),size(imPrep,2),length(rs));
% r_ind = 1;
% [X,Y] = find(imEdge == 1);
% for r = rs
%     tmp_c = gen_circle(r);
%     for i = 1:length(X)
%         c1 = X(i);
%         c2 = Y(i);
%         if c1 > r && c1< (size(imPrep,1) - r)
%             if c2 > r && c2< (size(imPrep,2) - r)
%                 HS((c1-r):(c1+r),(c2-r):(c2+r),r_ind) = HS((c1-r):(c1+r),(c2-r):(c2+r),r_ind)+tmp_c;
%             end
%         end
%     end
%     r_ind = r_ind + 1;
% end
% % imshow5(HS)
% [linInd] = find(HS == max(HS,[],'all'));
% [y,x,r] = ind2sub(size(HS),linInd);
% 
% imshow(imPrep);hold on; for i = 1:length(x);h = images.roi.Circle(gca,'Center',[x(i) y(i)],'Radius',r(i));end
% roiMask = createMask(h);
%% parametricke kontury
% stanoveni pocatecnich hranic pro vystup HT - zatim nevyuzite
% X = bwboundaries(roiMask == 1, 8); 

% stanoveni pocatecnich hranic pro vystup region growing
sumRegion = reshape(sum(sum(segIm)),[3 1 1]);
[~,smallObjChannel] = min(sumRegion);
X = bwboundaries(imdilate(segIm(:,:,smallObjChannel),[1 1 1; 1 1 1; 1 1 1]), 8); 

% switch x and y coordinates
X = X{1};
X = X(:,[2,1]);
figure;
imshow(imPrep,[]); hold on
plot(X(:,1), X(:,2), 'r')

% Parameters
W = 3.0;    
alpha = 1;
beta = 1;
step = 2;
num_ite = 100;


% Optimization
G = imgaussfilt(rgb2gray(imPrep),5);
[aGoG,pGoG] = imgradient(G);
[Gx,Gy] = imgradientxy(aGoG.^2);
% [Gx,Gy] = vpocitat_gvf(aGoG.^2,100,1.5); % eliminate the problem of contour in the background
figure;
imshow(imPrep,[])
hold on
quiver(Gx,Gy)
hold off

% Gx = Gx.^2;
% Gy = Gy.^2;

for ite = 1 : num_ite
    h = sqrt(diff([X(:,1);X(1,1)]).^2 + diff([X(:,2);X(1,2)]).^2);
    h(h==0) = 1; % due to prevent problems with 2 same points in contour
    % Internal Force:
    a = beta./((h.^4));
    b = -(4*beta./(h.^4)) - (alpha./(h.^2));
    c = (6*beta./(h.^4)) + (2*alpha./(h.^2));
  
    %Create A matrix
    numP = size(X,1);
    A = zeros(numP);
    A = diag(c);
    A = A + diag(b(1:numP-1),1);
    A = A + diag(a(3:numP),-2);
    A = A + diag(b(2:numP),-1);
    A = A + diag(a(1:numP-2),2);
    A(1,end)=b(1);
    A(1,end-1)=a(1);
    A(2,end)=a(2);
    A(end,1)=b(end);
    A(end-1,1)=a(end-1);
    A(end,2)=a(end);
    
    % External forces
    Fext = [];
    for i = 1:size(X,1)
        Fext(i,1) = Gx(floor(X(i,2)),floor(X(i,1)));
        Fext(i,2) = Gy(floor(X(i,2)),floor(X(i,1)));
    end
    
    % Balloon force
    N = comp_normal(X);
    w_bal = 0.07;
%     
    % Pohyb kontury
%     X = (eye(numP)+step.*A)^(-1) * (X + step.*W.*Fext); % basic optimalization
%     X = (eye(numP)+step.*A)^(-1) * (X); % the contour will collapse => we use balloon force
    X = (eye(numP)+step.*A)^(-1) * (X + step.*W.*Fext + step.*w_bal.*N);
    
    imshow(imPrep,[]);
    hold on
    plot(X(:,1), X(:,2), 'r')
    hold off
    drawnow;

    % interpolation of distant contour points
    if ite==num_ite
        D =[0; cumsum(sum(abs(diff(X)),2))];
        X = interp1(D,X,D(1):5:D(end)); % ...to close the gaps
    end
    
    if any(X(:,1)>(size(imPrep,2)-2)) || any(X(:,1)<2) || any(X(:,2)>(size(imPrep,1)-2)) || any(X(:,2)<2)
        break
    end

end
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