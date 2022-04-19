function imPrep = FLight(inputImage)
% uprava osvetleni
[m,n,o] = size(inputImage);
mm = zeros(m,n,o);
N = 20;
meanMask = 1/(N^2).*ones(N,N);
for j = 1:o
    mm(:,:,j) = 0.3.*conv2(inputImage(:,:,j),meanMask,'same'); % vaha 0.3 podle Sanchez2018 - lehce zmenena
end
imPrep = inputImage - mm;
% imPrepGray = rgb2gray(imPrep);
% imPrepHSV = rgb2hsv(imPrep);
% imPrepLab = rgb2lab(imPrep);
end