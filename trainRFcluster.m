%% ABO - Projekt c.10 Polypy - tvorba random forest klasifikatoru, shlukovani
% @JanSima,@OndrejNantl,@TerezieDobrovolna
clear all; clc;
%% urceni cesty
pathCVC_Orig = 'D:\andyn\OneDrive - Vysoké učení technické v Brně\materialy_4r_moje\MPA-ABO\projekt\CVC-ClinicDB\Original\';
pathCVC_Sort = 'D:\andyn\OneDrive - Vysoké učení technické v Brně\materialy_4r_moje\MPA-ABO\projekt\CVC-ClinicDB-sorted\Original\';
pathCVC_Mask = 'D:\andyn\OneDrive - Vysoké učení technické v Brně\materialy_4r_moje\MPA-ABO\projekt\CVC-ClinicDB\Ground Truth\';
%% zisk kategorii obrazu
contOrig = dir([pathCVC_Orig '*.tif']);
contOrig = {contOrig.name};
cats = zeros(length(contOrig),1);
catNames ={'nezaraditelne','primo','prurez','zboku'};
for i = 1:length(catNames)
        contSorted{i,1} = dir([pathCVC_Sort catNames{i} '\*.tif']);
        contSorted{i,1} = {contSorted{i}.name}';
end
for j = 1:length(cats)
    for k = 1:length(catNames)
        if any(ismember(contSorted{k},contOrig{j}))
            cats(j) = k;
        end
    end
end

%% tvorba a ulozeni stromu
load('Parametric_Field_Original.mat')
Mdl = TreeBagger(100,ParametricFieldOrig,cats,'Method','classification');

%% shlukovani skupin podle parametru
[ParFieldN,mu,sigma] = zscore(ParametricFieldOrig);
% tvorba boxplotu - nestandardizovane
figure
for i = 1:size(ParametricFieldOrig,2)
    subplot(4,3,i)
    boxplot(ParametricFieldOrig(:,i),cats,'Labels',catNames);
    title(Labels{i})
end
% tvorba boxplotu - standardizovane
figure
for i = 1:size(ParametricFieldOrig,2)
    subplot(4,3,i)
    boxplot(zscore(ParametricFieldOrig(:,i)),cats,'Labels',catNames);
    title(Labels{i})
end

%shlukovani
k = 4;
[imClusters,centroids] = kmeans(ParametricFieldOrig,k);
centroidsN = (centroids - repmat(mu,k,1))./repmat(sigma,k,1);

% centroidy
figure
subplot 211
plot(centroids'); 
legend({'Cluster 1','Cluster 2','Cluster 3','Cluster 4'});
title('Hodnoty jednotlivých příznaků pro centroidy - nenormalizovano');
xticks(1:12)
xticklabels(Labels)
subplot 212
plot(centroidsN'); 
legend({'Cluster 1','Cluster 2','Cluster 3','Cluster 4'});
title('Hodnoty jednotlivých příznaků pro centroidy - normalizovano');
xticks(1:12)
xticklabels(Labels)


% extrahovani obsahu clusteru
contClust = cell(length(unique(imClusters)),1);
for m = 1:length(unique(imClusters))
    contClust{m,1} = contOrig(imClusters == m)';
end
% ulozeni vsech vystupu
save('RFandClust.mat')