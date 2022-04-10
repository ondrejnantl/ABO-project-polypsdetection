%% ABO - Projekt c.10 Polypy - tvorba random forest klasifikatoru
% @JanSima,@OndrejNantl
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
        contSorted{i} = dir([pathCVC_Sort catNames{i} '\*.tif']);
        contSorted{i} = {contSorted{i}.name}';
end
for j = 1:length(cats)
    for i = 1:length(catNames)
        if any(ismember(contSorted{i},contOrig{j}))
            cats(j) = i;
        end
    end
end
%% nacteni parametrickeho pole
load('Parametric_Field.mat')
%% tvorba stromu
Tree = TreeBagger(100,ParametricField,cats,'OOBPrediction','On','Method','classification');
save('RF.mat','Tree','cat')