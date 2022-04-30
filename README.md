# ABO-project-polypsdetection
-----------------------------------------------------------------------------------------------------------------
## Authors: Ondřej Nantl, Jan Šíma, Terezie Dobrovolná
-----------------------------------------------------------------------------------------------------------------
This repository stores code of the algorithm designed for polyps detection and segmentation in colonoscopy images.
This repository was created as a semester project at FEEC at the Brno University of Technology.
The code for polyp detection and segmentation is implemented in MATLAB. 

For evaluation of detection with our algorithm run function evalPolyps whose input argument is the directory with your dataset 
(for more details see function documentation)

Method of segmentation can be changed in the evalPolyps function at line 87 (see detectPolyps help for valid inputs) 

All functions have description in function head

# List of scripts and functions:

detectPolyps - performs detection of polyps in images after cropping the black edge using other support functions

FClear - performs removing of specular highlights in the image

FHouTrans - performs search for the seed for region growing using Hough transform for circles

FHysThres - performs search for the seed for region growing using hysteresis thresholding

FLight - performs correction of variant illumination in the image

FRegionGrow - performs final segmentation using region growing technique

gen_circle - support function for FHouTrans

hyst_thresh - support function for FHysThres 

ParametricMap - script for extracting features for statistical analysis of polyp and non-polyp region

ParametricOrig - script for extracting features for statistical analysis of the whole image

polypsEval - function for evaluation of detection with our algorithm using your dataset

trainClassificator - script for creating classificator for determining confidence of our segmentation using features from detection

# List of data files:
10_presentation - presentation for project defence

2Groups - data for trainClassificator script

Parametric_Field - extracted features for statistical analysis of polyp and non-polyp region

Parametric_Field_Orig - extracted features for statistical analysis of the whole image
