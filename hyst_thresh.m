% HYSTHRESH - Hysteresis thresholding
% Arguments:
%             im  - image to be thresholded (assumed to be non-negative)
%             T1  - upper threshold value
%             T2  - lower threshold value
%                   (T1 and T2 can be entered in any order, the larger of the
%                   two values is used as the upper threshold)
% Function performs hysteresis thresholding of an image.
% All pixels with values above threshold T1 are marked as edges
% All pixels that are connected to points that have been marked as edges
% and with values above threshold T2 are also marked as edges. Eight
% connectivity is used.
% Reference:Peter Kovesi
% July 2005

function bw = hyst_thresh(imPre, T1, T2)
    if T1 >  T2    % T1 and T2 reversed - swap values 
	tmp = T1;
	T1 = T2; 
	T2 = tmp;
    end
    
  aboveT2 = imPre > T2;                  % Edge points above lower threshold. 
  [aboveT1r, aboveT1c] = find(imPre > T1);  % Row and colum coords of points
                                           % above upper threshold.			   
    % Obtain all connected regions in aboveT2 that include a point that has a
    % value above T1 
  bw = bwselect(aboveT2, aboveT1c, aboveT1r, 8);
