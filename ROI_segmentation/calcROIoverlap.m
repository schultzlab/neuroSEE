% Written by Ann Go
% 
% This script determines whether 2 ROIs overlap and if so, quantifies the
% overlap.
%
% INPUTS
%   mask1
%   mask2
% OUTPUTS
%   TF : overlap?
%   overlap1 : overlap region as percentage of area of mask1
%   overlap2 : overlap region as percentage of area of mask2

function [TF, overlap1, overlap2] = calcROIoverlap(mask1, mask2)

mask1 = double(mask1);
mask2 = double(mask2);

sum = mask1+mask2;
if max(max(sum))>1
    TF = true;
    c1 = regionprops(mask1,'PixelIdxList');
    c2 = regionprops(mask2,'PixelIdxList');
    pix1 = c1.PixelIdxList;
    pix2 = c2.PixelIdxList;
    com = intersect(pix1,pix2);
    overlap1 = length(com)/length(pix1);
    overlap2 = length(com)/length(pix2);
else
    TF = false;
    overlap1 = 0;
    overlap2 = 0;
end