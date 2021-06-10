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

mask1 = masks(:,:,1);
mask2 = masks(:,:,2);
TF = overlaps(mask1, mask2)
