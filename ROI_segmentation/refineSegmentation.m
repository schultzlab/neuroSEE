% Written by Ann Go
%
% This function refines segementation by rejecting some ROIs based on area,
% saturation, and noise 
%
% INPUTS
%   cell_tsG    : time series from green channel
%   cell_tsR    : time series from red channel
%   cell_masks  : ROI masks for cells
%   settings.
%       satThresh   : ROI is accepted if its fluorescence is below satThresh
%       satTime     :   satTime fraction of the time
%       areaThresh  : min acceptable ROI area (pixels)
%       noiseThresh : max acceptable ROI noise level
%
% OUTPUTS
%   tsG         : green channel time series of accepted ROIs
%   tsR         : red channel time series of accepted ROIs
%   masks       : masks of accepted ROIs
%   fig         : fig handle of mean image with ROIs


function [tsG, masks] = refineSegmentation(cell_tsG, cell_masks, settings)
    
    str = sprintf('%s: Refining ROIs\n', file);
    cprintf(str);
   
    % Area thresholding
    Numcells = size(cell_masks,3);
    maskArea = zeros(1,Numcells);
    for i = 1:Numcells
        maskArea(i) = bwarea(cell_masks(:,:,i));
    end
    i_area = find( maskArea >= settings.areaThresh );
    
    % Noise thresholding
    maskNoise = zeros(1,Numcells);
    for i = 1:Numcells
        y = cell_tsG(i,:);
        maskNoise(i) = GetSn(y,[0.25,0.5],'logmexp');
    end
    i_snr = find( maskNoise <= settings.noiseThresh );
    
    % Saturation thresholding
    ROIsat = (cell_tsG >= settings.satThresh);
    i_sat = find( mean(ROIsat,ndims(cell_tsG)) < settings.satTime );
        
    ind = intersect(intersect(i_area,i_snr),i_sat);
    tsG = cell_tsG(ind,:);
    masks = cell_masks(:,:,ind);    
    
    currstr = sprintf( '%s: Final ROIs identified. %g ROIs rejected\n', file, size(cell_masks,3)-size(masks,3) );
    refreshdisp( currstr, str );

end