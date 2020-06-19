roiarea_thr = 85;

masks_all = masks; clear masks 
elim_masks_orig = elim_masks; elim_masks

clear masks elim_masks
% Eliminate very small rois and rois touching image border
area = zeros(size(masks_all,3),1);
borderpix = 3;
for j = 1:size(masks_all,3)
    mask = masks_all(borderpix:size(masks_all,1)-borderpix,borderpix:size(masks_all,2)-borderpix,j);
    im = imclearborder(mask);
    c = regionprops(im,'area');
    if ~isempty(c)
        area(j) = c.Area;                    % area of each ROI
    end
end
masks = masks_all(:,:,area>roiarea_thr);
elim_masks = masks_all(:,:,area<roiarea_thr);
tsG = tsG_all(area>roiarea_thr,:);
df_f = df_f_all(area>roiarea_thr,:);

% ROIs overlayed on correlation image
plotopts.plot_ids = 1; % set to 1 to view the ID number of the ROIs on the plot
fig = plotContoursOnSummaryImage(corr_image, masks, plotopts);
savefig(fig, fname_fig1(1:end-4));
saveas(fig, fname_fig1(1:end-4), 'png');
close(fig);

% eliminated ROIs overlayed on correlation image
if ~isempty(elim_masks)
    plotopts.plot_ids = 1; % set to 1 to view the ID number of the ROIs on the plot
    fig = plotContoursOnSummaryImage(corr_image, elim_masks, plotopts);
end