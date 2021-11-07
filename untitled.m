masks_all = masks;
roiarea_min = 70;
roiarea_max = 560;
invcirc_max = 4;
overlap_thr = 0.25;
borderpix = 4;

area = zeros(size(masks_all,3),1);
inc = []; exc = []; 
for j = 1:size(masks_all,3)
    mask = masks_all(borderpix:size(masks_all,1)-borderpix,borderpix:size(masks_all,2)-borderpix,j);
    im = imclearborder(mask);
    c = regionprops(im,'area','perimeter');
    try
        if ~isempty(c)
            area(j) = c.Area;                    % area of each ROI
            invcirc(j) = (c.Perimeter.^2)/(4*pi*c.Area);
            if all([area(j)>roiarea_min,...
                    area(j)<roiarea_max,...
                    invcirc<invcirc_max])
                inc = [inc; j];
            else
                exc = [exc; j];
            end
        else
            exc = [exc; j];
        end
    catch
        exc = [exc; j];
    end
end
% eliminate highly overlapping rois
exc2 = [];
for j = 1:length(inc)
    for k = 1:length(inc)
        if j~=k
            [~, overlap1, overlap2] = calcROIoverlap(masks_all(:,:,j), masks_all(:,:,k));
            if overlap1>overlap_thr || overlap2>overlap_thr
                exc2 = [exc2; j; k];
                inc(inc == j) = [];
                inc(inc == k) = [];
            end
        end
    end
end

exc = unique([exc; exc2]);
            
            
frun_pipeline_imreg( 'list_m62_fov2_fam1nov-fam1.txt', '20181015_09_37_54', true, [0;0;0;0;0;0], [1;1;1;1;1;1], 5, 87, 99, 0.05, 0.5, 3.0 )