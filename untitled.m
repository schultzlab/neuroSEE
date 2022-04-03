masks_all = masks;
% masks_all(:,:,size(masks_all,3)+1:size(masks,3)+size(elim_masks,3)) = elim_masks;
elim_masks0 = elim_masks;
tsG_all(1:size(masks,3),:) = tsG;
%tsG_all(size(masks,3)+1:size(masks,3)+size(elim_masks,3),:) = elim_tsG;
df_f_all(1:size(masks,3),:) = df_f;
%df_f_all(size(masks,3)+1:size(masks,3)+size(elim_masks,3),:) = elim_df_f;

roiarea_min = 70; %50;
roiarea_max = 560; %400;
invcirc_max = 4;
overlap_thr = 0.25;
borderpix = 4;

area = zeros(size(masks_all,3),1);
invcirc = zeros(size(masks_all,3),1);
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
                    invcirc(j)<invcirc_max])
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
            
masks = masks_all(:,:,inc);     elim_masks = masks_all(:,:,exc);    elim_masks(:,:,size(elim_masks,3)+1:size(elim_masks,3)+size(elim_masks0,3)) = elim_masks0;
tsG = tsG_all(inc,:);           elim_tsG = tsG_all(exc,:);
df_f = df_f_all(inc,:);         elim_df_f = df_f_all(exc,:);

% Save output
output.incmasks = inc;  output.excmasks = exc;
output.tsG = tsG;       output.elim_tsG = elim_tsG;
output.df_f = df_f;     output.elim_df_f = elim_df_f;
output.masks = masks;   output.elim_masks = elim_masks;
output.corr_image = corr_image;
output.F0 = F0;
output.A = A;
output.params = params;       
save('m62_fov2_fam1_s1-4_ref20181015_09_37_54_segment_output.mat','-struct','output');

%%%%%%%%
a = 5*7420+1; b = 10*7420;
df_f = df_f(:,a:b);
elim_df_f = elim_df_f(:,a:b);
tsG = tsG(:,a:b);
elim_tsG = elim_tsG(:,a:b);
clear a b
save('m77_fam1fam2-fam2_ref20190302_15_38_59_segment_output.mat')

a = 12*7420+1; b = 17*7420;
ddf_f = ddf_f(:,a:b);
dtsG = dtsG(:,a:b);
clear a b
save('m125_fam1fam2fam1-fam1r2_ref20210907_18_11_49_fissa_output.mat')

frun_pipeline_imreg( 'list_m111_fam1fam2-fam1.txt', '20201208_14_55_28', true, [0;0;0;0;0;0], [1;1;1;1;1;1], 5, 89, 99, 0.03, 0.3, 2.5 )

frun_pipeline_imreg( 'list_m125_fam1fam2fam1-fam2.txt', '20210907_18_11_49', true, [0;0;0;0;0;0], [1;1;1;1;1;1], 5, 87, 99, 0.03, 0.3, 2.5, true )

frun_pipeline_imreg( 'list_m77_fam1fam2-fam1.txt', '20190302_15_38_59', true, [0;0;0;1;0;0], [1;1;1;1;1;1], 5, 85, 99, 0.03, 0.3 )
frun_pipeline_imreg( 'list_m70_fam1nov-fam1.txt', '20181101_13_09_55', true, [0;0;0;0;0;0], [1;1;1;1;1;1], 5, 90, 99, 0.05, 0.5 )

[ activity, meanspeed ] = getActivityMeanspeed( 'list_m62_BT.txt', 13, true )
frun_collate_indivproc_results( 'list_m129_fov1_fam1novfam1-fam1.txt' )
frun_mcorr_batch( 3, 'list_m129_fov1_fam1novfam1.txt', 'normcorre', true, '20211009_14_54_57', 'green', 30, 30 )
frun_mcorr_batch( 6, 'list_m134_fov2_fam1fam1revfam1.txt', 'normcorre', true, [], 'green', 30, 30, 3, false, 2  )