%% USER INPUT
list = 'list_m66_fam1fam2-fam2.txt';
reffile = '20181013_14_50_18';
imreg_method = 'normcorre';
mcorr_method = 'normcorre';
segment_method = 'CaImAn';
dofissa = true;
groupreg_method = 'imreg';

roiarea_thr = 70;           % roi area smaller than this will be eliminated
                            % neuroSEE_segment already filtered areas <70
borderpix = 4;              % thickness (in pix) of image border to be cleared of rois
removerois = [];            % specific rois to eliminate (e.g. overlapping rois)

%% Load module folders and define data directory
[data_locn,~,err] = load_neuroSEEmodules;
if ~isempty(err)
    beep
    cprintf('Errors',err);    
    return
end

%% MouseID and experiment name
[ mouseid, expname ] = find_mouseIDexpname(list);

%% Location of processed group data for list
if dofissa, str_fissa = 'FISSA'; else, str_fissa = 'noFISSA'; end

if strcmpi(imreg_method, mcorr_method)
    grp_sdir = [data_locn 'Analysis/' mouseid '/' mouseid '_' expname '/group_proc/'...
                groupreg_method '_' imreg_method '_' segment_method '_' str_fissa '/'...
                mouseid '_' expname '_imreg_ref' reffile '/'];
else
    grp_sdir = [data_locn 'Analysis/' mouseid '/' mouseid '_' expname '/group_proc/'...
                groupreg_method '_' imreg_method '_' segment_method '_' str_fissa '/'...
                mouseid '_' expname '_imreg_ref' reffile '_' mcorr_method '/'];
end

%% Manually eliminate rois from ROI segmentation data
load([grp_sdir mouseid '_' expname '_ref' reffile '_segment_output.mat'])    

% Eliminate rois with area < roiarea_thr and rois touching image border
% within borderpix pixels
area = zeros(size(masks_orig,3),1);
for j = 1:size(masks_orig,3)
    mask = masks(borderpix:size(masks_orig,1)-borderpix,borderpix:size(masks,2)-borderpix,j);
    im = imclearborder(mask);
    c = regionprops(im,'area');
    if ~isempty(c)
        area(j) = c.Area;                    % area of each ROI
    end
end
exc_idx = [removerois; find(area<roiarea_thr)];
elim_masks = cat(3,elim_masks,masks(:,:,exc_idx));

inc_idx = [];
for k = 1:size(masks,3)
    if ~ismember(k,exc_idx)
        inc_idx = [inc_idx; k];
    end
end
masks = masks(:,:,inc_idx);
tsG = tsG(inc_idx,:);
df_f = df_f(inc_idx,:);

% Save data
output.corr_image = corr_image;
output.elim_masks = elim_masks;
output.masks = masks;
output.tsG = tsG;
output.df_f = df_f;
output.F0 = F0;
output.params = params;
if exist(GUIdata,'var'), output.GUIdata = GUIdata; end
save([grp_sdir mouseid '_' expname '_ref' reffile '_segment_output.mat'],'-struct','output');

% Save revised figures
plotROIsegmentdata(corr_image, masks, elim_masks, tsG, df_f, [grp_sdir mouseid '_' expname '_ref' reffile])
clear corr_image elim_masks masks tsG df_f F0 GUIdata params

%% Revise fissa data
if exist([grp_sdir mouseid '_' expname '_ref' reffile '_fissa_output.mat'],'file')
    load([grp_sdir mouseid '_' expname '_ref' reffile '_fissa_output.mat']);
    
    dtsG = dtsG(inc_idx,:);
    ddf_f = ddf_f(inc_idx,:);
    output.dtsG = dtsG;
    output.ddf_f = ddf_f;
    output.params = params;
    save([grp_sdir mouseid '_' expname '_ref' reffile '_fissa_output.mat'],'-struct','output');
    
    multiplot_ts(dtsG, [grp_sdir mouseid '_' expname '_ref' reffile '_fissa_result.fig'], 'Fissa-corrected raw timeseries');
    multiplot_ts(ddf_f, [grp_sdir mouseid '_' expname '_ref' reffile '_fissa_df_f.fig'], 'Fissa-corrected dF/F');
    clear dtsG ddf_f params
end

%% Revise spike data
if exist([grp_sdir mouseid '_' expname '_ref' reffile '_spikes.mat'],'file')
    load([grp_sdir mouseid '_' expname '_ref' reffile '_spikes.mat']);
    
    spikes = spikes(inc_idx,:);
    output.spikes = spikes;
    output.params = params;
    save([grp_sdir mouseid '_' expname '_ref' reffile '_spikes.mat'],'-struct','output');
    
    plotSpikes(spikes, [grp_sdir mouseid '_' expname '_ref' reffile '_spikes'])
    clear spikes params
end
