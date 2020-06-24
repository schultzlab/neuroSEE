function [ Nrois, Npcs ] = frun_getExpNroisNpcs( list, reffile )

mcorr_method = 'normcorre';  % image registration method [normcorre, normcorre-r, normcorre-nr, fftRigid] 
groupreg_method = 'imreg';      % method for concatenating file data (either register images or rois)
segment_method = 'CaImAn';      % [ABLE,CaImAn]    
dofissa = true;                 % flag to implement FISSA (when false, overrides force(3) setting)
    if dofissa, str_fissa = 'FISSA'; else, str_fissa = 'noFISSA'; end
    
% Load module folders and define data directory
[data_locn, ~, err] = load_neuroSEEmodules;
if ~isempty(err)
    beep
    cprintf('Errors',err);    
    return
end

params = neuroSEE_setparams;
roiarea_thr = params.ROIsegment.roiarea_thr;

% Mouseid, Experiment name, files
[ mouseid, expname ] = find_mouseIDexpname(list);

% Location of processed group data for list
grp_sdir = [data_locn 'Analysis/' mouseid '/' mouseid '_' expname '/group_proc/'...
            groupreg_method '_' mcorr_method '_' segment_method '/'...
            mouseid '_' expname '_imreg_ref' reffile '/'];
    
if exist([grp_sdir mouseid '_' expname '_ref' reffile '_segment_output.mat'],'file')
    M = load([grp_sdir mouseid '_' expname '_ref' reffile '_segment_output.mat']);
    masks_all = M.masks;
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
    Nrois = length(area>roiarea_thr);
else
    Nrois = [];
end
if exist([grp_sdir '/' str_fissa '/' mouseid '_' expname '_ref' reffile '_PFmap_output.mat'],'file')
    N = load([grp_sdir '/' str_fissa '/' mouseid '_' expname '_ref' reffile '_PFmap_output.mat']);
    Npcs = numel(N.hist.SIspk.pcIdx);
else
    Npcs = [];
end