function [ Nrois, Npcs ] = frun_getExpNroisNpcs( list, reffile )

mcorr_method = 'normcorre';  % motion correction method for individual image files
imreg_method = 'normcorre';  % image registration method [normcorre, normcorre-r, normcorre-nr, fftRigid] 
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

% Mouseid, Experiment name, files
[ mouseid, expname ] = find_mouseIDexpname(list);

% Location of processed group data for list
if strcmpi(imreg_method, mcorr_method)
    grp_sdir = [data_locn 'Analysis/' mouseid '/' mouseid '_' expname '/group_proc/'...
                groupreg_method '_' imreg_method '_' segment_method '_' str_fissa '/'...
                mouseid '_' expname '_imreg_ref' reffile '/'];
else
    grp_sdir = [data_locn 'Analysis/' mouseid '/' mouseid '_' expname '/group_proc/'...
                groupreg_method '_' imreg_method '_' segment_method '_' str_fissa '/'...
                mouseid '_' expname '_imreg_ref' reffile '_' mcorr_method '/'];
end
    
if exist([grp_sdir mouseid '_' expname '_ref' reffile '_segment_output.mat'],'file')
    M = load([grp_sdir mouseid '_' expname '_ref' reffile '_segment_output.mat']);
    Nrois = size(M.masks,3);
else
    Nrois = [];
end
if exist([grp_sdir mouseid '_' expname '_ref' reffile '_PFmap_output.mat'],'file')
    N = load([grp_sdir mouseid '_' expname '_ref' reffile '_PFmap_output.mat']);
    Npcs = numel(N.hist.SIspk.pcIdx);
else
    Npcs = [];
end