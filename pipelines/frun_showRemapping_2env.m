% Written by Ann Go
% 
% This script registers the ROIs from two environments and plots the place
% field tuning curves for each sorted according to their own cells and according
% to the cells of the other environment.
%
% INPUTS
% mouseid   : 'm##' e.g. 'm62'
% env1      : environment 1
% env2      : environment 2
%   e.g. 'fam1', 'fam2', 'nov', 'fam1rev'
% ref1      : image registration template file for env1
% ref2      : image registration template file for env2
%   format: 'YYYYMMDD_hh_mm_ss'
% force     : flag to force generation of comparison figures even though they
%               already exist
% figclose    : flag to close figures that will be generated (they are
%               automatically saved)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The section labeled "USER-DEFINED INPUT" requires user input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [A_union, matched_ROIs, env1PF, env2PF] = frun_showRemapping_2env_regallROIs( mouseid, env1, env2, ref1, ref2, bl1, bl2, force, fsave, figclose, histsmoothWin )

if nargin<6, bl1 = 85; end
if nargin<7, bl2 = 85; end
if nargin<8, force = false; end
if nargin<9, fsave = true; end
if nargin<10, figclose = true; end
if nargin<11, histsmoothWin = 7; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USER-DEFINED INPUT                         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Basic settings
slacknotify = false;                    % flag to send Ann slack notification re start and end of processing
mcorr_method = 'normcorre';
segment_method = 'CaImAn';
dofissa = true;
    if dofissa, str_fissa = 'FISSA'; else, str_fissa = 'noFISSA'; end

% ROI registration parameters
params.ROIreg.maxthr = [];                     
params.ROIreg.dist_maxthr = 0.1;        % threshold for turning spatial components into binary masks [default: 0.1]
params.ROIreg.dist_exp = 0.8;           % power n for distance between masked components: dist = 1 - (and(m1,m2)/or(m1,m2))^n [default: 1]
params.ROIreg.dist_thr = 0.7;           % threshold for setting a distance to infinity    [default: 0.5]
params.ROIreg.dist_overlap_thr = 0.7;   % overlap threshold for detecting if one ROI is a subset of another [default: 0.8]
params.ROIreg.plot_reg = true;
params.ROIreg.print_msg = false;

options = neuroSEE_setparams('mcorr_method', mcorr_method, 'dofissa', dofissa); 
params.ROIreg_mc.r = options.mcorr.normcorre_r;
params.ROIreg_mc.nr = options.mcorr.normcorre_nr;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load module folders and define data directory
[data_locn,~,err] = load_neuroSEEmodules(false);
if ~isempty(err)
    beep
    cprintf('Errors',err);    
    return
end


%% Pre-processing
% Check if data exist for mouseID in env1 and env2. Quit if data does not exist
dir_env1 = [data_locn 'Analysis/' mouseid '/' mouseid '_' env1 env2 '-' env1 ...
            '/group_proc/imreg_' mcorr_method '_' segment_method '/' mouseid '_' env1 env2 '-' env1 '_imreg_ref' ref1 '/'];
dir_env2 = [data_locn 'Analysis/' mouseid '/' mouseid '_' env1 env2 '-' env2 ...
            '/group_proc/imreg_' mcorr_method '_' segment_method '/' mouseid '_' env1 env2 '-' env2 '_imreg_ref' ref2 '/'];

data_env1 = [dir_env1  str_fissa '/bl_prctile' num2str(bl1) '/' mouseid '_' env1 env2 '-' env1 '_ref' ref1 '_PFmap_output.mat'];
data_env2 = [dir_env2  str_fissa '/bl_prctile' num2str(bl2) '/' mouseid '_' env1 env2 '-' env2 '_ref' ref2 '_PFmap_output.mat'];
       
if ~and(exist(data_env1,'file'), exist(data_env2,'file'))       
    beep
    err = sprintf('No processed data for %s in %s & %s\n', mouseid, env1, env2);
    cprintf('Errors',err);    
    return
end

% Send Ann slack message if processing has started
if slacknotify
    slacktext = [mouseid '_' env1 env2 ': processing started'];
    neuroSEE_slackNotify( slacktext );
end


%% ROI registration across sessions
tic
fdir = [data_locn 'Analysis/' mouseid '/' mouseid '_' env1 env2 '/remapping/imreg_' ...
       mcorr_method '_' segment_method '_' str_fissa '/' mouseid '_' env1 env2 '_imreg_ref' ref1 '-' ref2 '/'];
fname_remap = [fdir  'bl_' num2str(bl1) '-' num2str(bl2) '/' mouseid '_' env1 env2 '_remapping_output.mat'];
fname_remapfig = [fdir  'bl_' num2str(bl1) '-' num2str(bl2) '/' mouseid '_' env1 env2 '_remapping_summary.fig'];

if ~exist(fname_remap,'file') || force 
    fprintf('%s: loading data\n',[mouseid '_' env1 env2]);
    
    M1 = load([dir_env1 mouseid '_' env1 env2 '-' env1 '_ref' ref1 '_segment_output.mat']);
    M2 = load([dir_env2 mouseid '_' env1 env2 '-' env2 '_ref' ref2 '_segment_output.mat']);
    
    PF1 = load(data_env1);
    PF2 = load(data_env2);

    masks1 = M1.masks;
    masks2 = M2.masks;
    A1 = zeros(size(masks1,1)*size(masks1,2),size(masks1,3));
    A2 = zeros(size(masks2,1)*size(masks2,2),size(masks2,3));
    for i = 1:size(masks1,3)
        mask1 = masks1(:,:,i);
        A1(:,i) = mask1(:);
    end
    for i = 1:size(masks2,3)
        mask2 = masks2(:,:,i);
        A2(:,i) = mask2(:);
    end
    A1 = sparse(A1);
    A2 = sparse(A2);

    params.ROIreg.d1 = size(masks1,1);
    params.ROIreg.d2 = size(masks1,2);
    params.ROIreg.plot_reg = true;
    params.ROIreg.print_msg = true;

    t1 = load([dir_env1 mouseid '_' env1 env2 '-' env1 '_ref' ref1 '_imreg_template.mat']);
    t2 = load([dir_env2 mouseid '_' env1 env2 '-' env2 '_ref' ref2 '_imreg_template.mat']);
    template1 = t1.template_g;
    template2 = t2.template_g;

    fprintf('%s: registering ROIs\n',[mouseid '_' env1 env2]);
    if fsave
        fname_fig = [fdir  'bl_' num2str(bl1) '-' num2str(bl2) '/' mouseid '_' env1 env2 '_regROIs_output'];
    else
        fname_fig = [];
    end
    [ ~, ~, ~, A2_shifted, ~, ~, template2_shifted ] = ...
       register_ROIs( A1, A2, params.ROIreg, template1, template2, params.ROIreg_mc.r, [], false );
    [ matched_ROIs, nonmatched_1, nonmatched_2, A2, R, A_union, template2_shifted ] = ...
       register_ROIs( A1, A2_shifted, params.ROIreg, template1, template2_shifted, params.ROIreg_mc.nr ,fname_fig, false );
%     [ matched_ROIs, nonmatched_1, nonmatched_2, A2, R, A_union ] = ...
%         register_ROIs( A1, A2, params.ROIreg, template1, template2, params.ROIreg_mc, fname_fig, true );
    masks_union = reshape(full(A_union), params.ROIreg.d1, params.ROIreg.d2, size(A_union,2));
    masks2_reg = reshape(full(A2), params.ROIreg.d1, params.ROIreg.d2, size(A2,2));
    for j = 1:size(masks_union,3)
        outlines{:,:,j} = bwboundaries(masks_union(:,:,j));    % boundary of each ROI
    end

    fprintf('%s: processing remapping\n',[mouseid '_' env1 env2]);
    % create image of all PC masks overlaid by fam1 maks
    im_nonmatched2 = zeros(512,512,3);
    c = 0.5;
    for ii = 1:numel(nonmatched_2)
        k = nonmatched_2(ii);
        imR = masks2(:,:,k)*c;
        imG = masks2(:,:,k)*c;
        imB = masks2(:,:,k)*c;
        im_nonmatched2 = im_nonmatched2 + cat(3,imR,imG,imB); 
        im_nonmatched2(im_nonmatched2 > c) = c;
    end

    im_env1masks = im_nonmatched2;
    for ii = 1:size(masks1,3)
        imR = masks1(:,:,ii)*1; 
        imG = masks1(:,:,ii)*0;
        imB = masks1(:,:,ii)*0;
        im_env1masks = im_env1masks + cat(3,imR,imG,imB); 
    end

    % create image of all PC masks overlaid by fam1 masks
    im_nonmatched1 = zeros(512,512,3);
    c = 0.5;
    for ii = 1:numel(nonmatched_1)
        k = nonmatched_1(ii);
        imR = masks1(:,:,k)*c;
        imG = masks1(:,:,k)*c;
        imB = masks1(:,:,k)*c;
        im_nonmatched1 = im_nonmatched1 + cat(3,imR,imG,imB); 
        im_nonmatched1(im_nonmatched1 > c) = c;
    end
    
    im_env2masks = im_nonmatched1;
    for ii = 1:size(masks2,3)
        imR = masks2_reg(:,:,ii)*0;
        imG = masks2_reg(:,:,ii)*0;
        imB = masks2_reg(:,:,ii)*1;
        im_env2masks = im_env2masks + cat(3,imR,imG,imB); 
    end

    % create image of all PC masks overlaid by masks common to fam1 & fam2
    im_nonmatched = im_nonmatched1 + im_nonmatched2;
    im_nonmatched(im_nonmatched > c) = c;
    
    im_env1env2masks = im_nonmatched;
    for ii = 1:size(matched_ROIs,1)
        k = matched_ROIs(ii,1);
        imR = masks1(:,:,k)*1;
        imG = masks1(:,:,k)*0;
        imB = masks1(:,:,k)*1;
        im_env1env2masks = im_env1env2masks + cat(3,imR,imG,imB); 
    end
    
    % PF maps for fam1 and fam2 in their own sorting
    env1PF = zeros(size(PF1.hist.SIspk.sort_pfMap));
    for j = 1:length(PF1.hist.SIspk.pcIdx)
        if PF1.params.histsmoothWin == histsmoothWin
            env1PF(j,:) = PF1.hist.normrateMap_sm(PF1.hist.SIspk.pcIdx(j),:);
        else
            ppp = PF1.hist.rateMap(PF1.hist.SIspk.pcIdx(j),:);
            ppp_sm = circularSmooth(ppp,histsmoothWin);
            env1PF(j,:) = ppp_sm./max(ppp_sm);
        end
    end
    
    env2PF = zeros(size(PF2.hist.SIspk.sort_pfMap));
    for j = 1:length(PF2.hist.SIspk.pcIdx)
        if PF2.params.histsmoothWin == histsmoothWin
            env2PF(j,:) = PF2.hist.normrateMap_sm(PF2.hist.SIspk.pcIdx(j),:);
        else
            ppp = PF2.hist.rateMap(PF2.hist.SIspk.pcIdx(j),:);
            ppp_sm = circularSmooth(ppp,histsmoothWin);
            env2PF(j,:) = ppp_sm./max(ppp_sm);
        end
    end
    
    % PF maps for fam2 in fam1 sorting
    env1PF_env2Sorting = zeros(size(env2PF));
    for j = 1:length(PF2.hist.SIspk.pcIdx)
       [matched,ind] = ismember(PF2.hist.SIspk.pcIdx(j),matched_ROIs(:,2));
       if matched
           if PF2.params.histsmoothWin == histsmoothWin
               env1PF_env2Sorting(j,:) = PF1.hist.normrateMap_sm(matched_ROIs(ind,1),:);
           else
               ppp = PF1.hist.rateMap(matched_ROIs(ind,1),:);
               ppp_sm = circularSmooth(ppp,histsmoothWin);
               env1PF_env2Sorting(j,:) = ppp_sm./max(ppp_sm);
           end
       end
    end
    
    % PF maps for fam2 in fam1 sorting
    env2PF_env1Sorting = zeros(size(env1PF));
    for j = 1:length(PF1.hist.SIspk.pcIdx)
       [matched,ind] = ismember(PF1.hist.SIspk.pcIdx(j),matched_ROIs(:,1));
       if matched
           if PF1.params.histsmoothWin == histsmoothWin
               env2PF_env1Sorting(j,:) = PF2.hist.normrateMap_sm(matched_ROIs(ind,2),:);
           else
               ppp = PF2.hist.rateMap(matched_ROIs(ind,2),:);
               ppp_sm = circularSmooth(ppp,histsmoothWin);
               env2PF_env1Sorting(j,:) = ppp_sm./max(ppp_sm);
           end
       end
    end
    clear matched ind
    
    % sort maps
    [ env1_pfLoc, ~, ~ ] = prefLoc_fieldSize_1d( env1PF );
    [ ~, env1_sortIdx ] = sort( env1_pfLoc );
    env1PF = env1PF(env1_sortIdx,:);
    env2PF_env1Sorting = env2PF_env1Sorting(env1_sortIdx,:);
    
    [ env2_pfLoc, ~, ~ ] = prefLoc_fieldSize_1d( env2PF );
    [ ~, env2_sortIdx ] = sort( env2_pfLoc );
    env2PF = env2PF(env2_sortIdx,:);
    env1PF_env2Sorting = env1PF_env2Sorting(env2_sortIdx,:);
    
    % save data
    if fsave
        remapping_output.masks = masks_union;
        remapping_output.outlines = outlines;
        remapping_output.matched_ROIs = matched_ROIs;
        remapping_output.nonmatched_env1 = nonmatched_1;
        remapping_output.nonmatched_env2 = nonmatched_2;
        remapping_output.alignment_matrix = R;
        remapping_output.aligned_env2_ROIs = A2;
        remapping_output.masks2_reg = masks2_reg;
        remapping_output.params = params.ROIreg;
        remapping_output.template_env1 = template1;
        remapping_output.template_env2 = template2_shifted;

        remapping_output.im_env1masks = im_env1masks;
        remapping_output.im_env2masks = im_env2masks;
        remapping_output.im_env1env2masks = im_env1env2masks;
        remapping_output.env1PF = env1PF;
        remapping_output.env2PF_env1Sorting = env2PF_env1Sorting;
        remapping_output.env1PF_env2Sorting = env1PF_env2Sorting;
        remapping_output.env2PF = env2PF;

        fprintf('%s: saving remapping data\n',[mouseid '_' env1 env2]);
        if ~exist(fdir,'dir'), mkdir(fdir); fileattrib fdir +w '' s; end
        save(fname_remap, '-struct', 'remapping_output')
    end
else
    if ~exist(fname_remapfig,'file') 
        fprintf('%s: loading remapping data\n',[mouseid '_' env1 env2]);
        c = load(fname_remap);
        % nonmatched_1 = c.nonmatched_env1;
        % nonmatched_2 = c.nonmatched_env2;
        im_env1masks = c.im_env1masks;
        im_env2masks = c.im_env2masks;
        im_env1env2masks = c.im_env1env2masks;
        env1PF = c.env1PF;
        env2PF_env1Sorting = c.env2PF_env1Sorting;
        env1PF_env2Sorting = c.env1PF_env2Sorting;
        env2PF = c.env2PF;
        clear c
    end
end

if ~exist(fname_remapfig,'file') || force
    fh = figure;
    fontsize = 16;
    Nbins = size(env1PF,2);
    subplot(2,12,1:4);
        imshow(im_env1masks); title(env1,'Fontweight','normal','Fontsize',fontsize);    
    subplot(2,12,5:8);
        imshow(im_env2masks); title(env2,'Fontweight','normal','Fontsize',fontsize);
    subplot(2,12,9:12);
        imshow(im_env1env2masks); title([env1 ' \cap ' env2],'Fontweight','normal','Fontsize',fontsize);
    subplot(2,12,13:15);
        cmap = viridisMap;
        imagesc(env1PF); 
        colormap(cmap); %colorbar
        title(env1,'Fontweight','normal','Fontsize',fontsize);
        yticks([1 size(env1PF,1)]); yticklabels([1 size(env1PF,1)]);
        xticks([1 Nbins]); xticklabels([1 100]);
        xlabel('Position (cm)'); %ylabel('Cell no.');
    subplot(2,12,16:18);
        imagesc(env2PF_env1Sorting); 
        title(env2,'Fontweight','normal','Fontsize',fontsize);
        yticks([]);
        xticks([1 Nbins]); xticklabels([1 100]);
        xlabel('Position (cm)');
    subplot(2,12,19:21);
        imagesc(env1PF_env2Sorting); 
        title(env1,'Fontweight','normal','Fontsize',fontsize); 
        yticks([1 size(env2PF,1)]); yticklabels([1 size(env2PF,1)]);
        xticks([1 Nbins]); xticklabels([1 100]);
        xlabel('Position (cm)');
    subplot(2,12,22:24);
        imagesc(env2PF); 
        title(env2,'Fontweight','normal','Fontsize',fontsize); 
        yticks([]);
        xticks([1 Nbins]); xticklabels([1 100]);
        xlabel('Position (cm)');

    if fsave    
        fprintf('%s: saving remapping summary figure\n',[mouseid '_' env1 env2]);
        savefig( fh, fname_remapfig(1:end-4) );
        saveas( fh, fname_remapfig(1:end-4), 'png' );
        if figclose, close( fh ); end   
    end
end

    
t = toc;
str = sprintf('%s: Processing done in %g hrs\n', [mouseid '_' env1 env2], round(t/3600,2));
cprintf(str)
end



