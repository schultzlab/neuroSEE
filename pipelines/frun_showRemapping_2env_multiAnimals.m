% Written by Ann Go
% 
% This script registers the ROIs from two environments and plots the place
% field tuning for each sorted according to their own cells and according
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The section labeled "USER-DEFINED INPUT" requires user input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function frun_showRemapping_2env_multiAnimals( mouseid_array, env1, env2, ref1_array, ref2_array, force, figclose )

if nargin<6, force = false; end
if nargin<7, figclose = true; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USER-DEFINED INPUT                         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Basic settings
default = false;                 % flag to use default parameters
mcorr_method = 'normcorre-nr';    % values: [normcorre, normcorre-r, normcorre-nr, fftRigid] 
                                    % CaImAn NoRMCorre method: 
                                    %   normcorre-r (rigid),
                                    %   normcorre-nr (nonrigid), 
                                    % fft-rigid method (Katie's)
slacknotify = false;              % flag to send Ann slack notification re start and end of processing

% Processing parameters (if not using default)
if ~default
    % motion correction
        % Katie's method
        if strcmpi(mcorr_method,'fftRigid')
            params.mcorr.fftRigid.imscale = 1;             % image downsampling factor                                             [default: 1]
            params.mcorr.fftRigid.Nimg_ave = 10;           % no. of images to be averaged for calculating pixel shift (zippering)  [default: 10]
            params.mcorr.fftRigid.refChannel = 'green';    % channel to be used for calculating image shift (green,red)            [default: 'green']
            params.mcorr.fftRigid.redoT = 300;             % no. of frames at start of file to redo motion correction for after 1st iteration [default: 300]
        end
        % NoRMCorre-rigid
        if strcmpi(mcorr_method,'normcorre-r')
            params.mcorr.normcorre_r = NoRMCorreSetParms(...
                'd1', 512,...
                'd2', 512,...
                'max_shift',30,...          % default: 30
                'bin_width',200,...         % default: 200
                'us_fac',50,...             % default: 50
                'init_batch',200);          % default: 200
        end
        % NoRMCorre-nonrigid
        if strcmpi(mcorr_method,'normcorre-nr')    
            params.mcorr.normcorre_nr = NoRMCorreSetParms(...
                'd1', 512,...
                'd2', 512,...
                'grid_size',[64,64],...     % default: [64,64]
                'overlap_pre',[64,64],...   % default: [64,64]
                'overlap_post',[64,64],...  % default: [64,64]
                'iter',1,...                % default: 1
                'use_parallel',false,...    % default: false
                'max_shift',20,...          % default: 20
                'mot_uf',4,...              % default: 4
                'bin_width',200,...         % default: 200
                'max_dev',3,...             % default: 3
                'us_fac',50,...             % default: 50
                'init_batch',200);          % default: 200
        end
    % roi registration
        params.ROIreg.maxthr = [];                     
        params.ROIreg.dist_maxthr = 0.1;        % threshold for turning spatial components into binary masks [default: 0.1]
        params.ROIreg.dist_exp = 0.8;           % power n for distance between masked components: dist = 1 - (and(m1,m2)/or(m1,m2))^n [default: 1]
        params.ROIreg.dist_thr = 0.7;           % threshold for setting a distance to infinity    [default: 0.5]
        params.ROIreg.dist_overlap_thr = 0.7;   % overlap threshold for detecting if one ROI is a subset of another [default: 0.8]
        params.ROIreg.plot_reg = true;
        params.ROIreg.print_msg = false;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load module folders and define data directory
[data_locn,~,err] = load_neuroSEEmodules(false);
if ~isempty(err)
    beep
    cprintf('Errors',err);    
    return
end

% Default parameters
if default
    C = load('default_params.mat');
    if strcmpi(mcorr_method,'fftRigid')
        params.mcorr.fftRigid = C.mcorr.fftRigid;
    elseif strcmpi(mcorr_method,'normcorre-r')
        params.mcorr.normcorre_r = C.mcorr.normcorre_r;
    elseif strcmpi(mcorr_method,'normcorre-nr')
        params.mcorr.normcorre_nr = C.mcorr.normcorre_nr;
    end
    clear C
end


%% Pre-processing
% Check if data exist for mouseID in env1 and env2. Quit if data does not exist
dir_env1 = [data_locn 'Analysis/' mouseid '/' mouseid '_' env1 env2 '-' env1 ...
           '/group_proc/normcorre-nr_CaImAn_FISSA/' mouseid '_' env1 env2 '-' env1 '_imreg_ref' ref1 '/'];
data_env1 = [dir_env1 mouseid '_' env1 env2 '-' env1 '_ref' ref1 '_PFmap_output.mat'];

dir_env2 = [data_locn 'Analysis/' mouseid '/' mouseid '_' env1 env2 '-' env2 ...
           '/group_proc/normcorre-nr_CaImAn_FISSA/' mouseid '_' env1 env2 '-' env2 '_imreg_ref' ref2 '/'];
data_env2 = [dir_env2 mouseid '_' env1 env2 '-' env2 '_ref' ref2 '_PFmap_output.mat'];
       
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
fdir = [data_locn 'Analysis/' mouseid '/' mouseid '_' env1 env2 ...
           '/remapping/normcorre-nr_CaImAn_FISSA/' mouseid '_' env1 env2 '_imreg_ref' ref1 '-' ref2 '/'];
if ~exist(fdir,'dir'), mkdir(fdir); end
fname_remap = [fdir  mouseid '_' env1 env2 '_remapping_output.mat'];
fname_remapfig = [fdir  mouseid '_' env1 env2 '_remapping_summary.fig'];

if ~exist(fname_remap,'file') || force 
    fprintf('%s: registering ROIs\n',[mouseid '_' env1 env2]);
    M1 = load([dir_env1 mouseid '_' env1 env2 '-' env1 '_ref' ref1 '_segment_output.mat']);
    M2 = load([dir_env2 mouseid '_' env1 env2 '-' env2 '_ref' ref2 '_segment_output.mat']);
    
    PF1 = load([dir_env1 mouseid '_' env1 env2 '-' env1 '_ref' ref1 '_PFmap_output.mat']);
    PF2 = load([dir_env2 mouseid '_' env1 env2 '-' env2 '_ref' ref2 '_PFmap_output.mat']);

    masks1 = M1.masks(:,:,PF1.hist.pcIdx_SIsec);
    masks2 = M2.masks(:,:,PF2.hist.pcIdx_SIsec);
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

    t1 = load([dir_env1 mouseid '_' env1 env2 '-' env1 '_ref' ref1 '_mcorr_template.mat']);
    t2 = load([dir_env2 mouseid '_' env1 env2 '-' env2 '_ref' ref2 '_mcorr_template.mat']);
    template1 = t1.template_g;
    template2 = t2.template_g;

    params.ROIreg_mc = NoRMCorreSetParms(...
                'd1',params.ROIreg.d1,...        % width of image [default: 512]  *Regardless of user-inputted value, neuroSEE_motioncorrect reads this 
                'd2',params.ROIreg.d2,...        % length of image [default: 512] *value from actual image    
                'grid_size',[64,64],...     % default: [64,64]
                'overlap_pre',[64,64],...   % default: [64,64]
                'overlap_post',[64,64],...  % default: [64,64]
                'iter',1,...                % default: 1
                'use_parallel',false,...    % default: false
                'max_shift',20,...          % default: 20
                'mot_uf',4,...              % default: 4
                'bin_width',200,...         % default: 200
                'max_dev',3,...             % default: 3
                'us_fac',50,...             % default: 50
                'init_batch',200);          % default: 200
    params.ROIreg_mc.print_msg = false;

    fname_fig = [fdir  mouseid '_' env1 env2 '_regROIs_output'];
    [matched_ROIs,nonmatched_1,nonmatched_2,A2,R,A_union] = register_ROIs(A1,A2,params.ROIreg,template1,template2,params.ROIreg_mc,fname_fig,true);
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
    
    % PF maps for fam2 in fam1 sorting
    env2PF_env1Sorting = zeros(size([PF1.hist.sort_normpfMap_SIsec_sm]));
    for i = 1:size(masks1,3)
       [matched,ind] = ismember(PF1.hist.sortIdx_SIsec(i),matched_ROIs(:,1));
       if matched
           env2PF_env1Sorting(i,:) = PF2.hist.sort_normpfMap_SIsec_sm(matched_ROIs(ind,2),:);
       end
    end
%     env2PF_env1Sorting = [env2PF_env1Sorting; PF2.hist.normpfMap_SIsec_sm(nonmatched_2,:)];
%     Nbins = numel(PF1.PFdata.occMap);
%     env1PF = [PF1.hist.sort_normpfMap_SIsec_sm; zeros(size(nonmatched_2,2),Nbins)];
%     env2PF = [PF2.hist.sort_normpfMap_SIsec_sm; zeros(size(nonmatched_1,2),Nbins)];
    
    env1PF = PF1.hist.sort_normpfMap_SIsec_sm;
    env2PF = PF2.hist.sort_normpfMap_SIsec_sm;
    
    % save data
    remapping_output.masks = masks_union;
    remapping_output.outlines = outlines;
    remapping_output.matched_ROIs = matched_ROIs;
    remapping_output.nonmatched_env1 = nonmatched_1;
    remapping_output.nonmatched_env2 = nonmatched_2;
    remapping_output.alignment_matrix = R;
    remapping_output.aligned_env2_ROIs = A2;
    remapping_output.masks2_reg = masks2_reg;
    remapping_output.params = params.ROIreg;

    remapping_output.im_env1masks = im_env1masks;
    remapping_output.im_env2masks = im_env2masks;
    remapping_output.im_env1env2masks = im_env1env2masks;
    remapping_output.env1PF = env1PF;
    remapping_output.env2PF_env1Sorting = env2PF_env1Sorting;
    remapping_output.env2PF = env2PF;
    
    fprintf('%s: saving remapping data\n',[mouseid '_' env1 env2]);
    save(fname_remap, '-struct', 'remapping_output')
else
    if ~exist(fname_remapfig,'file') 
        fprintf('%s: loading remapping data\n',[mouseid '_' env1 env2]);
        c = load(fname_remap);
        nonmatched_1 = c.nonmatched_env1;
        nonmatched_2 = c.nonmatched_env2;
        im_env1masks = c.im_env1masks;
        im_env2masks = c.im_env2masks;
        im_env1env2masks = c.im_env1env2masks;
        env1PF = c.env1PF;
        env2PF_env1Sorting = c.env2PF_env1Sorting;
        env2PF = c.env2PF;
        clear c
    end
end

if ~exist(fname_remapfig,'file')
    fh = figure;
    fontsize = 16;
    Nbins = size(env1PF,2);
    subplot(231);
        imshow(im_env1masks); title(env1,'Fontweight','normal','Fontsize',fontsize);    
    subplot(232);
        imshow(im_env2masks); title(env2,'Fontweight','normal','Fontsize',fontsize);
    subplot(233);
        imshow(im_env1env2masks); title([env1 ' \cap ' env2],'Fontweight','normal','Fontsize',fontsize);
    subplot(234);
        cmap = jet(256);
        imagesc(env1PF); 
        colormap(cmap(1:225,:)); %colorbar
        title(env1,'Fontweight','normal','Fontsize',fontsize);
        yticks([]); %yticks([1 159]); yticklabels([159 1]);
        xticks([1 Nbins]); xticklabels([1 100]);
        xlabel('Position (cm)'); %ylabel('Cell no.');
    subplot(235);
        imagesc(env2PF_env1Sorting); 
        title(env2,'Fontweight','normal','Fontsize',fontsize);
        yticks([]);
        xticks([1 Nbins]); xticklabels([1 100]);
        xlabel('Position (cm)');
    subplot(236);
        imagesc(env2PF); 
        title(env2,'Fontweight','normal','Fontsize',fontsize); 
        yticks([]);
        xticks([1 Nbins]); xticklabels([1 100]);
        xlabel('Position (cm)');
        
    fprintf('%s: saving remapping summary figure\n',[mouseid '_' env1 env2]);
    savefig( fh, fname_remapfig(1:end-4) );
    saveas( fh, fname_remapfig(1:end-4), 'png' );
    if figclose, close( fh ); end   
end

    
t = toc;
str = sprintf('%s: Processing done in %g hrs\n', [mouseid '_' env1 env2], round(t/3600,2));
cprintf(str)
end



