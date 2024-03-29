% Written by Ann Go
% INPUTS
%   list:   in the format explist_mXX_exp... .txt
%           *The first experiment on the list must correspond to the reference file  

function frun_ROIreg_multisession( list, ref_array, bl_prctile_array, pfactivet_thr, sessionind_array, force, figclose )

if nargin<7, figclose = true; end
if nargin<6, force = false; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USER-DEFINED INPUT                         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Basic settings
groupreg_method = 'ROIreg';
mcorr_method = 'normcorre';
segment_method = 'CaImAn';
dofissa = true;
    if dofissa, str_fissa = 'FISSA'; else, str_fissa = 'noFISSA'; end

% ROI registration parameters
params.ROIreg.maxthr = [];                     
params.ROIreg.dist_maxthr = 0.1;        % threshold for turning spatial components into binary masks [default: 0.1]
params.ROIreg.dist_exp = 0.5;           % power n for distance between masked components: dist = 1 - (and(m1,m2)/or(m1,m2))^n [default: 1]
params.ROIreg.dist_thr = 0.6;           % threshold for setting a distance to infinity    [default: 0.5]
params.ROIreg.dist_overlap_thr = 0.8;   % overlap threshold for detecting if one ROI is a subset of another [default: 0.8]
params.ROIreg.plot_reg = true;
params.ROIreg.print_msg = false;

% options = neuroSEE_setparams('mcorr_method', mcorr_method, 'dofissa', dofissa); 
g = 64;
options = neuroSEE_setparams(...
            'mcorr_method', mcorr_method, ...
            'segment_method', segment_method,...
            'dofissa', dofissa, ...
            'max_shift_r', 30,...       
            'max_shift_nr', 30,...
            'grid_size_nr', [g,g],...
            'iter',1,...
            'max_dev',8,...     
            'overlap_pre', [g/4,g/4],... 
            'min_patch_size', [g/4,g/4],...      
            'min_diff', [g/8,g/8]);             

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

tic
% Output directory
[ mouseid, expname, fov ] = find_mouseIDexpname( list );
if ~isempty(fov)
    sdir = [data_locn 'Analysis/' mouseid '/' fov '/' mouseid '_' expname '/group_proc/'...
        groupreg_method '_' mcorr_method '_' segment_method '/' str_fissa '/'];
else
    sdir = [data_locn 'Analysis/' mouseid '/' mouseid '_' expname '/group_proc/'...
        groupreg_method '_' mcorr_method '_' segment_method '/' str_fissa '/'];
end
fname_mat_ROIreg = [sdir mouseid '_' expname '_multisessionROIreg_output.mat'];

fname_mat_PFdata = [sdir 'bl' bl_str '_pfTthr' num2str(pfactivet_thr) '/' mouseid '_' expname '_multisessionROIreg_PFdata.mat'];
figdir = [sdir 'bl' bl_str '_pfTthr' num2str(pfactivet_thr) '/'];

%% Load individual session data
if any([ ~exist(fname_mat_ROIreg,'file'), ~exist(fname_mat_PFdata,'file'), force ])
    fprintf('%s: loading individual session PF data\n',[mouseid '_' expname]);
    listfile = [data_locn 'Digital_Logbook/lists_imaging/' list];
    exps = extractExpnamesFromTxtfile( listfile );
    Nexps = numel(exps);

    % initialise variables
    for n = 1:Nexps
        [ mouseid_n, exp_n ] = find_mouseIDexpname( exps{n} );
        if ~strcmpi(mouseid_n, mouseid)
            beep
            cprintf('Errors','Invalid list. Experiments must all be for a single animal.');   
            return
        end
        if ~isempty(fov)
            fname_pref = [data_locn 'Analysis/' mouseid '/' fov '/' mouseid '_' exp_n ...
                       '/group_proc/imreg_' mcorr_method '_' segment_method '/' ...
                       mouseid '_' exp_n '_imreg_ref' ref_array(n,:) '/' ...
                       mouseid '_' exp_n '_ref' ref_array(n,:)];
        else 
            fname_pref = [data_locn 'Analysis/' mouseid '/' mouseid '_' exp_n ...
                       '/group_proc/imreg_' mcorr_method '_' segment_method '/' ...
                       mouseid '_' exp_n '_imreg_ref' ref_array(n,:) '/' ...
                       mouseid '_' exp_n '_ref' ref_array(n,:)];
        end
        
        % load roi segmentation output
        M = load([fname_pref '_segment_output.mat']);
        masks{n} = M.masks;
        AA = zeros(size(masks{n},1)*size(masks{n},2),size(masks{n},3));
        for i = 1:size(masks{n},3)
            mask = masks{n}(:,:,i);
            AA(:,i) = mask(:);
        end
        A{n} = sparse(AA); 
        
        % load motion correction templates
        t = load([fname_pref '_imreg_template.mat']);
        templates{n} = t.template_g;
        % templatesR{n} = t.template_r;
        
        % load pf mapping data
        if ~isempty(fov)
            fname = [data_locn 'Analysis/' mouseid '/' fov '/' mouseid '_' exp_n ...
                       '/group_proc/imreg_' mcorr_method '_' segment_method '/' ...
                       mouseid '_' exp_n '_imreg_ref' ref_array(n,:) '/' str_fissa ...
                       '/bl_prctile' num2str(bl_prctile_array(n)) '/' ...
                       mouseid '_' exp_n '_ref' ref_array(n,:) '_PFmap_output.mat'];
        else
            fname = [data_locn 'Analysis/' mouseid '/' mouseid '_' exp_n ...
                       '/group_proc/imreg_' mcorr_method '_' segment_method '/' ...
                       mouseid '_' exp_n '_imreg_ref' ref_array(n,:) '/' str_fissa ...
                       '/bl_prctile' num2str(bl_prctile_array(n)) '/' ...
                       mouseid '_' exp_n '_ref' ref_array(n,:) '_PFmap_output.mat'];
        end
        
        pfMat = load(fname);
        normrMap_sm{n} = pfMat.hist.normrMap_sm;
        if size(pfMat.hist.normrMap_sm,3) == 1
            params.mode_dim = '1D';    
            normrMap_sm{n} = pfMat.hist.normrMap_sm;
            sortpcIdx{n} = pfMat.hist.SIspk.sortpcIdx;
            pfLoc{n} = pfMat.hist.pfLoc;
            try
                normspkRaster{n} = pfMat.pfData.normspkRaster;
                ytick_files{n} = pfMat.pfData.ytick_files;
            catch
                normspkRaster{n} = pfMat.PFdata.normspkRaster;
                ytick_files{n} = pfMat.PFdata.ytick_files;
            end
        else
            params.mode_dim = '2D';    
            rMap_sm{n} = pfMat.hist.normrMap_sm;
            normrMap_sm{n} = pfMat.hist.normrMap_sm;
            pfLoc{n} = pfMat.hist.centroid;
            activeData{n} = pfMat.activeData;
        end
        
        % PCs only
        pcIdx{n} = pfMat.hist.SIspk.pcIdx;
        masks_pcs{n} = M.masks(:,:,pcIdx{n});
        AA = zeros(size(masks_pcs{n},1)*size(masks_pcs{n},2),size(masks_pcs{n},3));
        for i = 1:size(masks_pcs{n},3)
            mask = masks_pcs{n}(:,:,i);
            AA(:,i) = mask(:);
        end
        A_pcs{n} = sparse(AA);
        
    end
    
    
    if ~exist(fname_mat_ROIreg,'file') || force
        %% Multisession registration
        params.ROIreg.d1 = size(masks{n},1);
        params.ROIreg.d2 = size(masks{n},2);
        params.ROIreg.plot_reg = true;
        params.ROIreg.print_msg = true;

        fprintf('%s: registering ROIs multisession\n',[mouseid '_' expname]);
        [A_union, A_shifted, assignments, matchings, templates_shifted, matched_ROIs, nonmatched_1, nonmatched_2, R ] = ...
            register_multisession_singleref(A, params.ROIreg, templates, params.ROIreg_mc, [], false);
        for nn = 1:Nexps-1
            masks_shifted{nn} = reshape(full(A_shifted{nn}), params.ROIreg.d1, params.ROIreg.d2, size(A_shifted{nn},2));
            masks_union{nn} = reshape(full(A_union{nn+1}), params.ROIreg.d1, params.ROIreg.d2, size(A_union{nn+1},2));
        end


        %% Multisession registration for PCs only
        params.ROIreg.d1 = size(masks_pcs{n},1);
        params.ROIreg.d2 = size(masks_pcs{n},2);
        params.ROIreg.plot_reg = true;
        params.ROIreg.print_msg = true;

        [~, A_shifted_pcs, assignments_pcs, matchings_pcs, ~, matched_ROIs_pcs, nonmatched_1_pcs, nonmatched_2_pcs, ~, A_union_pcs ] = ...
            register_multisession(A_pcs, params.ROIreg, templates, params.ROIreg_mc, [], false);
        for nn = 1:Nexps-1
            masks_shifted_pcs{nn} = reshape(full(A_shifted_pcs{nn}), params.ROIreg.d1, params.ROIreg.d2, size(A_shifted_pcs{nn},2));
            masks_union_pcs{nn} = reshape(full(A_union_pcs{nn}), params.ROIreg.d1, params.ROIreg.d2, size(A_union_pcs{nn},2));
        end
    else
        fprintf('%s: loading multisession ROI registration\n',[mouseid '_' expname]);
        c1 = load(fname_mat_ROIreg);
        templates = c1.templates;
        templates_shifted = c1.templates_shifted;
        matched_ROIs = c1.matched_ROIs;
        nonmatched_1 = c1.nonmatched_1;
        nonmatched_2 = c1.nonmatched_2;
        masks_shifted = c1.masks_shifted;
            for n = 1:size(masks_shifted,2)
            AA = zeros(size(masks_shifted{n},1)*size(masks_shifted{n},2),size(masks_shifted{n},3));
                for i = 1:size(masks_shifted{n},3)
                    mask = masks_shifted{n}(:,:,i);
                    AA(:,i) = mask(:);
                end
                A_shifted{n} = sparse(AA); 
            end
        matched_ROIs_pcs = c1.matched_ROIs_pcs;
        nonmatched_1_pcs = c1.nonmatched_1_pcs;
        nonmatched_2_pcs = c1.nonmatched_2_pcs;
        masks = c1.masks;
            for n = 1:size(masks,2)
                AA = zeros(size(masks{n},1)*size(masks{n},2),size(masks{n},3));
                for i = 1:size(masks{n},3)
                    mask = masks{n}(:,:,i);
                    AA(:,i) = mask(:);
                end
                A{n} = sparse(AA); 
            end
        masks_pcs = c1.masks_pcs;
            for n = 1:size(masks_pcs,2)
                AA = zeros(size(masks_pcs{n},1)*size(masks_pcs{n},2),size(masks_pcs{n},3));
                for i = 1:size(masks_pcs{n},3)
                    mask = masks_pcs{n}(:,:,i);
                    AA(:,i) = mask(:);
                end
                A_pcs{n} = sparse(AA); 
            end
        masks_union = c1.masks_union;
            for n = 1:size(masks_union,2)
                AA = zeros(size(masks_union{n},1)*size(masks_union{n},2),size(masks_union{n},3));
                for i = 1:size(masks_union{n},3)
                    mask = masks_union{n}(:,:,i);
                    AA(:,i) = mask(:);
                end
                A_union{n} = sparse(AA); 
            end
        masks_union_pcs = c1.masks_union_pcs;
            for n = 1:size(masks_union_pcs,2)
                AA = zeros(size(masks_union_pcs{n},1)*size(masks_union_pcs{n},2),size(masks_union_pcs{n},3));
                for i = 1:size(masks_union_pcs{n},3)
                    mask = masks_union_pcs{n}(:,:,i);
                    AA(:,i) = mask(:);
                end
                A_union_pcs{n} = sparse(AA); 
            end
        masks_shifted_pcs = c1.masks_shifted_pcs;
            for n = 1:size(masks_shifted_pcs,2)
                AA = zeros(size(masks_shifted_pcs{n},1)*size(masks_shifted_pcs{n},2),size(masks_shifted_pcs{n},3));
                for i = 1:size(masks_shifted_pcs{n},3)
                    mask = masks_shifted_pcs{n}(:,:,i);
                    AA(:,i) = mask(:);
                end
                A_shifted_pcs{n} = sparse(AA); 
            end
        assignments = c1.assignments;
        matchings = c1.matchings;
    end    
    
    % Save output
    % ROIreg output
    output1.templates = templates;
    output1.templates_shifted = templates_shifted;
    output1.params = params;
%    output1.A_shifted = A_shifted; % exclude, output file becomes too large when included
%    output1.A_union = A_union;
    output1.matched_ROIs = matched_ROIs;
    output1.nonmatched_1 = nonmatched_1;
    output1.nonmatched_2 = nonmatched_2;
    output1.masks_shifted = masks_shifted;
    
    output1.matched_ROIs_pcs = matched_ROIs_pcs;
    output1.nonmatched_1_pcs = nonmatched_1_pcs;
    output1.nonmatched_2_pcs = nonmatched_2_pcs;
    output1.masks = masks;
    output1.masks_pcs = masks_pcs;
    output1.masks_union = masks_union;
    output1.masks_union_pcs = masks_union_pcs;
    output1.masks_shifted_pcs = masks_shifted_pcs;
    output1.assignments = assignments;
    % output1.assignments_pcs = assignments_pcs;
    output1.matchings = matchings;
    % output1.matchings_pcs = matchings_pcs;
    if ~exist(sdir,'dir'); mkdir(sdir); fileattrib(sdir,'+w','g','s'); end
    save(fname_mat_ROIreg,'-struct','output1');

    % PF data
    output2.rMap_sm = rMap_sm;
    output2.normrMap_sm = normrMap_sm;
    output2.pfLoc = pfLoc;
    output2.pcIdx = pcIdx;
    output2.activeData = activeData;
    if strcmp(params.mode_dim,'1D')
        output2.normspkRaster = normspkRaster;
        output2.ytick_files = ytick_files;
        output2.sortpcIdx = sortpcIdx;
    end
    fprintf('%s: saving multisession roi registration output\n',[mouseid '_' expname]);
    if ~exist(figdir,'dir'), mkdir(figdir); fileattrib(figdir,'+w','g','s'); end
    save(fname_mat_PFdata,'-struct','output2');
        
else
    % load ROIreg data
    c1 = load(fname_mat_ROIreg);
    templates = c1.templates;
    templates_shifted = c1.templates_shifted;
    matched_ROIs = c1.matched_ROIs;
    nonmatched_1 = c1.nonmatched_1;
    nonmatched_2 = c1.nonmatched_2;
    masks_shifted = c1.masks_shifted;
        for n = 1:size(masks_shifted,2)
        AA = zeros(size(masks_shifted{n},1)*size(masks_shifted{n},2),size(masks_shifted{n},3));
            for i = 1:size(masks_shifted{n},3)
                mask = masks_shifted{n}(:,:,i);
                AA(:,i) = mask(:);
            end
            A_shifted{n} = sparse(AA); 
        end
    matched_ROIs_pcs = c1.matched_ROIs_pcs;
    nonmatched_1_pcs = c1.nonmatched_1_pcs;
    nonmatched_2_pcs = c1.nonmatched_2_pcs;
    masks = c1.masks;
        for n = 1:size(masks,2)
            AA = zeros(size(masks{n},1)*size(masks{n},2),size(masks{n},3));
            for i = 1:size(masks{n},3)
                mask = masks{n}(:,:,i);
                AA(:,i) = mask(:);
            end
            A{n} = sparse(AA); 
        end
    masks_pcs = c1.masks_pcs;
        for n = 1:size(masks_pcs,2)
            AA = zeros(size(masks_pcs{n},1)*size(masks_pcs{n},2),size(masks_pcs{n},3));
            for i = 1:size(masks_pcs{n},3)
                mask = masks_pcs{n}(:,:,i);
                AA(:,i) = mask(:);
            end
            A_pcs{n} = sparse(AA); 
        end
    masks_union = c1.masks_union;
        for n = 1:size(masks_union,2)
            AA = zeros(size(masks_union{n},1)*size(masks_union{n},2),size(masks_union{n},3));
            for i = 1:size(masks_union{n},3)
                mask = masks_union{n}(:,:,i);
                AA(:,i) = mask(:);
            end
            A_union{n} = sparse(AA); 
        end
    masks_union_pcs = c1.masks_union_pcs;
        for n = 1:size(masks_union_pcs,2)
            AA = zeros(size(masks_union_pcs{n},1)*size(masks_union_pcs{n},2),size(masks_union_pcs{n},3));
            for i = 1:size(masks_union_pcs{n},3)
                mask = masks_union_pcs{n}(:,:,i);
                AA(:,i) = mask(:);
            end
            A_union_pcs{n} = sparse(AA); 
        end
    masks_shifted_pcs = c1.masks_shifted_pcs;
        for n = 1:size(masks_shifted_pcs,2)
            AA = zeros(size(masks_shifted_pcs{n},1)*size(masks_shifted_pcs{n},2),size(masks_shifted_pcs{n},3));
            for i = 1:size(masks_shifted_pcs{n},3)
                mask = masks_shifted_pcs{n}(:,:,i);
                AA(:,i) = mask(:);
            end
            A_shifted_pcs{n} = sparse(AA); 
        end
    assignments = c1.assignments;
    matchings = c1.matchings;
    % assignments_pcs = c2.assignments_pcs;
    % matchingspcs = c2.matchings_pcs;

    % load data for paper plots
    c2 = load(fname_mat_PFdata);
    normrMap_sm = c2.normrMap_sm;
    if size(normrMap_sm{1},3) == 1
        params.mode_dim = '1D';   
        normspkRaster = c2.normspkRaster;
        ytick_files = c2.ytick_files;
        sortpcIdx = c2.sortpcIdx;
    else
        params.mode_dim = '2D';
        rMap_sm = c2.rMap_sm;
        activeData = c2.activeData;
    end
    pfLoc = c2.pfLoc;
    pcIdx = c2.pcIdx;
      
end

%% Report plots
% comparison of registered templates   
Nexps = size(A,2); 
nCol = Nexps-1;
if ~exist([sdir mouseid '_' expname '_regtemplates.fig'],'file') || force
    fh1 = figure;
    for n = 1:Nexps-1
        subplot(3,nCol,n); 
            imshow(imfuse( templates{1}, templates{n+1}, 'falsecolor', 'Scaling', 'joint', 'ColorChannels', [2 1 0]) );
            str = sprintf('%g to %g: Before reg', sessionind_array(n+1), sessionind_array(1));
            title( str );

        subplot(3,nCol,nCol+n); 
            imshow( imfuse( templates{1}, templates_shifted{n}, 'falsecolor', 'Scaling', 'joint', 'ColorChannels', [2 1 0]) );
            str = sprintf('%g to %g: After reg', sessionind_array(n+1), sessionind_array(1));
            title( str );

        subplot(3,nCol,2*nCol+n); 
            imshow( imfuse( templates{1}, templates_shifted{n}, 'falsecolor', 'Scaling', 'joint', 'ColorChannels', [1 2 0]) );
            str = sprintf('Reversed colours');
            title( str );
    end
    fprintf('%s: saving multisession template registration summary figure\n',[mouseid '_' expname]);
    savefig( fh1, [sdir mouseid '_' expname '_regtemplates'] );
    saveas( fh1, [sdir mouseid '_' expname '_regtemplates'], 'png' );
    if figclose, close( fh1 ); end   
end

% comparison of registered rois across sessions
if ~exist([sdir mouseid '_' expname '_regROIs'],'file') || force
    fh2 = figure;
    for n = 1:Nexps-1
        subplot(2,nCol,n);
            imagesc(templates{1}); colormap(gray); axis off; axis square;
            options.plot_bck_image = false;
            plot_contours( A{n}(:,matched_ROIs{n}(:,1)), templates{1}, options,0,[],[],'m'); hold on;
            plot_contours( A_shifted{n}(:,matched_ROIs{n}(:,2)), templates{1}, options,0,[],[],'c'); hold on;
            plot_contours( A{n}(:,nonmatched_1{n}), templates{1}, options, 0,[],[],'r'); hold on;
            plot_contours( A_shifted{n}(:,nonmatched_2{n}), templates{1}, options, 0,[],[],'b'); hold on;
            h = zeros(4, 1);
            h(1) = plot(NaN,NaN,'m');
            h(2) = plot(NaN,NaN,'c');
            h(3) = plot(NaN,NaN,'r');
            h(4) = plot(NaN,NaN,'b');
            legend(h, ['M#' num2str(sessionind_array(n+1))], ['M#' sessionind_array(1)],...
                ['NoM#' num2str(sessionind_array(n+1))], ['NoM#' sessionind_array(1)]); hold off;
            str = sprintf('%g aligned to %g', sessionind_array(n+1), sessionind_array(1));
            title( str );
        subplot(2,nCol,nCol+n);
            imagesc(templates{n+1}); colormap(gray); axis off; axis square;
            options.plot_bck_image = false;
            plot_contours( A_union{n}, templates{n+1}, options,0,[],[],'w'); hold on;
            num = num2str(1:n); num = num(~isspace(num));
            title( [sessionind_array(1) '\cup' sessionind_array(n+1)] );
    end    
    fprintf('%s: saving multisession roi registration summary figure\n',[mouseid '_' expname]);
    savefig( fh2, [sdir mouseid '_' expname '_regROIs'] );
    saveas( fh2, [sdir mouseid '_' expname '_regROIs'], 'png' );
    if figclose, close( fh2 ); end 
end

% comparison of registered rois across sessions for PLACE CELLS only
%     fh3 = figure;
%     for n = 1:Nexps-1
%         subplot(2,nCol,n);
%             imagesc(templates{n+1}); colormap(gray); axis off; axis square;
%             options.plot_bck_image = false;
%             plot_contours( A_pcs{n}(:,matched_ROIs_pcs{n}(:,1)), templates{1}, options,0,[],[],'m'); hold on;
%             plot_contours( A_shifted_pcs{n}(:,matched_ROIs_pcs{n}(:,2)), templates{1}, options,0,[],[],'c'); hold on;
%             plot_contours( A_pcs{n}(:,nonmatched_1_pcs{n}), templates{n+1}, options, 0,[],[],'r'); hold on;
%             plot_contours( A_shifted_pcs{n}(:,nonmatched_2_pcs{n}), templates{1}, options, 0,[],[],'b'); hold on;
%             h = zeros(4, 1);
%             h(1) = plot(NaN,NaN,'m');
%             h(2) = plot(NaN,NaN,'c');
%             h(3) = plot(NaN,NaN,'r');
%             h(4) = plot(NaN,NaN,'b');
%             legend(h, ['M#' num2str(sessionind_array(n+1))], ['M#' sessionind_array(1)],...
%                 ['NoM#' num2str(sessionind_array(n+1))], ['NoM#' sessionind_array(1)]); hold off;
%             str = sprintf('%g aligned to %g', sessionind_array(n+1), sessionind_array(1));
%             title( str );
%         subplot(2,nCol,nCol+n);
%             imagesc(templates{n+1}); colormap(gray); axis off; axis square;
%             options.plot_bck_image = false;
%             plot_contours( A_union_pcs{n}, templates{n+1}, options,0,[],[],'w'); hold on;
%             title( [sessionind_array(1) '\cup' sessionind_array(n+1)] );
%     end    
%     fprintf('%s: saving multisession roi registration summary figure\n',[mouseid '_' expname]);
%     savefig( fh3, [figdir mouseid '_' expname '_regROIs_PCs'] );
%     saveas( fh3, [figdir mouseid '_' expname '_regROIs_PCs'], 'png' );
%     if figclose, close( fh3 ); end 

% spike trial raster plot across sessions
Ncells = size(A_union{end},2);
if strcmp(params.mode_dim,'1D')
    if ~exist([figdir 'rasterplots/'],'dir') || force
        nRow = 10;
        nCol = 4;
        cmap0 = [0.9 0.9 0.9];
        cmap1 = [0 0 1];
        cmap = zeros(50,3);
        for j=1:3
            cmap(:,j) = linspace(cmap0(j),cmap1(j),50);
        end
        colormap(cmap);

        for ii=0:ceil(Ncells/nRow)-1 
            fh = figure('Position',[680 316 535 782]); 
            ha = tight_subplot(nRow,nCol,[.05 0.04],[.03 .05],[.06 .02]);
            for jj=0:nRow-1
                if (ii*nRow+jj+1) <= Ncells
                    for n = 1:Nexps
                        axes(ha(jj*nCol+n));
                        idx = assignments(ii*nRow+jj+1, n);
                        if n == 1
                            ylabel(['Cell ' num2str(ii*nRow+jj+1)]);
                        end
                        if ~isnan(idx)
                            imagesc( normspkRaster{n}{idx} ); colormap(cmap);
                            yticks([1 ytick_files{n}(end)]); %yticklabels(ytick_files{n}); 
                            xticks([]); 
                            if n == 1
                                ylabel(['Cell ' num2str(ii*nRow+jj+1)]);
                            end
                            if numel(find(pcIdx{n} == idx)) > 0
                                title_str = ['PC ' num2str(idx)];
                            else
                                title_str = ['nonPC ' num2str(idx)];
                            end
                            title( title_str, 'fontsize', 10);
                        end
                    end
                end
            end
            if ii == 0
                fprintf('%s: saving spike trial raster plots\n',[mouseid '_' expname]);
                if ~exist([figdir 'rasterplots/'],'dir')
                    mkdir([figdir 'rasterplots/']); 
                    fileattrib([figdir 'rasterplots/'],'+w','g','s');
                end
            end
            savefig( fh, [figdir 'rasterplots/' mouseid '_' expname '_normspkRaster_' num2str(ii+1)] );
            saveas( fh, [figdir 'rasterplots/' mouseid '_' expname '_normspkRaster_' num2str(ii+1)], 'png' );
            if figclose, close( fh ); end 
        end
    end
    
    % place tuning across sessions
    if ~exist([figdir mouseid '_' expname '_sortPFmaps.fig'],'file') || force
        Nbins = size(normrMap_sm{1},2);
        for n = 1:Nexps
            for j = 1:Nexps
                if n ~= j
                    pf_sort{n}{j} = zeros(numel(sortpcIdx{n}),Nbins);
                    for i = 1:size(sortpcIdx{n},1)
                        iIdx = sortpcIdx{n}(i);
                        kIdx = matchings{n}(iIdx);
                        jIdx = assignments(kIdx,j);
                        if ~isnan(jIdx)
                            pf_sort{n}{j}(i,:) = normrMap_sm{j}(jIdx,:);
                        end
                    end
                else
                   pf_sort{n}{j} = normrMap_sm{n}(sortpcIdx{n},:); 
                end
            end
        end

        fh4 = figure('Position', [680 547 446 551]);
        nRow = Nexps; nCol = Nexps; 
        for n = 0:Nexps-1
            for j = 0:Nexps-1
                if n~=j, k = 1; else, k = 15; end
                cmap = viridisMap;
                ax = subplot(nRow,nCol,n*nRow+j+1);
                imagesc( pf_sort{n+1}{j+1} );
                colormap(ax,cmap(k:end,:));
                xticks([]); yticks([1 size(pf_sort{n+1}{j+1},1)]); 
                if n*nRow+j+1 <= nCol
                    title(['Session ' num2str(j+1)], 'fontsize', 10);
                end
                if j==0
                    ylabel(['S' num2str(n+1) ' sorting']);
                end
            end
        end
        fprintf('%s: saving place tuning across sessions summary figure\n',[mouseid '_' expname]);
        savefig( fh4, [figdir mouseid '_' expname '_sortPFmaps'] );
        saveas( fh4, [figdir mouseid '_' expname '_sortPFmaps'], 'png' );
        if figclose, close( fh4 ); end 
    end
else % 2D
    % spike locations, place field maps across sessions
    if ~exist([figdir 'pfmaps/'],'dir') || force
        nRow = 10;
        nCol = 2*Nexps;

        for ii=0:ceil(Ncells/nRow)-1 
            fh = figure('Position',[680 678 500 600]); 
            ha = tight_subplot(nRow,nCol,[.04 .02],[.02 .05],[.02 .02]);

            for jj=0:nRow-1
                c = ii*nRow+jj+1;
                if c <= Ncells
                    for s = 1:Nexps
                        idx = assignments(ii*nRow+jj+1, s);
                        % find the delineations for the video: find t = 0
                        idx_file = find(diff(activeData{s}.t) < 0);
                        idx_file = [0; idx_file; numel(activeData{s}.t)] +1;
                        max_spksize = 18;
                        spk_col = 'b';

                        axes(ha(jj*nCol+1+(s-1)*2));
                        hold on; axis([-140 140 -140 140]); 
                        if s >1, axis off; end
                        % no mistake here in the order of plotting x&y, this is
                        % to match the image pixel indexing to matrix row and
                        % column indexing
                        for kk = 1: numel(idx_file) - 1
                            plot(activeData{s}.y(idx_file(kk):idx_file(kk+1)-1),-activeData{s}.x(idx_file(kk):idx_file(kk+1)-1),...
                                'Color',[0.8 0.8 0.8], 'LineWidth',1); axis square; 
                            if s == 1
                                ylabel(['Cell ' num2str(c)]);
                            end
                        end
                        if ~isnan(idx)
                            z = activeData{s}.spikes(idx,:);
                            ind = find(z>0);
                            if ~isempty(ind)
                                x = activeData{s}.x(ind);
                                y = activeData{s}.y(ind);
                                spikes = z(ind);
                                [spkampl_sorted,sort_ind] = sort(spikes);
                                spk_size = spkampl_sorted/max(spkampl_sorted)*max_spksize;
                                scatter(y(sort_ind),-x(sort_ind),spk_size,spk_col,'filled');
                                title_str = sprintf('Cell %g', c); 
                                title(title_str,'fontsize',12);
                                hold off; 
                            end

                            if numel(find(pcIdx{s} == c)) > 0
                                title_str = 'PC';
                            else
                                title_str = 'nonPC';
                            end
                            title( title_str, 'fontsize', 10);

                            axes(ha(jj*nCol+(s-1)*2+2));
                            cmap = viridisMap_whitelowest;
                            colormap(cmap);
                            imagesc(squeeze(normrMap_sm{s}(:,:,idx))');
                            axis off; axis square; % colorbar; 
                            %title_str = sprintf('Max %.2f events/s', max(max(rMap_sm{idx}(:,:,c)))); 
                            %title(title_str,'fontsize',12);
                        end
                    end
                end
            end
            if ii == 0
                fprintf('%s: saving multisession firing locations and pf maps\n',[mouseid '_' expname]);
                if ~exist([figdir 'pfmaps/'],'dir')
                    mkdir([figdir 'pfmaps/']); 
                    fileattrib([figdir 'pfmaps/'],'+w','g','s'); 
                end
            end
            savefig( fh, [figdir 'pfmaps/' mouseid '_' expname '_pfmaps_' num2str(ii+1)] );
            saveas( fh, [figdir 'pfmaps/' mouseid '_' expname '_pfmaps_' num2str(ii+1)], 'png' );
            if figclose, close( fh ); end 
        end
    end
end

% no. of active cells, pcs, stability across sessions
if ~exist([figdir mouseid '_' expname '_stability.fig'],'file') || force
    Nactivecells = zeros(Nexps,1);
    Npcs = zeros(Nexps,1);
    for n = 1:Nexps
        Nactivecells(n) = 100*(size(A{n},2)/size(A_union{end},2));
        Npcs(n) = 100*(numel(pcIdx{n})/size(A{n},2));        
    end
    fh5 = figure('Position', [1340 932 560 162]);
    subplot(131); plot(1:Nexps,Nactivecells); xlabel('Session no.'); ylabel('Active cells (%)');
        hold on; plot(1:Nexps,Nactivecells,'b.','markersize',12); 
        hold off; axis([1 4 0 100]); box off;
    subplot(132); plot(1:Nexps,Npcs); xlabel('Session no.'); ylabel('No. of PCs (%)');
        hold on; plot(1:Nexps,Npcs,'b.','markersize',12); 
        hold off; axis([1 sessionind_array(end) 1 100]); box off;
    subplot(133); xlabel('Elapsed time'); ylabel('Stability');

    fprintf('%s: saving population summary on stability\n',[mouseid '_' expname]);
    savefig( fh5, [figdir mouseid '_' expname '_stability'] );
    saveas( fh5, [figdir mouseid '_' expname '_stability'], 'png' );
    if figclose, close( fh5 ); end 
end
    
t = toc;
str = sprintf('%s: Processing done in %g hrs\n', [mouseid '_' expname], round(t/3600,2));
cprintf(str)
end




