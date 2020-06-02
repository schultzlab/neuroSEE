function frun_ROIreg_multisession( list, ref_array, force, figclose )

if nargin<4, figclose = true; end
if nargin<3, force = false; end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USER-DEFINED INPUT                         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Basic settings
groupreg_method = 'imreg';
imreg_method = 'normcorre';
mcorr_method = 'normcorre';
segment_method = 'CaImAn';
dofissa = true;
    if dofissa, str_fissa = 'FISSA'; else, str_fissa = 'noFISSA'; end

% ROI registration parameters
params.ROIreg.maxthr = [];                     
params.ROIreg.dist_maxthr = 0.1;        % threshold for turning spatial components into binary masks [default: 0.1]
params.ROIreg.dist_exp = 0.7;           % power n for distance between masked components: dist = 1 - (and(m1,m2)/or(m1,m2))^n [default: 1]
params.ROIreg.dist_thr = 0.7;           % threshold for setting a distance to infinity    [default: 0.5]
params.ROIreg.dist_overlap_thr = 0.6;   % overlap threshold for detecting if one ROI is a subset of another [default: 0.8]
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

tic
% Output directory
[ mouseid, expname ] = find_mouseIDexpname( list );
sdir = [data_locn 'Analysis/' mouseid '/' mouseid '_' expname '/group_proc/'...
        groupreg_method '_' imreg_method '_' segment_method '_' str_fissa '/'];

fname_mat1 = [sdir mouseid '_' expname '_multisessionROIreg_output1.mat'];
fname_mat2 = [sdir mouseid '_' expname '_multisessionROIreg_output2.mat'];
figdir = [sdir 'multisession_pfmaps/'];

%% Load individual session data
if any([ ~exist(fname_mat1,'file'), ~exist(fname_mat2,'file'), force ])
    fprintf('%s: registering ROIs multisession\n',[mouseid '_' expname]);
    listfile = [data_locn 'Digital_Logbook/lists/' list];
    exps = extractExpnamesFromTxtfile( listfile );
    Nexps = numel(exps);

    % initialise variables
    masks{1:Nexps} = [];
    A{1:Nexps} = [];
    templates{1:Nexps} = [];
    normrMap_sm{1:Nexps} = [];
    normpfMap_sm{1:Nexps} = [];
    sortIdx{1:Nexps} = [];
    pfLoc{1:Nexps} = [];
    normspkRaster{1:Nexps} = [];
    ytick_files{1:Nexps} = [];
    pcIdx{1:Nexps} = [];
    masks_pcs{1:Nexps} = [];
    A_pcs{1:Nexps} = [];
    for n = 1:Nexps
        [ mouseid_n, exp_n ] = find_mouseIDexpname( exps{n} );
        if ~strcmpi(mouseid_n, mouseid)
            beep
            cprintf('Errors','Invalid list. Not all experiments for a single animal.');   
            return
        end
        if strcmpi(imreg_method, mcorr_method)
            fname_pref = [data_locn 'Analysis/' mouseid '/' mouseid '_' exp_n ...
                       '/group_proc/imreg_' imreg_method '_' segment_method '_' str_fissa '/' ...
                       mouseid '_' exp_n '_imreg_ref' ref_array(n,:) '/' ...
                       mouseid '_' exp_n '_ref' ref_array(n,:)];
        else
            fname_pref = [data_locn 'Analysis/' mouseid '/' mouseid '_' exp_n ...
                       '/group_proc/imreg_' imreg_method '_' segment_method '_' str_fissa '/' ...
                       mouseid '_' exp_n '_imreg_ref' ref_array(n,:) '_' mcorr_method '/' ...
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
        t = load([fname_pref '_mcorr_template.mat']);
        templates{n} = t.template_g;
        % templatesR{n} = t.template_r;
        
        % load pf mapping data
        pfdata = load([fname_pref '_PFmap_output.mat']);
        normrMap_sm{n} = pfdata.hist.normrMap_sm;
        normpfMap_sm{n} = pfdata.hist.SIsec.normpfMap_sm;
        sortIdx{n} = pfdata.hist.SIsec.sortIdx;
        pfLoc{n} = pfdata.hist.maxLoc;
        normspkRaster{n} = pfdata.pfData.normspkRaster;
        ytick_files{n} = pfdata.pfData.ytick_files;
        
        % PCs only
        pcIdx{n} = pfdata.hist.SIsec.pcIdx;
        masks_pcs{n} = M.masks(:,:,pcIdx{n});
        AA = zeros(size(masks_pcs{n},1)*size(masks_pcs{n},2),size(masks_pcs{n},3));
        for i = 1:size(masks_pcs{n},3)
            mask = masks_pcs{n}(:,:,i);
            AA(:,i) = mask(:);
        end
        A_pcs{n} = sparse(AA);
        
    end
    
    
    %% Choose reference template for multisession registration
%     figure;
%     [nRow, nCol] = getnRownCol(Nexps);
%     for n = 1:Nexps
%         subplot(2*nRow,nCol,n);
%             imagesc(templates{n}); colormap(gray); axis off; axis square;
%             title(['Session ' num2str(n)]);
%         subplot(2*nRow,nCol,nCol+n);
%             imagesc(templates{n}); colormap(gray); axis off; axis square;
%             options.plot_bck_image = false;
%             plot_contours( A{n}, templates{n}, options,0,[],[],'w'); hold on;
%     end


    %% Multisession registration
    params.ROIreg.d1 = size(masks{n},1);
    params.ROIreg.d2 = size(masks{n},2);
    params.ROIreg.plot_reg = true;
    params.ROIreg.print_msg = true;

    [~, A_shifted, assignments, matchings, templates_shifted, matched_ROIs, nonmatched_1, nonmatched_2, ~, A_union ] = ...
        register_multisession(A, params.ROIreg, templates, params.ROIreg_mc, [], false);
    masks_shifted{1:Nexps-1} = [];
    masks_union{1:Nexps-1} = [];
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
    masks_shifted_pcs{1:Nexps-1} = [];
    masks_union_pcs{1:Nexps-1} = [];
    for nn = 1:Nexps-1
        masks_shifted_pcs{nn} = reshape(full(A_shifted_pcs{nn}), params.ROIreg.d1, params.ROIreg.d2, size(A_shifted_pcs{nn},2));
        masks_union_pcs{nn} = reshape(full(A_union_pcs{nn}), params.ROIreg.d1, params.ROIreg.d2, size(A_union_pcs{nn},2));
    end
    
    
    % Save output
    if ~exist(sdir,'dir'), mkdir(sdir); end
    % output needed for paper plots
    output1.normrMap_sm = normrMap_sm;
    output1.normpfMap_sm = normpfMap_sm;
    output1.sortIdx = sortIdx;
    output1.pfLoc = pfLoc;
    output1.normspkRaster = normspkRaster;
    output1.ytick_files = ytick_files;
    output1.pcIdx = pcIdx;
    output1.masks = masks;
    output1.masks_union = masks_union;
    output1.assignments = assignments;
    output1.matchings = matchings;
    fprintf('%s: saving multisession roi registration output\n',[mouseid '_' expname]);
    save(fname_mat1,'-struct','output1');
    
    % output needed for other report plots
    output2.templates = templates;
    output2.templates_shifted = templates_shifted;
    output2.params = params;
%    output2.A_shifted = A_shifted; % exclude, output file becomes too large when included
%    output2.A_union = A_union;
    output2.matched_ROIs = matched_ROIs;
    output2.nonmatched_1 = nonmatched_1;
    output2.nonmatched_2 = nonmatched_2;
    output2.masks_shifted = masks_shifted;
    
    output2.matched_ROIs_pcs = matched_ROIs_pcs;
    output2.nonmatched_1_pcs = nonmatched_1_pcs;
    output2.nonmatched_2_pcs = nonmatched_2_pcs;
    output2.masks_pcs = masks_pcs;
    output2.masks_union_pcs = masks_union_pcs;
    output2.masks_shifted_pcs = masks_shifted_pcs;
    output2.assignments_pcs = assignments_pcs;
    output2.matchings_pcs = matchings_pcs;
    save(fname_mat2,'-struct','output2');
    
else
    if ~exist(figdir,'dir')
        % load data for paper plots
        c1 = load(fname_mat1);
        normrMap_sm = c1.normrMap_sm;
        normpfMap_sm = c1.normpfMap_sm;
        sortIdx = c1.sortIdx;
        % pfLoc = c1.pfLoc;
        normspkRaster = c1.normspkRaster;
        ytick_files = c1.ytick_files;
        pcIdx = c1.pcIdx;
        masks = c1.masks;
            A{1:size(masks,2)} = [];
            for n = 1:size(masks,2)
                AA = zeros(size(masks{n},1)*size(masks{n},2),size(masks{n},3));
                for i = 1:size(masks{n},3)
                    mask = masks{n}(:,:,i);
                    AA(:,i) = mask(:);
                end
                A{n} = sparse(AA); 
            end
        masks_union = c1.masks_union;
            A_union{1:size(masks_union,2)} = [];
            for n = 1:size(masks_union,2)
                AA = zeros(size(masks_union{n},1)*size(masks_union{n},2),size(masks_union{n},3));
                for i = 1:size(masks_union{n},3)
                    mask = masks_union{n}(:,:,i);
                    AA(:,i) = mask(:);
                end
                A_union{n} = sparse(AA); 
            end
        assignments = c1.assignments;
        matchings = c1.matchings;

        
        % load data for other report plots
        c2 = load(fname_mat2);
        templates = c2.templates;
        templates_shifted = c2.templates_shifted;
        matched_ROIs = c2.matched_ROIs;
        nonmatched_1 = c2.nonmatched_1;
        nonmatched_2 = c2.nonmatched_2;
        masks_shifted = c2.masks_shifted;
            AA = zeros(size(masks_shifted{n},1)*size(masks_shifted{n},2),size(masks_shifted{n},3));
            A_shifted{1:size(masks_shifted,2)} = [];
            for n = 1:size(masks_shifted,2)
                for i = 1:size(masks_shifted{n},3)
                    mask = masks_shifted{n}(:,:,i);
                    AA(:,i) = mask(:);
                end
                A_shifted{n} = sparse(AA); 
            end
        matched_ROIs_pcs = c2.matched_ROIs_pcs;
        nonmatched_1_pcs = c2.nonmatched_1_pcs;
        nonmatched_2_pcs = c2.nonmatched_2_pcs;
        masks_pcs = c2.masks_pcs;
            A_pcs{1:size(masks_pcs,2)} = [];
            for n = 1:size(masks_pcs,2)
                AA = zeros(size(masks_pcs{n},1)*size(masks_pcs{n},2),size(masks_pcs{n},3));
                for i = 1:size(masks_pcs{n},3)
                    mask = masks_pcs{n}(:,:,i);
                    AA(:,i) = mask(:);
                end
                A_pcs{n} = sparse(AA); 
            end
        masks_union_pcs = c2.masks_union_pcs;
            A_union_pcs{1:size(masks_union_pcs,2)} = [];
            for n = 1:size(masks_union_pcs,2)
                AA = zeros(size(masks_union_pcs{n},1)*size(masks_union_pcs{n},2),size(masks_union_pcs{n},3));
                for i = 1:size(masks_union_pcs{n},3)
                    mask = masks_union_pcs{n}(:,:,i);
                    AA(:,i) = mask(:);
                end
                A_union_pcs{n} = sparse(AA); 
            end
        masks_shifted_pcs = c2.masks_shifted_pcs;
            A_shifted_pcs{1:size(masks_shifted_pcs,2)} = [];
            for n = 1:size(masks_shifted_pcs,2)
                AA = zeros(size(masks_shifted_pcs{n},1)*size(masks_shifted_pcs{n},2),size(masks_shifted_pcs{n},3));
                for i = 1:size(masks_shifted_pcs{n},3)
                    mask = masks_shifted_pcs{n}(:,:,i);
                    AA(:,i) = mask(:);
                end
                A_shifted_pcs{n} = sparse(AA); 
            end
        % assignments_pcs = c2.assignments_pcs;
        % matchingspcs = c2.matchings_pcs;
    end
end

if ~exist(figdir,'dir') || force
    if ~exist(figdir,'dir'), mkdir(figdir); end
    % comparison of registered templates
    Nexps = size(A,2);
    nCol = Nexps-1;
    fh1 = figure;
    for n = 1:Nexps-1
        subplot(2,nCol,n); 
            imshow(imfuse( templates{n}, templates{n+1}, 'falsecolor', 'Scaling', 'joint', 'ColorChannels', [1 2 0]) );
            str = sprintf('%g to %g: Before reg', n, n+1);
            title( str );
        
        subplot(2,nCol,nCol+n); 
            imshow( imfuse( templates_shifted{n}, templates{n+1}, 'falsecolor', 'Scaling', 'joint', 'ColorChannels', [1 2 0]) );
            str = sprintf('%g to %g: After reg', n, n+1);
            title( str );
    end
    fprintf('%s: saving multisession template registration summary figure\n',[mouseid '_' expname]);
    savefig( fh1, [figdir mouseid '_' expname '_regtemplates'] );
    saveas( fh1, [figdir mouseid '_' expname '_regtemplates'], 'png' );
    if figclose, close( fh1 ); end   
    
    % comparison of registered rois across sessions
    fh2 = figure;
    for n = 1:Nexps-1
        subplot(2,nCol,n);
            imagesc(templates{n+1}); colormap(gray); axis off; axis square;
            options.plot_bck_image = false;
            plot_contours( A{n+1}(:,matched_ROIs{n}(:,1)), templates{n+1}, options,0,[],[],'m'); hold on;
            plot_contours( A_shifted{n}(:,matched_ROIs{n}(:,2)), templates{n+1}, options,0,[],[],'c'); hold on;
            plot_contours( A{n+1}(:,nonmatched_1{n}), templates{n+1}, options, 0,[],[],'r'); hold on;
            plot_contours( A_shifted{n}(:,nonmatched_2{n}), templates{n+1}, options, 0,[],[],'b'); hold on;
            h = zeros(4, 1);
            h(1) = plot(NaN,NaN,'m');
            h(2) = plot(NaN,NaN,'c');
            h(3) = plot(NaN,NaN,'r');
            h(4) = plot(NaN,NaN,'b');
            legend(h, ['M#' num2str(n+1)], ['M#' num2str(n)], ['NoM#' num2str(n+1)], ['NoM#' num2str(n)]); hold off;
            str = sprintf('%g aligned to %g', n, n+1);
            title( str );
        subplot(2,nCol,nCol+n);
            imagesc(templates{n+1}); colormap(gray); axis off; axis square;
            options.plot_bck_image = false;
            plot_contours( A_union{n}, templates{n+1}, options,0,[],[],'w'); hold on;
            title( [num2str(n) '\cup' num2str(n+1)] );
    end    
    fprintf('%s: saving multisession roi registration summary figure\n',[mouseid '_' expname]);
    savefig( fh2, [figdir mouseid '_' expname '_regROIs'] );
    saveas( fh2, [figdir mouseid '_' expname '_regROIs'], 'png' );
    if figclose, close( fh2 ); end 
    
    % comparison of registered rois across sessions for PLACE CELLS only
    fh3 = figure;
    for n = 1:Nexps-1
        subplot(2,nCol,n);
            imagesc(templates{n+1}); colormap(gray); axis off; axis square;
            options.plot_bck_image = false;
            plot_contours( A_pcs{n+1}(:,matched_ROIs_pcs{n}(:,1)), templates{n+1}, options,0,[],[],'m'); hold on;
            plot_contours( A_shifted_pcs{n}(:,matched_ROIs_pcs{n}(:,2)), templates{n+1}, options,0,[],[],'c'); hold on;
            plot_contours( A_pcs{n+1}(:,nonmatched_1_pcs{n}), templates{n+1}, options, 0,[],[],'r'); hold on;
            plot_contours( A_shifted_pcs{n}(:,nonmatched_2_pcs{n}), templates{n+1}, options, 0,[],[],'b'); hold on;
            h = zeros(4, 1);
            h(1) = plot(NaN,NaN,'m');
            h(2) = plot(NaN,NaN,'c');
            h(3) = plot(NaN,NaN,'r');
            h(4) = plot(NaN,NaN,'b');
            legend(h, ['M#' num2str(n+1)], ['M#' num2str(n)], ['NoM#' num2str(n+1)], ['NoM#' num2str(n)]); hold off;
            str = sprintf('%g aligned to %g', n, n+1);
            title( str );
        subplot(2,nCol,nCol+n);
            imagesc(templates{n+1}); colormap(gray); axis off; axis square;
            options.plot_bck_image = false;
            plot_contours( A_union_pcs{n}, templates{n+1}, options,0,[],[],'w'); hold on;
            title( [num2str(n) '\cup' num2str(n+1)] );
    end    
    fprintf('%s: saving multisession roi registration summary figure\n',[mouseid '_' expname]);
    savefig( fh3, [figdir mouseid '_' expname '_regROIs_PCs'] );
    saveas( fh3, [figdir mouseid '_' expname '_regROIs_PCs'], 'png' );
    if figclose, close( fh3 ); end 
    
    % spike trial raster plot across sessions
    if ~exist([figdir 'rasterplots/'],'dir'), mkdir([figdir 'rasterplots/']); end
    Ncells = size(A_union{end},2);
    nRow = 10;
    nCol = 4;
    for ii=0:ceil(Ncells/nRow)-1 
        fh = figure; 
        ha = tight_subplot(nRow,nCol,[.01 .005],[.01 .07],[.01 .01]);
        for jj=0:nRow-1
            if (ii*nRow+jj+1) <= Ncells
                map = viridisMap;
                for n = 1:Nexps
                    axes(ha(jj*nCol+n));
                    idx = assignments(ii*nRow+jj+1, n);
                    if ~isnan(idx)
                        imagesc( normspkRaster{n}{idx} ); colormap(map);
                        yticks(ytick_files{n}); yticklabels(ytick_files{n}); 
                        if n == 1
                            ylabel(['Cell ' num2str(ii*nRow+jj+1)], 'fontsize', 12);
                        end
                        xticks([]); 
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
        end
        savefig( fh, [figdir 'rasterplots/' mouseid '_' expname '_normspkRaster_' num2str(ii+1)] );
        saveas( fh, [figdir 'rasterplots/' mouseid '_' expname '_normspkRaster_' num2str(ii+1)], 'png' );
        if figclose, close( fh ); end 
    end

    % place tuning across sessions
    Nbins = size(normrMap_sm{1},2);
    pf_sort{1:Nexps} = [];
    for n = 1:Nexps
        for j = 1:Nexps
            if n ~= j
                pf_sort{n}{j} = zeros(numel(sortIdx{n}),Nbins);
                for i = 1:size(sortIdx{n})
                    iIdx = sortIdx{n}(i);
                    kIdx = matchings{n}(iIdx);
                    jIdx = assignments(kIdx,j);
                    if ~isnan(jIdx)
                        pf_sort{n}{j}(i,:) = normrMap_sm{j}(jIdx,:);
                    end
                end
            else
               pf_sort{n}{j} = normpfMap_sm{n}(sortIdx{n},:); 
            end
        end
    end
    
    fh4 = figure;
    nRow = Nexps; nCol = Nexps; 
    for n = 0:Nexps-1
        for j = 0:Nexps-1
            if n~=j, k = 1; else, k = 25; end
            cmap = viridisMap;
            ax = subplot(nRow,nCol,n*nRow+j+1);
            imagesc( pf_sort{n+1}{j+1} );
            colormap(ax,cmap(k:end,:));
            xticks([]); yticks([]); 
            title(['Session ' num2str(j+1)], 'fontsize', 10);
            if j==0
                ylabel(['Session ' num2str(n+1) ' sorting']);
            end
        end
    end
    fprintf('%s: saving place tuning across sessions summary figure\n',[mouseid '_' expname]);
    savefig( fh4, [figdir mouseid '_' expname '_sortPFmaps'] );
    saveas( fh4, [figdir mouseid '_' expname '_sortPFmaps'], 'png' );
    if figclose, close( fh4 ); end 
    
    % no. of active cells, pcs, stability across sessions
    Nactivecells = zeros(Nexps,1);
    Npcs = zeros(Nexps,1);
    for n = 1:Nexps
        Nactivecells(n) = 100*(size(A{n},2)/size(A_union{end},2));
        Npcs(n) = 100*(numel(sortIdx{n})/size(A{n},2));
        
    end
    fh5 = figure;
    subplot(131); plot(1:Nexps,Nactivecells); xlabel('Session no.'); ylabel('Active cells (%)');
        hold on; plot(1:Nexps,Nactivecells,'b.','markersize',12); hold off;
    subplot(132); plot(1:Nexps,Npcs); xlabel('Session no.'); ylabel('No. of PCs (%)');
        hold on; plot(1:Nexps,Npcs,'b.','markersize',12); hold off;
    subplot(133); xlabel('Elapsed time'); ylabel('Stability');
        
    fprintf('%s: saving population summary on stability\n',[mouseid '_' expname]);
    savefig( fh5, [figdir mouseid '_' expname '_stability'] );
    saveas( fh5, [figdir mouseid '_' expname '_stability'], 'png' );
    if figclose, close( fh5 ); end 
end
    
t = toc;
str = sprintf('%s: Processing done in %g hrs\n', [mouseid(1,:) '_' env '_s1-' size(list_array,1)], round(t/3600,2));
cprintf(str)
end




