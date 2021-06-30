% Written by Ann Go
% This script compares rois across sessions. There is no need to register
% the rois as they have been segmented from an image file composed of all
% images in the file list (temporally concatenated).

function frun_ROIreg_multisession_tempcat( list, reffile, bl_prctile_array, useind, pfactivet_thr, activetrials_thr,...
            force, figclose, fname_append, fsave )
if nargin<10, fsave = true; end
if nargin<9, fname_append = []; end
if nargin<8, figclose = true; end
if nargin<7, force = false; end
if nargin<6, activetrials_thr = 0.5; end
if nargin<5, pfactivet_thr = 0.05; end

% Load module folders and define data directory
[data_locn,~,err] = load_neuroSEEmodules(false);
if ~isempty(err)
    beep
    cprintf('Errors',err);    
    return
end

% Basic settings
mcorr_method = 'normcorre';            
segment_method = 'CaImAn';
dofissa = true;
    if dofissa, str_fissa = 'FISSA'; else, str_fissa = 'noFISSA'; end
params = neuroSEE_setparams(...
            'mcorr_method',mcorr_method,...
            'segment_method',segment_method,... 
            'dofissa',dofissa,...
            'pfactivet_thr',pfactivet_thr,...
            'activetrials_thr',activetrials_thr);

tic
% mouse id and experiment name
[ mouseid, expname, fov ] = find_mouseIDexpname( list );
listfile = [data_locn 'Digital_Logbook/lists_imaging/' list];
files = extractFilenamesFromTxtfile( listfile );
if nargin<4, useind = 1:size(files,1); end
if nargin<2, reffile = files(1,:); end

% Location of processed group data for list
if ~isempty(fov)
    sdir = [data_locn 'Analysis/' mouseid '/' fov '/' mouseid '_' expname ...
       '/group_proc/imreg_' mcorr_method '_' segment_method '/' ...
       mouseid '_' expname '_imreg_ref' reffile '/'];
else
    sdir = [data_locn 'Analysis/' mouseid '/' mouseid '_' expname ...
       '/group_proc/imreg_' mcorr_method '_' segment_method '/' ...
       mouseid '_' expname '_imreg_ref' reffile '/'];
end
if ~isempty(fname_append)
    fname_mat = [sdir str_fissa '/multisessionROIs_' fname_append '/' mouseid '_' expname '_ref' ...
            reffile '_multisessionROIreg_tempcat_output.mat'];
    figdir = [sdir str_fissa '/multisessionROIs_' fname_append '/multisession_tempcat_PFdata/'];
else
    fname_mat = [sdir str_fissa '/multisessionROIs/' mouseid '_' expname '_ref' ...
                reffile '_multisessionROIreg_tempcat_output.mat'];
    figdir = [sdir str_fissa '/multisessionROIs/multisession_tempcat_PFdata/'];
end

if force || ~exist(fname_mat,'file')
    %% load segmentation output, spike and position data
    M = load([sdir mouseid '_' expname '_ref' reffile '_segment_output.mat'],'masks','tsG');
    masks = M.masks;
    
    if dofissa
        M = load([sdir str_fissa '/' mouseid '_' expname '_ref' reffile '_fissa_output.mat'],'dtsG');
        MdtsG = M.dtsG;
    else
        MtsG = M.tsG;
    end

    if ~isempty(fov)
        posdata = load([data_locn 'Analysis/' mouseid '/' fov '/' mouseid '_' expname ...
            '/group_proc/' mouseid '_' expname '_downTrackdata.mat']);
    else
        posdata = load([data_locn 'Analysis/' mouseid '/' mouseid '_' expname ...
            '/group_proc/' mouseid '_' expname '_downTrackdata.mat']);
    end
    if any(posdata.r < 100)
        params.mode_dim = '2D';                     % open field
        params.PFmap.Nbins = params.PFmap.Nbins_2D; % number of location bins in [x y]               
    else 
        params.mode_dim = '1D';                     % circular linear track
        params.PFmap.Nbins = params.PFmap.Nbins_1D; % number of location bins  
    end


    %% divide spike and position data into sessions
    % specify which files to use: all or selected?
    % if all
%     useind = 1:size(files,1);
    
    % if selected
%     useind = [1:10,15:16,17,18];
    
    % find demarcations of sessions
    t = [];
    for i = 1:numel(useind)
        t = [t (useind(i)-1)*7420+1:useind(i)*7420]; 
    end
    SdownTrackdata = structfun(@(x) x(t), posdata, 'Un', 0);
    usefiles = files(useind,:);
    Nfiles = size(usefiles,1);
    day = datetime(usefiles(:,1:8),'InputFormat','yyyyMMdd');
    ddiff = split(caldiff(day,'days'),'days');
    sdays = find(ddiff'); sdays = [1 sdays+1 Nfiles+1];
    n_sessions = numel(sdays)-1;  % number of sessions
    daylabels = cumsum([1; ddiff(find(ddiff))]);
    Vthr = params.PFmap.Vthr;
    
    for s = 1:n_sessions
        M = load([sdir str_fissa '/bl_prctile' num2str(bl_prctile_array(s)) '/' ...
        mouseid '_' expname '_ref' reffile '_spikes.mat'],'spikes');
        Sspikes{s} = M.spikes(:,t);
        if dofissa
            StsG{s} = MdtsG(:,t);
        else
            StsG{s} = MtsG(:,t);
        end
            

        % divide spike and position data into sessions
        start = (sdays(s)-1)*7420+1; last = (sdays(s+1)-1)*7420;
        spikes{s} = Sspikes{s}(:,start:last);
        tsG{s} = StsG{s}(:,start:last);
        downTrackdata{s} = structfun(@(x) x(start:last), SdownTrackdata, 'Un', 0);
    
        % for each session
        %   1. find active times (when animal was active)
        %   2. bin active spike and position data
        %   3. find active cells
        
        activephi{s}   = downTrackdata{s}.phi(downTrackdata{s}.speed > Vthr);
        activespk{s}   = spikes{s}(:,downTrackdata{s}.speed > Vthr);
        activetsG{s}   = tsG{s}(:,downTrackdata{s}.speed > Vthr);
        activet{s}     = downTrackdata{s}.time(downTrackdata{s}.speed > Vthr);
        
        [bin_phi{s},~] = discretize(activephi{s},params.PFmap.Nbins);
        activecells{s} = []; inactivecells{s} = [];
        for c = 1:size(masks,3)
            sn = GetSn(activetsG{s}(c,:));
            if max(activetsG{s}(c,:))>10*sn
                activecells{s} = [activecells{s}; c ];
            else
                inactivecells{s} = [inactivecells{s}; c ];
            end
        end
    
        % generate place field maps
        [ hist{s}, ~, PFdata{s}, ~, ~, ~ ] = generatePFmap_1d( spikes{s}, downTrackdata{s}, params );
        normspkRaster{s} = PFdata{s}.normspkRaster;
        normrMap_sm{s} = hist{s}.normrateMap_sm;
        pcIdx{s} = hist{s}.SIspk.pcIdx;
        if isfield(hist{s}.SIspk,'sortpcIdx')
            sortpcIdx{s} = hist{s}.SIspk.sortpcIdx;
        else
            sortpcIdx{s} = [];
        end
    end
    
    % make plots
    % report plots
    %   1. registered templates per session

    %   2. raster plots across sessions
    Ncells = size(masks,3);
    nRow = 10;
    nCol = n_sessions;
    cmap0 = [0.9 0.9 0.9];
    cmap1 = [0 0 1];
    cmap = zeros(50,3);
    for j=1:3
        cmap(:,j) = linspace(cmap0(j),cmap1(j),50);
    end
    colormap(cmap);
    
    for ii=0:ceil(Ncells/nRow)-1 
        fh3 = figure('Position',[680 316 535 782]); 
        ha = tight_subplot(nRow,nCol,[.05 0.04],[.03 .05],[.06 .02]);
        for jj=0:nRow-1
            c = ii*nRow+jj+1;
            if c <= Ncells
                for s = 1:n_sessions
                    axes(ha(jj*nCol+s));
                    if s == 1
                        ylabel(['Cell ' num2str(c)]);
                    end
                    if ~isnan(c)
                        imagesc( normspkRaster{s}{c} ); colormap(cmap);
                        % yticks([1 ytick_files{n}(end)]); %yticklabels(ytick_files{n}); 
                        xticks([]); 
                        if s == 1
                            ylabel(['Cell ' num2str(c)]);
                        end
                        if numel(find(pcIdx{s} == c)) > 0
                            title_str = 'PC';
                        else
                            title_str = 'nonPC';
                        end
                        title( title_str, 'fontsize', 10);
                    end
                end
            end
        end
        if fsave
            if ii == 0
                fprintf('%s: saving spike trial raster plots\n',[mouseid '_' expname]);
                if ~exist([figdir 'rasterplots/'],'dir')
                    mkdir([figdir 'rasterplots/']); 
                    fileattrib([figdir 'rasterplots/'],'+w','g','s');
                end
            end
            savefig( fh3, [figdir 'rasterplots/' mouseid '_' expname '_normspkRaster_' num2str(ii+1)] );
            saveas( fh3, [figdir 'rasterplots/' mouseid '_' expname '_normspkRaster_' num2str(ii+1)], 'png' );
        end
        if figclose, close( fh3 ); end 
    end
    
    %   3. place tuning map across sessions (n_sessions x n_sessions)
    Nbins = size(normrMap_sm{1},2);
    for n = 1:n_sessions
        for j = 1:n_sessions
            if n ~= j
                pf_sort{n}{j} = zeros(numel(sortpcIdx{n}),Nbins);
                for i = 1:size(sortpcIdx{n},1)
                    iIdx = sortpcIdx{n}(i);
                    pf_sort{n}{j}(i,:) = normrMap_sm{j}(iIdx,:);
                end
            else
               pf_sort{n}{j} = normrMap_sm{n}(sortpcIdx{n},:); 
            end
        end
    end
    
    fh3 = figure('Position', [680 547 446 551]);
    nRow = n_sessions; nCol = n_sessions; 
    for n = 0:n_sessions-1
        for j = 0:n_sessions-1
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
    if fsave
        fprintf('%s: saving place tuning across sessions summary figure\n',[mouseid '_' expname]);
        savefig( fh3, [figdir mouseid '_' expname '_sortPFmaps'] );
        saveas( fh3, [figdir mouseid '_' expname '_sortPFmaps'], 'png' );
    end
    if figclose, close( fh3 ); end 
    
    %   4. recurrence, field correlation, place field shift
    % find all session intervals
    dd = [];
    for s1 = n_sessions:-1:2
        for s2 = s1-1:-1:1
            dd = [dd; daylabels(s1)-daylabels(s2)];
        end
    end
    dintervals = sort(unique(dd));
    
    for i = 1:numel(dintervals)
        RR_ac{i} = [];
        RR_pc{i} = [];
        FC{i} = [];
        PFS{i} = [];
    end
    for s1 = 1:n_sessions-1
        for s2 = s1+1:n_sessions
            % recurrence
            % find all PCs in s1 that are also PCs in s2
            dint_ind = dintervals == abs(daylabels(s2)-daylabels(s1));
            ACs12 = activecells{s1}(ismember(activecells{s1},activecells{s2}));
            RR_ac{dint_ind} = [RR_ac{dint_ind} mean([numel(ACs12)/numel(activecells{s1}), numel(ACs12)/numel(activecells{s2})])];
            PCs12 = pcIdx{s1}(ismember(pcIdx{s1},pcIdx{s2}));
            RR_pc{dint_ind} = [RR_pc{dint_ind} mean([numel(PCs12)/numel(pcIdx{s1}), numel(PCs12)/numel(pcIdx{s2})])];

            % field correlation
            if length(PCs12) > 0
                for i = 1:length(PCs12)
                    fieldcorr(i) = corr(normrMap_sm{s1}(PCs12(i),:)', normrMap_sm{s2}(PCs12(i),:)');
                end
                FC{dint_ind} = [FC{dint_ind}; fieldcorr'];

                % place field shift
                for i = 1:length(PCs12)
                    shift(i) = (hist{s1}.pfLoc(PCs12(i)) - hist{s2}.pfLoc(PCs12(i)));
                end
                PFS{dint_ind} = [PFS{dint_ind}; shift'];
            end
        end
    end
    for i = 1:numel(dintervals)
        R_ac(i) = mean(RR_ac{i});
        R_pc(i) = mean(RR_pc{i});
        FC_mean(i) = mean(FC{i});
        PFS_mean(i) = mean(PFS{i});
    end
    
    % find cells that were active in all sessions
    for c = 1:size(masks,3)
        alwaysactivestat(c) = 1;
        for s = 1:n_sessions
            if ismember(c,activecells{s})
                alwaysactivestat(c) = 0;
            end
        end
    end
    alwaysactive = numel(find(alwaysactivestat)); 

    % find cells that were active in all sessions
    for c = 1:size(masks,3)
        alwaysplaceystat(c) = 1;
        for s = 1:n_sessions
            if ~ismember(c,pcIdx{s})
                alwaysplaceystat(c) = 0;
            end
        end
    end
    alwaysplacey = numel(find(alwaysplaceystat));
    
    fh4 = figure;
    subplot(221); 
        for s = 1:n_sessions
            Nactivecells(s) = length(activecells{s});
            frac_pcs(s) = length(pcIdx{s})/Nactivecells(s);
        end
        frac_activecells = Nactivecells./size(masks,3);
        plot(1:n_sessions,frac_activecells,'k.-',1:n_sessions,frac_pcs,'r.-','markersize',12); axis([0.95,n_sessions,0,1]);
        ylabel('Fraction of cells'); legend('Active','Place'); legend boxoff;
        box off; 
        xticks(1:n_sessions); xticklabels(daylabels);
        xlabel('Day');

    subplot(222);
        plot(1:numel(dintervals),R_ac,'k.-',1:numel(dintervals),R_pc,'r.-','markersize',12); axis([0.95,numel(dintervals),0,1]);
        ylabel('Recurrence'); % legend('Active','Place'); legend boxoff;
        box off; 
        xticks(1:numel(dintervals)); xticklabels(num2str(dintervals));
        xlabel('Days apart')

    subplot(223);
        edges = -50:5:50;
        for i = numel(dintervals):-1:1
            PFSbin(i,:) = histcounts(PFS{i},edges,'Normalization','probability');
            plot(edges(1:end-1),PFSbin(i,:),'.-','Markersize',12); hold on;
        end
        hold off; legend(num2str(dintervals)); legend boxoff;
        xlabel('Place field shift (cm)'); ylabel('Fraction of place cells');
        box off;
        axis([-55 53 0 1]);

    subplot(224);
        g = []; mFC = [];
        for i = 1:numel(dintervals)
            g = [g; repmat(num2str(dintervals(i)),length(FC{i}),1)];
            mFC = [mFC; FC{i}];
        end
        boxplot(mFC,g);
        xlabel('Days apart'); ylabel('PF correlation');
        box off;
        
    if fsave
        fprintf('%s: saving recurrence, field correlation & place field shift plots\n',[mouseid '_' expname]);
        savefig( fh4, [figdir 'recurrennce_fieldcorr_shift'] );
        saveas( fh4, [figdir 'recurrennce_fieldcorr_shift'], 'png' );
    end
    if figclose, close( fh4 ); end 

    
    %   5. subset of ROIs active in each session
%     for j = 1:size(masks,3)
%         outline{:,:,j} = bwboundaries(masks(:,:,j));    % boundary of each ROI
%     end
% 
%     M = load([sdir mouseid '_' expname '_ref' reffile '_imreg_template.mat'],'template_g');
%     template = M.template_g;
%     
%     fh5 = figure; 
%     for s = 1:n_sessions
%         subplot(2,n_sessions,s); 
%         imagesc(template); colormap(gray); axis square; axis off; hold on
%         for j = 1:length(activecells{s})
%             plot(outline{1,1,j}{1}(:,2),outline{1,1,j}{1}(:,1),'w','Linewidth',1);
%         end
%         hold off;
%         title(['Active in session' num2str(s)]);
%         
%         subplot(2,n_sessions,n_sessions+s); 
%         imagesc(template); colormap(gray); axis square; axis off; hold on
%         for j = 1:length(activecells{s})
%             plot(outline{1,1,j}{1}(:,2),outline{1,1,j}{1}(:,1),'y','Linewidth',1);
%         end
%         hold off;
%         title(['Inactive in session' num2str(s)]);
%     end
%     if fsave
%         if ~exist(figdir,'dir'), mkdir(figdir); end
%         savefig( fh5, [figdir 'activecells'] );
%         saveas( fh5, [figdir 'activecells'], 'png' );
%     end
%     if figclose, close(fh2); end

    % save output
    if fsave
        if ~exist([sdir str_fissa '/multisessionROIs_' fname_append '/'],'dir') 
            mkdir([sdir str_fissa '/multisessionROIs_' fname_append '/']);
            fileattrib([sdir str_fissa '/multisessionROIs_' fname_append '/'],'+w','g','s');
        end
        save(fname_mat,'params','daylabels','usefiles','normspkRaster','pcIdx','sortpcIdx','normrMap_sm',...
                       'R_ac','R_pc','FC','FC_mean','PFS','PFS_mean','alwaysactive','alwaysplacey',...
                       'n_sessions','frac_activecells','frac_pcs','dintervals','PFSbin');
    end
    

else
%     if ~exist(figdir,'dir')
%         % make plots
%     end
end


t = toc;
str = sprintf('%s: Processing done in %g hrs\n', [mouseid '_' expname], round(t/3600,2));
cprintf(str)

end