% Written by Ann Go
% This script compares rois across sessions. There is no need to register
% the rois as they have been segmented from on image file composed of all
% images in the file list (temporally concatenated).

function frun_ROIreg_multisession_tempcat( list, reffile, bl_prctile, force, figclose )
if nargin<5, figclose = true; end
if nargin<4, force = false; end
if nargin<3, bl_prctile = 85; end

% Basic settings
mcorr_method = 'normcorre';            
segment_method = 'CaImAn';
dofissa = true;
    if dofissa, str_fissa = 'FISSA'; else, str_fissa = 'noFISSA'; end
params = neuroSEE_setparams('mcorr_method',mcorr_method,...
            'segment_method',segment_method, 'dofissa',dofissa);

% Load module folders and define data directory
[data_locn,~,err] = load_neuroSEEmodules(false);
if ~isempty(err)
    beep
    cprintf('Errors',err);    
    return
end

% mouse id and experiment name
[ mouseid, expname ] = find_mouseIDexpname( list );
listfile = [data_locn 'Digital_Logbook/lists_imaging/' list];
files = extractFilenamesFromTxtfile( listfile );
if nargin<2, reffile = files(1,:); end
Nfiles = size(files,1);

% Location of processed group data for list
sdir = [data_locn 'Analysis/' mouseid '/' mouseid '_' expname ...
       '/group_proc/imreg_' mcorr_method '_' segment_method '/' ...
       mouseid '_' expname '_imreg_ref' reffile '/'];
fname_mat = [sdir str_fissa '/multisessionROIreg/' mouseid '_' expname '_ref' ...
            reffile '_multisessionROIreg_tempcat_output.mat'];
figdir = [sdir str_fissa '/multisessionROIreg/multisession_tempcat_PFdata/'];

if force || ~exist(fname_mat,'file')
    %% load segmentation output, spike and position data
    M = load([sdir mouseid '_' expname '_ref' reffile '_segment_output.mat'],'masks');
    masks = M.masks;

    M = load([sdir str_fissa '/bl_prctile' num2str(bl_prctile) '/' ...
            mouseid '_' expname '_ref' reffile '_spikes.mat'],'spikes');
    Sspikes = M.spikes;
    
    SdownTrackdata = load([data_locn 'Analysis/' mouseid '/' mouseid '_' expname ...
            '/group_proc/' mouseid '_' expname '_downTrackdata.mat']);

    %% divide spike and position data into sessions
    % find demarcations of sessions
    for f = 1:Nfiles
        day(s) = str2double(files(s,7:8));
    end
    ddiff = [0 diff(day)];
    sdays = find(ddiff); sdays = [1 sdays Nfiles];
    n_sessions = numel(sdays)-1;  % number of sessions
    
    % divide spike and position data into sessions
    for s = 1:n_sessions
        start = (sdays(s)-1)*7420+1; last = sdays(s+1)*7420-1;
        spikes{s} = Sspikes(:,start:last);
        downTrackdata{s}.phi = SdownTrackdata.phi(start:last);
        downTrackdata{s}.x = SdownTrackdata.x(start:last);
        downTrackdata{s}.y = SdownTrackdata.y(start:last);
        downTrackdata{s}.speed = SdownTrackdata.speed(start:last);
        downTrackdata{s}.r = SdownTrackdata.r(start:last);
        downTrackdata{s}.time = SdownTrackdata.time(start:last);
    end
    
    % for each session
    %   1. find active times (when animal was active)
    %   2. bin active spike and position data
    %   3. find active cells
    Vthr = params.PFmap.Vthr;
    if any(SdownTrackdata.r < 100)
        params.mode_dim = '2D';                     % open field
        params.PFmap.Nbins = params.PFmap.Nbins_2D; % number of location bins in [x y]               
    else 
        params.mode_dim = '1D';                     % circular linear track
        params.PFmap.Nbins = params.PFmap.Nbins_1D; % number of location bins  
    end

    for s = 1:n_sessions
        activephi{s}   = downTrackdata{s}.phi(downTrackdata{s}.speed > Vthr);
        activespk{s}   = spikes{s}(:,downTrackdata{s}.speed > Vthr);
        activet{s}     = downTrackdata{s}.time(downTrackdata{s}.speed > Vthr);
        
        [bin_phi{s},~] = discretize(activephi{s},params.PFmap.Nbins);
        spk = activespk{s};
        activecells{s} = find(sum(spk,2));
        inactivecells{s} = find(sum(spk,2)==0);
    end
    
    %   3. generate place field maps
    for s = 1:n_sessions
        [ hist{s}, ~, PFdata{s}, ~, ~, ~ ] = generatePFmap_1d( spikes{s}, downTrackdata{s}, params );
        normspkRaster{s} = PFdata{s}.normspkRaster;
        pcIdx{s} = hist{s}.SIspk.pcIdx;
        sortpcIdx{s} = hist{s}.SIspk.sortpcIdx;
        normrMap_sm{s} = hist{s}.normrateMap_sm;
    end
    
    % make plots
    % report plots
    %   1. registered templates per session
    %   5. # active cells, # pcs, stability across sessions

    %   2. subset of ROIs active in each session
    for j = 1:size(masks,3)
        outline{:,:,j} = bwboundaries(masks(:,:,j));    % boundary of each ROI
    end

    M = load([sdir mouseid '_' expname '_ref' reffile '_imreg_template.mat'],'template_g');
    template = M.template_g;
    
    figure; 
    for s = 1:n_sessions
        subplot(2,n_sessions,s); 
        imagesc(template); colormap(gray); axis square; axis off; hold on
        for j = 1:length(activecells{s})
            plot(outline{1,1,j}{1}(:,2),outline{1,1,j}{1}(:,1),'y','Linewidth',1);
        end
        hold off;
        title(['Active in session' num2str(s)]);
        
        subplot(2,n_sessions,n_sessions+s); 
        imagesc(template); colormap(gray); axis square; axis off; hold on
        for j = 1:length(activecells{s})
            plot(outline{1,1,j}{1}(:,2),outline{1,1,j}{1}(:,1),'y','Linewidth',1);
        end
        hold off;
        title(['Inactive in session' num2str(s)]);
        
%         if ~isempty(figname_pref)
%             savefig(fig, [figname_pref '_ROIs']);
%             saveas(fig, [figname_pref '_ROIs'], 'png');
%         end
        if figclose, close(fig); end
    end
    
    % 3. raster plots across sessions
    if ~exist([figdir 'rasterplots/'],'dir'), mkdir([figdir 'rasterplots/']); end
    Ncells = size(masks,3);
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
%         if ii == 0
%             fprintf('%s: saving spike trial raster plots\n',[mouseid '_' expname]);
%         end
%         savefig( fh, [figdir 'rasterplots/' mouseid '_' expname '_normspkRaster_' num2str(ii+1)] );
%         saveas( fh, [figdir 'rasterplots/' mouseid '_' expname '_normspkRaster_' num2str(ii+1)], 'png' );
%         if figclose, close( fh ); end 
    end
    
    %   4. place tuning map across sessions (n_sessions x n_sessions)
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
    
    fh4 = figure('Position', [680 547 446 551]);
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
%     fprintf('%s: saving place tuning across sessions summary figure\n',[mouseid '_' expname]);
%     savefig( fh4, [figdir mouseid '_' expname '_sortPFmaps'] );
%     saveas( fh4, [figdir mouseid '_' expname '_sortPFmaps'], 'png' );
%     if figclose, close( fh4 ); end 
    
else
    if ~exist(figdir,'dir')
        % make plots
    end
end


t = toc;
str = sprintf('%s: Processing done in %g hrs\n', [mouseid '_' expname], round(t/3600,2));
cprintf(str)

end