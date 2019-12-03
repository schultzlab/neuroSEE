% Written by Ann Go
% list   :  name of text file containing filenames of files to be compared.
%           Typically in the format 'list_m##_expname.txt' or 'list_m##.txt'
%           (when wanting to include all files for mouse).
%
% Comparison figures will be saved in
%           /...thefarm2/CrazyEights/AD_2PCa/Analysis/m##/m##_expname/
%   e.g.    'list_m62_fam1nov-fam1.txt'         - all fam1 files in fam1nov experiment
%           'list_m62_fam1nov_fam1fam1rev.txt'  - all files in fam1nov and
%                                                   fam1fam1rev experiments
%           'list_m79_fam1_s1-5.txt'            - all fam1 files across 5 sessions           
%           'list_m86_open_s1-2.txt'            - all open field files across 2 sessions

function frun_collate_indivproc_results_batch( list, force )

if nargin<2, force = false; end
    
%% Load module folders and define data directory
test = false;   % flag to use test file directory

[data_locn,~,err] = load_neuroSEEmodules(test);
if ~isempty(err)
    beep
    cprintf('Errors',err);    
    return
end

mcorr_method = 'normcorre-nr';
segment_method = 'CaImAn';
dofissa = true;
    if dofissa, str_fissa = 'FISSA'; else, str_fissa = 'noFISSA'; end

%% MouseID and experiment name
[ mouseid, expname ] = find_mouseIDexpname(list);

%% Files
listfile = [data_locn 'Digital_Logbook/lists/' list];
files = extractFilenamesFromTxtfile( listfile );
Nfiles = size(files,1);

[nRow, nCol] = getnRownCol(Nfiles);
nPlot = nCol*nRow;

%% Load image data for each recording
sdir1 = [data_locn 'Analysis/' mouseid '/' mouseid '_' expname '/' mcorr_method '_' segment_method '_' str_fissa '/indiv_mcorr/'];
    if ~exist(sdir1,'dir'), mkdir(sdir1); end
    
if any([ force, ~exist([sdir1 mouseid '_' expname '_GREEN_mcorr.fig'],'file'),...
                ~exist([sdir1 mouseid '_' expname '_RED_mcorr.fig'],'file') ])
                
    for i = 1:Nfiles
        file = files(i,:);
        fname = [data_locn 'Data/' file(1:8) '/Processed/' file '/mcorr_' mcorr_method '/' file '_mcorr_output.mat'];
        if exist(fname,'file')
            c = load(fname);
            M(i).green = c.green;
            M(i).red = c.red;
        end
    end
    clear c

    %% Compare green and red images to see which ones are most similar

    for ii=0:(Nfiles/nPlot)-1
        fh = figure;
        ha = tight_subplot(nRow,nCol,[.04 .01],[.01 .07],[.01 .01]);
        for jj=0:nPlot-1
            if (ii*nPlot+jj+1) <= Nfiles
                if ~isempty(M(ii*nPlot+jj+1).green)
                    axes(ha(ii*nPlot+jj+1));
                    imagesc(M(ii*nPlot+jj+1).green.meanregframe./max(max(M(ii*nPlot+jj+1).green.meanregframe))); 
                    axis off; colormap(gray);
                    str = files(ii*nPlot+jj+1,:);
                    title([str(5:6) '-' str(7:8) ' ' str(10:11) ':' str(13:14)]);
                end
            end
        end
        axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0  1],'Box','off',...
            'Visible','off','Units','normalized', 'clipping' , 'off');
            titletext = [mouseid '-' expname ': GREEN channel'];
            ind = strfind(titletext,'_');
            titletext(ind) = '-';
            text('Position',[0.25 0.99], 'FontSize',14, 'String',titletext);
    
        fname_fig = [sdir1 mouseid '_' expname '_GREEN_mcorr.fig'];
        savefig( fh, fname_fig );
        saveas( fh, fname_fig(1:end-4), 'png' );
        %close( fh );
    end

    for ii=0:(Nfiles/nPlot)-1
        fh = figure;
        ha = tight_subplot(nRow,nCol,[.04 .01],[.01 .07],[.01 .01]);
        for jj=0:nPlot-1
            if (ii*nPlot+jj+1) <= Nfiles
                if ~isempty(M(ii*nPlot+jj+1).red)
                    axes(ha(ii*nPlot+jj+1));
                    imagesc(M(ii*nPlot+jj+1).red.meanregframe./max(max(M(ii*nPlot+jj+1).red.meanregframe))); 
                    axis off; colormap(gray);
                    str = files(ii*nPlot+jj+1,:);
                    title([str(5:6) '-' str(7:8) ' ' str(10:11) ':' str(13:14)]);
                end
            end
        end
        axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0  1],'Box','off',...
            'Visible','off','Units','normalized', 'clipping' , 'off');
            titletext = [mouseid '-' expname ': RED channel'];
            ind = strfind(titletext,'_');
            titletext(ind) = '-';
            text('Position',[0.25 0.99], 'FontSize',14, 'String',titletext);
        fname_fig = [sdir1 mouseid '_' expname '_RED_mcorr.fig'];
        savefig( fh, fname_fig );
        saveas( fh, fname_fig(1:end-4), 'png' );
        close( fh );
    end 
end


%% Load tracking data for each recording
sdir2 = [data_locn 'Analysis/' mouseid '/' mouseid '_' expname '/' mcorr_method '_' segment_method '_' str_fissa '/indiv_trajectories/'];
    if ~exist(sdir2,'dir'), mkdir(sdir2); end

if force || ~exist([sdir2 mouseid '_' expname '_traj.fig'],'file')
    for i = 1:Nfiles
        trackfile = findMatchingTrackingFile(data_locn, files(i,:), 0);
        c = load_trackfile(data_locn, files(i,:), trackfile, 0);
        M(i).trackdata.x = c.x;
        M(i).trackdata.y = c.y;
        M(i).trackdata.r = c.r;
        M(i).trackdata.phi = c.phi;
        M(i).trackdata.speed = c.speed;
        M(i).trackdata.time = c.time;
    end
    clear c


    %% Compare trajectories
    for ii=0:(Nfiles/nPlot)-1
        fh = figure;
        ha = tight_subplot(nRow,nCol,[.04 .01],[.01 .07],[.01 .01]);
        for jj=0:nPlot-1
            if (ii*nPlot+jj+1) <= Nfiles
                if ~isempty(M(ii*nPlot+jj+1).trackdata)
                    axes(ha(ii*nPlot+jj+1));
                    plot(M(ii*nPlot+jj+1).trackdata.x,M(ii*nPlot+jj+1).trackdata.y); 
                    axis off; 
                    str = files(ii*nPlot+jj+1,:);
                    title([str(5:6) '-' str(7:8) ' ' str(10:11) ':' str(13:14)]);
                end
            end
        end
        axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0  1],'Box','off',...
            'Visible','off','Units','normalized', 'clipping' , 'off');
            titletext = [mouseid '-' expname ': trajectories'];
            ind = strfind(titletext,'_');
            titletext(ind) = '-';
            text('Position',[0.4 0.99], 'FontSize',14, 'String',titletext);
        fname_fig = [sdir2 mouseid '_' expname '_traj.fig'];
        savefig( fh, fname_fig );
        saveas( fh, fname_fig(1:end-4), 'png' );
        close( fh );
    end 
else
    trackfile = findMatchingTrackingFile(data_locn, files(1,:), 0);
    c = load_trackfile(data_locn, files(1,:), trackfile, 0);
    M(1).trackdata.r = c.r;
end


%% Load PF maps
% Determine mode_dim
if any(M(1).trackdata.r < 100)
    mode_dim = '2D'; % open field
else 
    mode_dim = '1D'; % circular linear track
end

sdir3 = [data_locn 'Analysis/' mouseid '/' mouseid '_' expname '/' mcorr_method '_' segment_method '_' str_fissa '/indiv_PFmaps/'];
    if ~exist(sdir3,'dir'), mkdir(sdir3); end

if (force || ~exist([sdir3 mouseid '_' expname '_PFmaps.fig'],'file')) && strcmpi(mode_dim,'1D')
    for i = 1:Nfiles
        file = files(i,:);
        PFmapfile = [data_locn 'Data/' file(1:8) '/Processed/' file '/mcorr_' mcorr_method '/' segment_method '/'...
                     str_fissa '/PFmaps/' file '_PFmap_output.mat'];
        if exist(PFmapfile,'file')
            c = load(PFmapfile);
            M(i).PFmap = c.hist.sort_normpfMap_sm(:,:,1);
        end
    end
    clear c

    %% Compare PF maps
    for ii=0:(Nfiles/nPlot)-1
        fh = figure;
        ha = tight_subplot(nRow,nCol,[.04 .01],[.01 .07],[.01 .01]);
        for jj=0:nPlot-1
            if (ii*nPlot+jj+1) <= Nfiles
                if ~isempty(M(ii*nPlot+jj+1).PFmap)
                    axes(ha(ii*nPlot+jj+1));
                    imagesc(M(ii*nPlot+jj+1).PFmap); 
                    axis off; 
                    str = files(ii*nPlot+jj+1,:);
                    title([str(5:6) '-' str(7:8) ' ' str(10:11) ':' str(13:14)]);
                end
            end
        end
        axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0  1],'Box','off',...
            'Visible','off','Units','normalized', 'clipping' , 'off');
            titletext = [mouseid '-' expname ': PF maps'];
            ind = strfind(titletext,'_');
            titletext(ind) = '-';
            text('Position',[0.25 0.99], 'FontSize',14, 'String',titletext);
        fname_fig = [sdir3 mouseid '_' expname '_PFmaps.fig'];
        savefig( fh, fname_fig );
        saveas( fh, fname_fig(1:end-4), 'png' );
        close( fh );
    end 
end


%% Load ROIs
sdir4 = [data_locn 'Analysis/' mouseid '/' mouseid '_' expname '/' mcorr_method '_' segment_method '_' str_fissa '/indiv_ROIs/'];
    if ~exist(sdir4,'dir'), mkdir(sdir4); end

if (force || ~exist([sdir4 mouseid '_' expname '_ROIs.fig'],'file')) 
    for i = 1:Nfiles
        file = files(i,:);
        ROIfile = [data_locn 'Data/' file(1:8) '/Processed/' file '/mcorr_' mcorr_method '/' segment_method '/'...
                   file '_segment_output.mat'];
        if exist(ROIfile,'file')
            c = load(ROIfile);
            M(i).corr_image = c.corr_image;
            M(i).masks = c.masks;
            for j = 1:size(c.masks,3)
                M(i).outline{:,:,j} = bwboundaries(c.masks(:,:,j));    % boundary of each ROI
            end
        end
    end
    clear c

    %% Compare ROIs
    plotopts.plot_ids = 0; % set to 1 to view the ID number of the ROIs on the plot
    for ii=0:(Nfiles/nPlot)-1
        fh = figure;
        ha = tight_subplot(nRow,nCol,[.04 .01],[.01 .07],[.01 .01]);
        for jj=0:nPlot-1
            if (ii*nPlot+jj+1) <= Nfiles
                if ~isempty(M(ii*nPlot+jj+1).outline)
                    axes(ha(ii*nPlot+jj+1));
                    imagesc(M(ii*nPlot+jj+1).corr_image); colormap(gray); hold on;
                    for j = 1:size(M(ii*nPlot+jj+1).masks,3)
                        plot(M(ii*nPlot+jj+1).outline{1,1,j}{1}(:,2),M(ii*nPlot+jj+1).outline{1,1,j}{1}(:,1),'w','Linewidth',1);
                    end
                    hold off; axis off; 
                    str = files(ii*nPlot+jj+1,:);
                    title([str(5:6) '-' str(7:8) ' ' str(10:11) ':' str(13:14)]);
                end
            end
        end
        axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0  1],'Box','off',...
            'Visible','off','Units','normalized', 'clipping' , 'off');
            titletext = [mouseid '-' expname ': ROIs'];
            ind = strfind(titletext,'_');
            titletext(ind) = '-';
            text('Position',[0.25 0.99], 'FontSize',14, 'String',titletext);
        fname_fig = [sdir4 mouseid '_' expname '_ROIs.fig'];
        savefig( fh, fname_fig );
        saveas( fh, fname_fig(1:end-4), 'png' );
        close( fh );
    end 
end



