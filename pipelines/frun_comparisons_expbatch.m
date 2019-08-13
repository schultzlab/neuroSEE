% Written by Ann Go
% expname:  Experiment name in the format 'm##_experiments'. This
%           determines the name of the folder to which comparison figures will be
%           saved in
%           /...thefarm2/CrazyEights/AD_2PCa/Analysis/m##/individual file comparisons/
%   e.g.    'm62_fam1nov_fam1fam1rev' - compares all files in fam1nov and
%                                           fam1fam1rev experiments
%           'm79_fam1_s1-5'
%           'm86_open_s1-2'
%
% list   :  File name of text file containing filenames of images to be compared.
%           Typically in the format 'list_expname.txt' but may be shortened to 'list_m##.txt'
%           when wanting to include all files for mouse.

function frun_comparisons_expbatch( expname, list, force )

if nargin<3, force = false; end
    
%% Load module folders and define data directory
test = false;                   % flag to use one of smaller files in test folder)
[data_locn,~,err] = load_neuroSEEmodules(test);
if ~isempty(err)
    beep
    cprintf('Errors',err);    
    return
end

%% Files
listfile = [data_locn 'Digital_Logbook/lists/' list];
files = extractFilenamesFromTxtfile( listfile );
Nfiles = size(files,1);
mouseid = expname(1:3);

[nRow, nCol] = getnRownCol(Nfiles);
nPlot = nCol*nRow;


%% Load image data for each recording
sdir1 = [data_locn 'Analysis/' mouseid '/individual_file_comparisons/' expname '_mcorr/'];
    if ~exist(sdir1,'dir'), mkdir(sdir1); end
    
if any([ force, ~exist([sdir1 expname '_GREEN_mcorr.fig'],'file'),...
                ~exist([sdir1 expname '_RED_mcorr.fig'],'file') ])
                
    for i = 1:Nfiles
        file = files(i,:);
        fname = [data_locn 'Data/' file(1:8) '/Processed/' file '/mcorr_normcorre/' file '_mcorr_output.mat'];
        if exist(fname,'file')
            c = load(fname);
            M(i).green = c.green;
            M(i).red = c.red;
        end
    end
    clear c

    %% Compare green and red images to see which ones are most similar

    for ii=0:Nfiles/nPlot
        fh = figure;
        ha = tight_subplot(nRow,nCol,[.01 .01],[.01 .05],[.01 .01]);
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
            titletext = [expname ': GREEN channel'];
            ind = strfind(titletext,'_');
            titletext(ind) = '-';
            text('Position',[0.35 0.99], 'FontSize',14, 'String',titletext);
    end 
    fname_fig = [sdir1 expname '_GREEN_mcorr.fig'];
        savefig( fh, fname_fig );
        saveas( fh, fname_fig(1:end-4), 'png' );
        close( fh );


    for ii=0:Nfiles/nPlot
        fh = figure;
        ha = tight_subplot(nRow,nCol,[.01 .01],[.01 .05],[.01 .01]);
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
            titletext = [expname ': RED channel'];
            ind = strfind(titletext,'_');
            titletext(ind) = '-';
            text('Position',[0.35 0.99], 'FontSize',14, 'String',titletext);
    end 
    fname_fig = [sdir1 expname '_RED_mcorr.fig'];
        savefig( fh, fname_fig );
        saveas( fh, fname_fig(1:end-4), 'png' );
        close( fh );
end


%% Load tracking data for each recording
sdir2 = [data_locn 'Analysis/' mouseid '/individual_file_comparisons/' expname '_trajectories/'];
    if ~exist(sdir2,'dir'), mkdir(sdir2); end

if force || ~exist([sdir2 expname '_traj.fig'],'file')
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
    for ii=0:Nfiles/nPlot
        fh = figure;
        ha = tight_subplot(nRow,nCol,[.01 .01],[.01 .05],[.01 .01]);
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
            titletext = [expname ': trajectories'];
            ind = strfind(titletext,'_');
            titletext(ind) = '-';
            text('Position',[0.4 0.99], 'FontSize',14, 'String',titletext);
    end 
    fname_fig = [sdir2 expname '_traj.fig'];
        savefig( fh, fname_fig );
        saveas( fh, fname_fig(1:end-4), 'png' );
        close( fh );
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


sdir3 = [data_locn 'Analysis/' mouseid '/individual_file_comparisons/' expname '_PFmaps/'];
    if ~exist(sdir3,'dir'), mkdir(sdir3); end

if (force || ~exist([sdir3 expname '_PFmaps.fig'],'file')) && strcmpi(mode_dim,'1D')
    for i = 1:Nfiles
        file = files(i,:);
        PFmapfile = [data_locn 'Data/' file(1:8) '/Processed/' file '/mcorr_normcorre/CaImAn/FISSA/PFmaps/' file '_PFmap_output.mat'];
        if exist(PFmapfile,'file')
            c = load(PFmapfile);
            M(i).PFmap = c.hist.sort_normpfMap_sm(:,:,1);
        end
    end
    clear c


    %% Compare PF maps
    for ii=0:Nfiles/nPlot
        fh = figure;
        ha = tight_subplot(nRow,nCol,[.01 .01],[.01 .05],[.01 .01]);
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
            titletext = [expname ': PF maps'];
            ind = strfind(titletext,'_');
            titletext(ind) = '-';
            text('Position',[0.35 0.99], 'FontSize',14, 'String',titletext);
    end 
    fname_fig = [sdir3 expname '_PFmaps.fig'];
        savefig( fh, fname_fig );
        saveas( fh, fname_fig(1:end-4), 'png' );
        close( fh );
end

%% Load ROIs
sdir4 = [data_locn 'Analysis/' mouseid '/individual_file_comparisons/' expname '_ROIs/'];
    if ~exist(sdir4,'dir'), mkdir(sdir4); end

if (force || ~exist([sdir4 expname '_ROIs.fig'],'file')) 
    for i = 1:Nfiles
        file = files(i,:);
        ROIfile = [data_locn 'Data/' file(1:8) '/Processed/' file '/mcorr_normcorre/CaImAn/' file '_segment_output.mat'];
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
    for ii=0:Nfiles/nPlot
        fh = figure;
        ha = tight_subplot(nRow,nCol,[.01 .01],[.01 .05],[.01 .01]);
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
            titletext = [expname ': ROIs'];
            ind = strfind(titletext,'_');
            titletext(ind) = '-';
            text('Position',[0.35 0.99], 'FontSize',14, 'String',titletext);
    end 
    fname_fig = [sdir4 expname '_ROIs.fig'];
        savefig( fh, fname_fig );
        saveas( fh, fname_fig(1:end-4), 'png' );
        close( fh );
end



