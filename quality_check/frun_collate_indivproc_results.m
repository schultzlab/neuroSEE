% Written by Ann Go
%
% This script collates the results for the complete pipeline processing of
% the individual files in 'list'. Figures comparing average motion
% corrected frames, ROIs and trajectories are generated and saved
% in
%   /...thefarm2/CrazyEights/AD_2PCa/Analysis/m##/m##_expname/individual_proc/
%
% INPUTS
% list   : name of text file containing filenames of files to be compared.
%           Typically in the format 'list_m##_expname.txt'.
%   e.g.    'list_m62_fam1nov-fam1.txt'         - all fam1 files in fam1nov experiment
%           'list_m62_fam1nov.txt'              - all files in fam1nov experiments
%           'list_m79_fam1_s1-5.txt'            - all fam1 files across 5 sessions           
%           'list_m86_open_s1-2.txt'            - all open field files across 2 sessions
%
% force  : (optional, default: false) flag to force generation of comparison figures even though they
%           already exist

function frun_collate_indivproc_results( list, force )

if nargin<2, force = false; end

addpath(genpath('../pipelines'));

%% Load module folders and define data directory
[data_locn,~,err] = load_neuroSEEmodules;
if ~isempty(err)
    beep
    cprintf('Errors',err);    
    return
end

% mcorr_method = 'normcorre';
mcorr_method = 'normcorre';
segment_method = 'CaImAn';
dofissa = true;
    if dofissa, str_fissa = 'FISSA'; else, str_fissa = 'noFISSA'; end

%% MouseID and experiment name
[ mouseid, expname, fov ] = find_mouseIDexpname(list);

%% Files
listfile = [data_locn 'Digital_Logbook/lists_imaging/' list];
files = extractFilenamesFromTxtfile( listfile );
Nfiles = size(files,1);

[nRow, nCol] = getnRownCol(Nfiles);
nPlot = nCol*nRow;
Nfig = (Nfiles/nPlot)-1;
if Nfig < 0, Nfig = 0; end


%% Load image data for each file/recording
if ~isempty(fov)
    sdir1 = [data_locn 'Analysis/' mouseid '/' fov '/' mouseid '_' expname '/individual_proc/indiv_' ...
         mcorr_method '_' segment_method '_' str_fissa '/indiv_mcorr/'];
else
    sdir1 = [data_locn 'Analysis/' mouseid '/' mouseid '_' expname '/individual_proc/indiv_' ...
         mcorr_method '_' segment_method '_' str_fissa '/indiv_mcorr/'];
end
    
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

    % Compare green and red images to see which ones are most similar

    if exist('M','var') && isfield(M,'green')
        for ii=0:Nfig
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
                text('Position',[0.01 0.99], 'FontSize',14, 'String',titletext);

            if ~exist(sdir1,'dir'), mkdir(sdir1); fileattrib(sdir1,'+w','g','s'); end    
            fname_fig = [sdir1 mouseid '_' expname '_GREEN_mcorr.fig'];
            savefig( fh, fname_fig );
            saveas( fh, fname_fig(1:end-4), 'png' );
            close( fh );
        end

        for ii=0:Nfig
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
                text('Position',[0.01 0.99], 'FontSize',14, 'String',titletext);
            fname_fig = [sdir1 mouseid '_' expname '_RED_mcorr.fig'];
            savefig( fh, fname_fig );
            saveas( fh, fname_fig(1:end-4), 'png' );
            close( fh );
        end 
    else
        fprintf('%s: No motion corrected files found\n', list)
    end
end


%% Load ROIs for each recording
if ~isempty(fov)
    sdir4 = [data_locn 'Analysis/' mouseid '/' fov '/' mouseid '_' expname '/individual_proc/indiv_' ...
         mcorr_method '_' segment_method '_' str_fissa '/indiv_ROIs/'];
else
    sdir4 = [data_locn 'Analysis/' mouseid '/' mouseid '_' expname '/individual_proc/indiv_' ...
         mcorr_method '_' segment_method '_' str_fissa '/indiv_ROIs/'];
end
    
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

    if exist('M','var') && isfield(M,'masks')
        %% Compare ROIs
        plotopts.plot_ids = 0; % set to 1 to view the ID number of the ROIs on the plot
        for ii=0:Nfig
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
                text('Position',[0.01 0.99], 'FontSize',14, 'String',titletext);
            if ~exist(sdir4,'dir'), mkdir(sdir4); fileattrib(sdir4,'+w','g','s'); end
            fname_fig = [sdir4 mouseid '_' expname '_ROIs.fig'];
            savefig( fh, fname_fig );
            saveas( fh, fname_fig(1:end-4), 'png' );
            close( fh );
        end 
    else
        fprintf('%s: ROIs not found\n', list)
    end
end

%% Load tracking data for each recording
if ~isempty(fov)
    sdir2 = [data_locn 'Analysis/' mouseid '/' fov '/' mouseid '_' expname '/individual_proc/indiv_' ...
         mcorr_method '_' segment_method '_' str_fissa '/indiv_trajectories/'];
else
    sdir2 = [data_locn 'Analysis/' mouseid '/' mouseid '_' expname '/individual_proc/indiv_' ...
         mcorr_method '_' segment_method '_' str_fissa '/indiv_trajectories/'];
end
    
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

    if exist('M','var') && isfield(M,'trackdata')
        % Compare trajectories
        for ii=0:Nfig
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
                text('Position',[0.01 0.99], 'FontSize',14, 'String',titletext);
            if ~exist(sdir2,'dir'), mkdir(sdir2); fileattrib(sdir2,'+w','g','s'); end
            fname_fig = [sdir2 mouseid '_' expname '_traj.fig'];
            savefig( fh, fname_fig );
            saveas( fh, fname_fig(1:end-4), 'png' );
            close( fh );
        end 
    else
        fprintf('%s: No tracking files found\n', list)
    end
else
    trackfile = findMatchingTrackingFile(data_locn, files(1,:), 0);
    c = load_trackfile(data_locn, files(1,:), trackfile, 0);
    M(1).trackdata.r = c.r;
end





