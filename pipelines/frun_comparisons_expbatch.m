function frun_comparisons_expbatch( expname, list )

%% Load module folders and define data directory
test = false;                   % flag to use one of smaller files in test folder)
[data_locn,~,err] = load_neuroSEEmodules(test);
if ~isempty(err)
    beep
    cprintf('Errors',err);    
    return
end

%% Files
files = extractFilenamesFromTxtfile( list );
Nfiles = size(files,1);
mouseid = expname(1:3);

%% Load image data for each recording
for i = 1:Nfiles
    file = files(i,:);
    fname = [data_locn 'Data/' file(1:8) '/Processed/' file '/mcorr_normcorre/' file '_mcorr_output.mat'];
    c = load(fname);
    M(i).green = c.green;
    M(i).red = c.red;
end
clear c

%% Compare green and red images to see which ones are most similar
 
[nRow, nCol] = getnRownCol(Nfiles);
nPlot = nCol*nRow;
for ii=0:Nfiles/nPlot
    fh = figure;
    ha = tight_subplot(nRow,nCol,[.01 .01],[.01 .05],[.01 .01]);
    for jj=0:nPlot-1
        if (ii*nPlot+jj+1) <= Nfiles
            axes(ha(ii*nPlot+jj+1));
            imagesc(M(ii*nPlot+jj+1).green.meanregframe./max(max(M(ii*nPlot+jj+1).green.meanregframe))); 
            axis off; colormap(gray);
            str = files(ii*nPlot+jj+1,:);
            title([str(5:6) '-' str(7:8) ' ' str(10:11) ':' str(13:14)]);
        end
    end
    axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0  1],'Box','off',...
        'Visible','off','Units','normalized', 'clipping' , 'off');
        titletext = [expname ': GREEN channel'];
        text('Position',[0.35 0.99], 'FontSize',14, 'String',titletext);
end 
sdir = [data_locn 'Analysis/' mouseid '/'];
if ~exist(sdir,'dir'), mkdir(sdir); end
fname_fig = [sdir expname '_comparison_GREEN_channel.fig'];
    savefig( fh, fname_fig );
    saveas( fh, fname_fig(1:end-4), 'png' );
    close( fh );


for ii=0:Nfiles/nPlot
    fh = figure;
    ha = tight_subplot(nRow,nCol,[.01 .01],[.01 .05],[.01 .01]);
    for jj=0:nPlot-1
        if (ii*nPlot+jj+1) <= Nfiles
            axes(ha(ii*nPlot+jj+1));
            imagesc(M(ii*nPlot+jj+1).red.meanregframe./max(max(M(ii*nPlot+jj+1).red.meanregframe))); 
            axis off; colormap(gray);
            str = files(ii*nPlot+jj+1,:);
            title([str(5:6) '-' str(7:8) ' ' str(10:11) ':' str(13:14)]);
        end
    end
    axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0  1],'Box','off',...
        'Visible','off','Units','normalized', 'clipping' , 'off');
        titletext = [expname ': RED channel'];
        text('Position',[0.35 0.99], 'FontSize',14, 'String',titletext);
end 
fname_fig = [sdir expname '_comparison_RED_channel.fig'];
    savefig( fh, fname_fig );
    saveas( fh, fname_fig(1:end-4), 'png' );
    close( fh );


%% Load tracking data for each recording
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
            axes(ha(ii*nPlot+jj+1));
            plot(M(ii*nPlot+jj+1).trackdata.x,M(ii*nPlot+jj+1).trackdata.y); 
            axis off; 
            str = files(ii*nPlot+jj+1,:);
            title([str(5:6) '-' str(7:8) ' ' str(10:11) ':' str(13:14)]);
        end
    end
    axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0  1],'Box','off',...
        'Visible','off','Units','normalized', 'clipping' , 'off');
        titletext = [expname ': trajectories'];
        text('Position',[0.4 0.99], 'FontSize',14, 'String',titletext);
end 
fname_fig = [sdir expname '_comparison_trajectories.fig'];
    savefig( fh, fname_fig );
    saveas( fh, fname_fig(1:end-4), 'png' );
    close( fh );

        
