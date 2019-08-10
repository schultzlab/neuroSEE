function frun_comparison_imreg_expbatch( expname, reffile, list, mcorr_method, force )

% e.g. expname: 'm62_fam1nov_fam1fam1rev'
%               'm62_fam1fam1rev'
%               'm70_fam1_day1-5'
%               'm82_open_day1-2'
%               'm82_open_day1'

if nargin<5, force = false; end
if nargin<4, mcorr_method = 'normcorre'; end

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

[nRow, nCol] = getnRownCol(Nfiles);
nPlot = nCol*nRow;

% template file
filedir = [data_locn 'Data/' reffile(1:8) '/Processed/' reffile '/mcorr_' mcorr_method '/'];
c = load([filedir reffile '_mcorr_output.mat']);
templateG = c.template./max(max(c.template));
templateR = c.red.meanregframe./(max(max(c.red.meanregframe)));
clear c


%% Load image data for each recording
sdir = [data_locn 'Analysis/' mouseid '/individual file comparisons/' expname '_imreg_ref' reffile '/'];
if ~exist(sdir,'dir'), mkdir(sdir); end
    
if any([ force, ~exist([sdir expname '_GREEN_imreg_ref' reffile '.fig'],'file'),...
                ~exist([sdir expname '_RED_imreg_ref' reffile '.fig'],'file') ])
    for i = 1:Nfiles
        file = files(i,:);
        if ~strcmpi(file,reffile)
            fname = [data_locn 'Data/' file(1:8) '/Processed/' file '/mcorr_normcorre_ref' reffile '/'...
                     file '_imreg_ref' reffile '_output.mat'];
            c = load(fname);
            M(i).green = c.green;
            if isfield(c,'red')
                M(i).red = c.red;
            end
        else
            M(i).green.meanregframe = templateG;
            M(i).red.meanregframe = templateR;
        end
    end
    clear c

    %% Compare green and red images to see which ones are most similar

    for ii=0:Nfiles/nPlot
        fh = figure;
        ha = tight_subplot(nRow,nCol,[.01 .01],[.01 .05],[.01 .01]);
        for jj=0:nPlot-1
            if (ii*nPlot+jj+1) <= Nfiles
                axes(ha(ii*nPlot+jj+1));
                C = imfuse(M(ii*nPlot+jj+1).green.meanregframe./max(max(M(ii*nPlot+jj+1).green.meanregframe)),...
                            templateG,'falsecolor','Scaling','joint','ColorChannels',[1 2 0]);
                imshow(C);
                axis off; 
                str = files(ii*nPlot+jj+1,:);
                title([str(5:6) '-' str(7:8) ' ' str(10:11) ':' str(13:14)]);
            end
        end
        axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0  1],'Box','off',...
            'Visible','off','Units','normalized', 'clipping' , 'off');
            titletext = [expname(5:end) ': GREEN channel registered to '...
                             reffile(1:8) ' ' reffile(10:11) ':' reffile(13:14) ':' reffile(16:17)];
            text('Position',[0.35 0.99], 'FontSize',14, 'String',titletext);
    end 
    fname_fig = [sdir expname '_GREEN_imreg_ref' reffile '.fig'];
        savefig( fh, fname_fig );
        saveas( fh, fname_fig(1:end-4), 'png' );
        close( fh );

    if isfield(M(1),'red')
        for ii=0:Nfiles/nPlot
            fh = figure;
            ha = tight_subplot(nRow,nCol,[.01 .01],[.01 .05],[.01 .01]);
            for jj=0:nPlot-1
                if (ii*nPlot+jj+1) <= Nfiles
                    axes(ha(ii*nPlot+jj+1));
                    C = imfuse(M(ii*nPlot+jj+1).red.meanregframe./max(max(M(ii*nPlot+jj+1).red.meanregframe)),...
                                templateR,'falsecolor','Scaling','joint','ColorChannels',[1 2 0]);
                    imshow(C);
                    axis off; 
                    str = files(ii*nPlot+jj+1,:);
                    title([str(5:6) '-' str(7:8) ' ' str(10:11) ':' str(13:14)]);
                end
            end
            axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0  1],'Box','off',...
                'Visible','off','Units','normalized', 'clipping' , 'off');
                titletext = [mouseid '_' expname(5:end) ': RED channel registered to '...
                                 reffile(1:8) ' ' reffile(10:11) ':' reffile(13:14) ':' reffile(16:17)];
                text('Position',[0.35 0.99], 'FontSize',14, 'String',titletext);
        end 
        fname_fig = [sdir expname '_RED_imreg_ref' reffile '.fig'];
            savefig( fh, fname_fig );
            saveas( fh, fname_fig(1:end-4), 'png' );
            close( fh );
    end
end


