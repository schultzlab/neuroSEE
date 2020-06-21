% Written by Ann Go
%
% This script collates the image registration results for the files in
% 'list'. Comparison figures are saved in
%   /...thefarm2/CrazyEights/AD_2PCa/Analysis/m##/m##_expname/individual_proc/
%
% INPUTS
% list      : name of text file containing filenames of files to be compared.
%               Typically in the format 'list_m##_expname.txt'.
%   e.g.    'list_m62_fam1nov-fam1.txt'         - all fam1 files in fam1nov experiment
%           'list_m62_fam1nov.txt'              - all files in fam1nov experiments
%           'list_m79_fam1_s1-5.txt'            - all fam1 files across 5 sessions           
%           'list_m86_open_s1-2.txt'            - all open field files across 2 sessions
% reffile   : (optional) file used as registration template. If none
%               specified, the first file in the list is used
% imreg_method : (optional) image registration method (default: normcorre)
%   e.g.    'normcorre'                         - rigid + non-rigid normcorre
%           'normcorre-nr'                      - non-rigid normcorre
%           'normcorre-r'                       - rigid normcorre
%           'fftRigid'                          - Katie's method
% force     : (optional) flag to force generation of comparison figures even though they
%               already exist (default: false)
% mcorr_method : motion correction method used for correcting reffile
%                (optional, default: same as imreg_method)

function frun_collate_imreg_results( list, imreg_method, reffile, force, mcorr_method )

if nargin<2, imreg_method = 'normcorre'; end
% if nargin<3 see line 44
if nargin<4, force = false; end
if nargin<5, mcorr_method = imreg_method; end

%% Load module folders and define data directory
[data_locn,~,err] = load_neuroSEEmodules;
if ~isempty(err)
    beep
    cprintf('Errors',err);    
    return
end

%% MouseID and experiment name
[ mouseid, expname ] = find_mouseIDexpname(list);

%% Files
listfile = [data_locn 'Digital_Logbook/lists_imaging/' list];
files = extractFilenamesFromTxtfile( listfile );
if nargin<2, reffile = files(1,:); end
Nfiles = size(files,1);

[nRow, nCol] = getnRownCol(Nfiles);
nPlot = nCol*nRow;
Nfig = (Nfiles/nPlot)-1;
if Nfig < 0, Nfig = 0; end

%% Template file
filedir = [data_locn 'Data/' reffile(1:8) '/Processed/' reffile '/mcorr_' mcorr_method '/'];
c = load([filedir reffile '_mcorr_output.mat']);
templateG = c.green.meanregframe./(max(max(c.green.meanregframe)));
templateR = c.red.meanregframe./(max(max(c.red.meanregframe)));
clear c

%% Load image data for each recording
% directory where registration results summary is saved: sdir
if strcmpi(imreg_method, mcorr_method)
    sdir = [data_locn 'Analysis/' mouseid '/' mouseid '_' expname '/individual_proc/imreg_' imreg_method '/indiv_imreg_ref' reffile '/'];
else
    sdir = [data_locn 'Analysis/' mouseid '/' mouseid '_' expname '/individual_proc/imreg_' imreg_method '/indiv_imreg_ref' reffile '_' mcorr_method '/'];
end
if ~exist(sdir,'dir'), mkdir(sdir); end
    
if any([ force, ~exist([sdir mouseid '_' expname '_GREEN_imreg_ref' reffile '.fig'],'file'),...
                ~exist([sdir mouseid '_' expname '_RED_imreg_ref' reffile '.fig'],'file') ])
    for i = 1:Nfiles
        file = files(i,:);
        if ~strcmpi(file,reffile)
            % individual image registration result
            if strcmpi(imreg_method, mcorr_method)
                fname = [data_locn 'Data/' file(1:8) '/Processed/' file '/imreg_' imreg_method '_ref' reffile '/'...
                         file '_imreg_ref' reffile '_output.mat'];
            else
                fname = [data_locn 'Data/' file(1:8) '/Processed/' file '/imreg_' imreg_method '_ref' reffile '_' mcorr_method '/'...
                         file '_imreg_ref' reffile '_output.mat'];
            end
            if exist(fname,'file')
                c = load(fname);
                M(i).green = c.green;
                if isfield(c,'red')
                    M(i).red = c.red;
                else
                    M(i).red = [];
                end
            else
                M(i).green = [];
                M(i).red = [];
            end
        else
            M(i).green.meanregframe = templateG;
            M(i).red.meanregframe = templateR;
        end
    end
    clear c

    %% Compare green and red images to see which ones are most similar

    for ii=0:Nfig
        fh = figure;
        ha = tight_subplot(nRow,nCol,[.01 .01],[.01 .05],[.01 .01]);
        for jj=0:nPlot-1
            if (ii*nPlot+jj+1) <= Nfiles
                if ~isempty(M(ii*nPlot+jj+1).green)
                    axes(ha(ii*nPlot+jj+1));
                    C = imfuse(M(ii*nPlot+jj+1).green.meanregframe./max(max(M(ii*nPlot+jj+1).green.meanregframe)),...
                                templateG,'falsecolor','Scaling','joint','ColorChannels',[1 2 0]);
                    imshow(C);
                    axis off; 
                    str = files(ii*nPlot+jj+1,:);
                    title([str(5:6) '-' str(7:8) ' ' str(10:11) ':' str(13:14)]);
                end
            end
        end
        axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0  1],'Box','off',...
            'Visible','off','Units','normalized', 'clipping' , 'off');
            titletext = [mouseid '_' expname ': GREEN channel registered to '...
                             reffile];
            ind = strfind(titletext,'_');
            titletext(ind) = '-';
            text('Position',[0.35 0.99], 'FontSize',14, 'String',titletext);
    end 
    fname_fig = [sdir mouseid '_' expname '_GREEN_imreg_ref' reffile '.fig'];
        savefig( fh, fname_fig );
        saveas( fh, fname_fig(1:end-4), 'png' );
        close( fh );

    if isfield(M(1),'red')
        for ii=0:Nfig
            fh = figure;
            ha = tight_subplot(nRow,nCol,[.01 .01],[.01 .05],[.01 .01]);
            for jj=0:nPlot-1
                if (ii*nPlot+jj+1) <= Nfiles
                    if ~isempty(M(ii*nPlot+jj+1).red)
                        axes(ha(ii*nPlot+jj+1));
                        C = imfuse(M(ii*nPlot+jj+1).red.meanregframe./max(max(M(ii*nPlot+jj+1).red.meanregframe)),...
                                    templateR,'falsecolor','Scaling','joint','ColorChannels',[1 2 0]);
                        imshow(C);
                        axis off; 
                        str = files(ii*nPlot+jj+1,:);
                        title([str(5:6) '-' str(7:8) ' ' str(10:11) ':' str(13:14)]);
                    end
                end
            end
            axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0  1],'Box','off',...
                'Visible','off','Units','normalized', 'clipping' , 'off');
                titletext = [mouseid '_' expname ': RED channel registered to ' reffile];
                ind = strfind(titletext,'_');
                titletext(ind) = '-';
                text('Position',[0.35 0.99], 'FontSize',14, 'String',titletext);
        end 
        fname_fig = [sdir mouseid '_' expname '_RED_imreg_ref' reffile '.fig'];
            savefig( fh, fname_fig );
            saveas( fh, fname_fig(1:end-4), 'png' );
            close( fh );
    end
end


