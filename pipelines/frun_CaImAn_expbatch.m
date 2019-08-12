function frun_CaImAn_expbatch( expname, list, reffile, slacknotify )

% Written by Ann Go

tic

%% Load module folders and define data directory

test = false;                % flag to use one of smaller files in test folder)
[data_locn,comp,err] = load_neuroSEEmodules(test);
if ~isempty(err)
    beep
    cprintf('Errors',err);    
    return
end
if strcmpi(comp,'hpc')
    maxNumCompThreads(32);      % max # of computational threads, must be the same as # of ncpus specified in jobscript (.pbs file)
end

% Send Ann slack message if processing has started
if slacknotify
    slacktext = [expname ': starting CaImAn'];
    SendSlackNotification('https://hooks.slack.com/services/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', ...
       slacktext, '@m.go', ...
       [], [], [], []);   
end


%% Image files

listfile = [data_locn 'Digital_Logbook/lists/' list];
files = extractFilenamesFromTxtfile( listfile );

for i = 1:size(files,1)
    file = files(i,:);
    if strcmpi(file,reffile)
        fname = [data_locn 'Data/' file(1:8) '/Processed/' file '/mcorr_normcorre/' file '_2P_XYT_green_mcorr.tif'];
    else
        fname = [data_locn 'Data/' file(1:8) '/Processed/' file '/mcorr_normcorre_ref' reffile '/' file '_2P_XYT_green_imreg_ref' reffile '.tif'];
    end
    fprintf('%s: reading %s\n',expname,file)
    Yii = read_file(fname);
    Y(:,:,(i-1)*size(Yii,3)+1:i*size(Yii,3)) = Yii;
end

% Downsample
if size(files,1) <= 7
    k = 5;
elseif size(files,1) <= 10
    k = 7;
elseif size(files,1) <= 13
    k = 9;
else
    k = round( size(files,1)*7420/11000 );
end
imG = Y(:,:,1:k:end);
clear Y


%% ROI segmentation

if str2double(file(1:4)) > 2018
    cellrad = 6;                  % expnameected radius of a cell (pixels)    
    maxcells = 300;     % estimated number of cells in FOV      
else
    cellrad = 9;            
    maxcells = 200;         
end

[~, masks, corr_image] = CaImAn( imG, file, maxcells, cellrad );
clear imG 


%% Output
% ROIs overlayed on correlation image
plotopts.plot_ids = 1; % set to 1 to view the ID number of the ROIs on the plot
fig = plotContoursOnSummaryImage(corr_image, masks, plotopts);

% save masks.mat
for i = 1:size(files,1)
    file = files(i,:);
    if strcmpi(file,reffile)
        sdir = [data_locn 'Data/' file(1:8) '/Processed/' file '/mcorr_normcorre/CaImAn_' expname '/'];
        if ~exist(sdir,'dir'), mkdir(sdir); end
        sname = [sdir expname '_ref' reffile '_masks.mat'];
    else
        sdir = [data_locn 'Data/' file(1:8) '/Processed/' file '/mcorr_normcorre_ref' reffile '/CaImAn_' expname '/'];
        if ~exist(sdir,'dir'), mkdir(sdir); end
        sname = [sdir expname '_ref' reffile '_masks.mat'];
    end
    save(sname,'masks','corr_image');
    savefig(fig,[sdir expname '_ref' reffile '_ROIs']);
    saveas(fig,[sdir expname '_ref' reffile '_ROIs'],'png');
end
sdir = [data_locn 'Analysis/' expname(1:3) '/experiment_ROIs/' expname '_ref' reffile '/'];
sname = [sdir expname '_ref' reffile '_masks.mat'];
save(sname,'masks','corr_image');
savefig(fig,[sdir expname '_ref' reffile '_ROIs']);
saveas(fig,[sdir expname '_ref' reffile '_ROIs'],'png');
close(fig);

% Send Ann slack message if processing has finished
if slacknotify
    slacktext = [expname ': CaImAn FINISHED. No errors!'];
    SendSlackNotification('https://hooks.slack.com/services/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', ...
       slacktext, '@m.go', ...
       [], [], [], []);   
end

t = toc;
str = sprintf('%s: Processing done in %g hrs\n', expname, round(t/3600,2));
cprintf(str)

end
