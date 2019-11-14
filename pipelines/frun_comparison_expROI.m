function frun_comparison_expROI( expname, list, reffile, force )

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

%% Load ROIs for each file
sdir = [data_locn 'Analysis/' mouseid '/experiment_ROIs/' expname '_ref' reffile '/'];
    if ~exist(sdir,'dir'), mkdir(sdir1); end
    
if force || ~exist([sdir expname '_compareROIs.fig'],'file')
    for i = 1:Nfiles
        file = files(i,:);
        fname = [data_locn 'Data/' file(1:8) '/Processed/' file '/mcorr_normcorre/' file '_segment_output.tif'];
        if exist(fname,'file')
            c = load(fname);
            M(i).green = c.green;
            M(i).red = c.red;
        end
    end
    clear c
end


