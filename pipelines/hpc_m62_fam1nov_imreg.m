% Written by Ann Go

clear; % close all;
tic

%% Experiment name
exp = 'm62_fam1nov';

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
SendSlackNotification('https://hooks.slack.com/services/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', ...
   'hpc: m62 fam1nov imreg processing started', '@m.go', ...
   [], [], [], []);      


%% Image files

files = extractFilenamesFromTxtfile('list_m62_fam1nov.txt');
templatefile = '20181016_09_49_06'; % reference template (in fam1fam1rev list)

%% Image registration
max_shift = 50;
try
nonrigid = NoRMCorreSetParms(...
                'd1',size(template,1),...
                'd2',size(template,2),...
                'grid_size',[32,32],...
                'overlap_pre',[32,32],...
                'overlap_post',[32,32],...
                'iter',1,...
                'use_parallel',false,...
                'max_shift',50,...
                'mot_uf',4,...
                'bin_width',200,...
                'max_dev',3,...
                'us_fac',50,...
                'init_batch',200);
[imG, imR, shifts, template, options, col_shift, green, red] = normcorre_2ch( imG, imR, params.nonrigid, template);


for i = 1:size(files,1)
    file = files(i,:);
    fname = [data_locn 'Data/' file(1:8) '/Processed/' file '/mcorr_normcorre_ref20181016_09_49_06/' file '_2P_XYT_green_mcorr_ref20181016_09_49_06.tif'];
    fprintf('%s: reading %s',exp,file)
    Yii = read_file(fname);
    Y(:,:,(i-1)*size(Yii,3)+1:i*size(Yii,3)) = Yii;
end

% Downsample
imG = Y(:,:,1:7:end);
clear Y


%% ROI segmentation

if str2double(file(1:4)) > 2018
    cellrad = 6;                  % expected radius of a cell (pixels)    
    maxcells = 300;     % estimated number of cells in FOV      
else
    cellrad = 9;            
    maxcells = 200;         
end

[~, masks, ~] = CaImAn( imG, file, maxcells, cellrad );
clear imG 


%% Output
for i = 1:size(files,1)
    file = files(i,:);
    if strcmpi(file,ref)
        sdir = [data_locn 'Data/' file(1:8) '/Processed/' file '/mcorr_normcorre/CaImAn_fam1fam1rev/'];
        if ~exist(sdir,'dir'), mkdir(sdir); end
        sname = [sdir 'masks.mat'];
    else
        sdir = [data_locn 'Data/' file(1:8) '/Processed/' file '/mcorr_normcorre_ref20181016_09_49_06/CaImAn_fam1fam1rev/'];
        if ~exist(sdir,'dir'), mkdir(sdir); end
        sname = [sdir 'masks.mat'];
    end
    save(sname,'masks')
end

% Send Ann slack message if processing has finished
SendSlackNotification('https://hooks.slack.com/services/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', ...
   'hpc: m62 fam1fam1rev CaImAn processing FINISHED', '@m.go', ...
   [], [], [], []);      

t = toc;
str = sprintf('%s: Processing done in %g hrs\n', file, round(t/3600,2));
cprintf(str)

