[data_locn,comp,err] = load_neuroSEEmodules;
list = 'list_m62_fov2_fam1fam2-fam1.txt';
listfile = [data_locn 'Digital_Logbook/lists_imaging/' list];
files = extractFilenamesFromTxtfile( listfile );

tifflist = [];
for n = 1:size(files,1)
    file = files(n,:);
    tiffile = [data_locn 'Data/' file(1:8) '/Processed/' file '/mcorr_normcorre/' file '_2P_XYT_green_mcorr.tif'];
    if n>1
        tifflist = [tifflist ',' tiffile];
    else
        tifflist = tiffile;
    end
end

load([data_locn 'Analysis/m62/fov2/m62_fov2_fam1fam2-fam1/group_proc/imreg_normcorre_CaImAn/m62_fov2_fam1fam2-fam1_imreg_ref20181013_10_53_51_concenvrois/m62_fov2_fam1fam2-fam1_ref20181013_10_53_51_segment_output.mat'])
fissadir = [data_locn 'Analysis/m62/fov2/m62_fov2_fam1fam2-fam1/group_proc/imreg_normcorre_CaImAn/m62_fov2_fam1fam2-fam1_imreg_ref20181013_10_53_51_concenvrois/FISSA/'];

runFISSA( masks, tifflist, fissadir )

result = load([fissadir 'FISSAout/matlab.mat'],'result');
deltaf_result = load([fissadir 'FISSAout/matlab.mat'],'deltaf_result');

dtsG = zeros(size(masks,3),size(result.result.cell0.trial0,2));
    for i = 1:numel(fieldnames(result.result))
        dtsG(i,:) = result.result.(['cell' num2str(i-1)]).trial0(1,:);
    end

ddf_f = zeros(size(masks,3),size(deltaf_result.deltaf_result.cell0.trial0,2));
    for i = 1:numel(fieldnames(deltaf_result.deltaf_result))
        ddf_f(i,:) = deltaf_result.deltaf_result.(['cell' num2str(i-1)]).trial0(1,:);
    end

output.dtsG = dtsG;
output.ddf_f = ddf_f;

save([fissadir 'm62_fov2_fam1fam2-fam1_ref20181013_10_53_51_fissa_output.mat'],'-struct','output');

% dtsG = tsG; ddf_f = df_f; spikes = tsG;
% GUI_viewtimeseries(tsG, df_f, dtsG, ddf_f, spikes)

%% fissa single file test
[data_locn,comp,err] = load_neuroSEEmodules;
file = '20181011_14_54_19';
tiffile = [data_locn 'Data/' file(1:8) '/Processed/' file '/mcorr_normcorre/' file '_2P_XYT_green_mcorr.tif'];
roizip = [data_locn 'Data/20181011/Processed/20181011_14_54_19/imreg_normcorre_ref20181011_15_10_39/CaImAn_m62_fov1_fam1fam2-fam1/FISSA/rois.zip'];
outdir = '/Users/mgo/Desktop/FISSAout/';

mydir  = pwd;
ind   = strfind(mydir,'/');
newdir = mydir(1:ind(end)-1);
folder = fullfile(newdir(1:ind(end)-1),'/neuropil_decontamination');
pyfun = [folder '/runFISSA.py' ' ' tiffile ' ' roizip ' ' outdir];
user = mydir(ind(2)+1:ind(3)-1);
if exist(['/Users/' user '/anaconda3/envs/neuroSEE/'],'dir')
    python_executable = ['/Users/' user '/anaconda3/envs/neuroSEE/bin/python'];
elseif exist(['/home/' user '/anaconda3/envs/neuroSEE/'],'dir')
    python_executable = ['/home/' user '/anaconda3/envs/neuroSEE/bin/python'];
else
    python_executable = ['/rds/general/user/' user '/home/anaconda3/envs/neuroSEE/bin/python'];
end
pystr = [python_executable ' ' pyfun];
system( pystr )

spikes = zeros(size(tsG));
GUI_viewtimeseries(tsG, df_f, dtsG, ddf_f, spikes)

dtsG = zeros(size(tsG,1),size(result.cell0.trial0,2));
    for i = 1:numel(fieldnames(result))
        dtsG(i,:) = result.(['cell' num2str(i-1)]).trial0(1,:);
    end

ddf_f = zeros(size(tsG,1),size(result.cell0.trial0,2));
    for i = 1:numel(fieldnames(result))
        ddf_f(i,:) = deltaf_result.(['cell' num2str(i-1)]).trial0(1,:);
    end

    ddf_f = zeros(size(dtsG));
    ddf_prctile = 5;
    for i = 1:size(dtsG,1)
        x = lowpass( medfilt1(dtsG(i,:), 17), 1, 30.9 );
        fo = ones(size(x)) * prctile(x,ddf_prctile);
        while fo == 0
            fo = ones(size(x)) * prctile(x,ddf_prctile+5);
            ddf_prctile = ddf_prctile+5;
        end
        ddf_f(i,:) = (x - fo) ./ fo;
    end
