[data_locn,comp,err] = load_neuroSEEmodules;
list = 'list_m62_fov1_fam1fam2-fam1.txt';
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

roizip = '/Users/mgo/Desktop/rois.zip';
fissadir = '/Users/mgo/Desktop/FISSAout/';

mydir  = pwd;
ind   = strfind(mydir,'/');
newdir = mydir(1:ind(end)-1);
folder = fullfile(newdir(1:ind(end)-1),'/neuropil_decontamination');
pyfun = [folder '/runFISSA.py' ' ' tifflist ' ' roizip ' ' fissadir];
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

result = load([fissadir 'FISSAout/matlab.mat'],'result');
deltaf_result = load([fissadir 'FISSAout/matlab.mat'],'deltaf_result');

dtsG = zeros(size(masks,3),size(result.result.cell0.trial0,2));

result.result = result;
dtsG = zeros(184,size(result.result.cell0.trial0,2));
t = size(result.result.cell0.trial0,2);
for i = 1:numel(fieldnames(result.result))
    for j = 1:numel(fieldnames(result.result.cell0))
        dtsG(i,(j-1)*t+1:j*t) = result.result.(['cell' num2str(i-1)]).(['trial' num2str(j-1)])(1,:);
    end
end

df_result.deltaf_result = deltaf_result;
ddf_f = zeros(184,size(df_result.deltaf_result.cell0.trial0,2));
for i = 1:numel(fieldnames(df_result.deltaf_result))
    for j = 1:numel(fieldnames(df_result.deltaf_result.cell0))
        ddf_f(i,(j-1)*t+1:j*t) = lowpass( medfilt1(df_result.deltaf_result.(['cell' num2str(i-1)]).(['trial' num2str(j-1)])(1,:), 17), 1, 30.9 );
    end
end

output.dtsG = dtsG;
output.ddf_f = ddf_f;

save([fissadir 'm62_fov1_fam1fam2-fam1_ref20181011_15_10_39_fissa_output.mat'],'-struct','output');

dtsG = tsG; ddf_f = df_f; spikes = tsG;
tsG = dtsG; df_f = ddf_f; 
spikes = zeros(size(tsG));
GUI_viewtimeseries(tsG, df_f, dtsG, ddf_f, spikes)

%% fissa single file test
[data_locn,comp,err] = load_neuroSEEmodules;
file = '20181011_14_54_19';
tiffile = [data_locn 'Data/' file(1:8) '/Processed/' file '/mcorr_normcorre/' file '_2P_XYT_green_mcorr.tif'];
roizip = '/Users/mgo/Desktop/rois.zip';
fissadir = '/Users/mgo/Desktop/FISSAout/';

mydir  = pwd;
ind   = strfind(mydir,'/');
newdir = mydir(1:ind(end)-1);
folder = fullfile(newdir(1:ind(end)-1),'/neuropil_decontamination');
pyfun = [folder '/runFISSA.py' ' ' tiffile ' ' roizip ' ' fissadir];
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
GUI_viewtimeseries(tsG(:,1:7420), df_f(:,1:7420), dtsG, ddf_f, spikes(:,1:7420))
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
