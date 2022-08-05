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
        dtsG(i,:) = deltaf_result.deltaf_result.(['cell' num2str(i-1)]).trial0(1,:);
    end

output.dtsG = dtsG;
output.ddf_f = ddf_f;

save([fissadir 'm62_fov2_fam1fam2-fam1_ref20181013_10_53_51_fissa_output.mat'],'-struct','output');

% dtsG = tsG; ddf_f = df_f; spikes = tsG;
% GUI_viewtimeseries(tsG, df_f, dtsG, ddf_f, spikes)
