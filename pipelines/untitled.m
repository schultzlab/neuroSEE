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
% dtsG = tsG; ddf_f = df_f; spikes = tsG;
% GUI_viewtimeseries(tsG, df_f, dtsG, ddf_f, spikes)

