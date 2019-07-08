% Written by Ann Go
% 
% This function checks whether motion corrected images and the related
% motion correction output mat file exist

function yn = checkfor_mcorrIm( data_locn, file, mcorr_method )
    if strcmpi(mcorr_method,'normcorre')
        filedir = fullfile( data_locn, 'Data/', file(1:8), '/Processed/', file, '/mcorr_normcorre/' );
    else
        filedir = fullfile( data_locn, 'Data/', file(1:8), '/Processed/', file, '/mcorr_fftRigid/' );
    end
    
    if ~exist( filedir, 'dir' ), mkdir( filedir ); end
    fname_tif_gr_mcorr = [filedir file '_2P_XYT_green_mcorr.tif'];
    fname_tif_red_mcorr = [filedir file '_2P_XYT_red_mcorr.tif'];
    fname_mat_mcorr = [filedir file '_mcorr_output.mat'];

    yn_gr_mcorr = exist(fname_tif_gr_mcorr,'file');
    yn_red_mcorr = exist(fname_tif_red_mcorr,'file');
    yn_mat_mcorr = exist(fname_mat_mcorr,'file');

    if all( [yn_gr_mcorr, yn_red_mcorr, yn_mat_mcorr] )
        yn = 1;
    else
        yn = 0;
    end

