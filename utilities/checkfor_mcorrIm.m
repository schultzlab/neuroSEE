% Written by Ann Go
% 
% This function checks whether motion correction OR image registration to a
% template has been done for a file. The check only returns true if the tif
% images and the output data matrix exist.

function check = checkfor_mcorrIm( data_locn, file, mcorr_method, reffile )
    if nargin<4, reffile = []; end
    if nargin<3, mcorr_method = 'normcorre-nr'; end
    
    if isempty(reffile)
        filedir = [ data_locn 'Data/' file(1:8) '/Processed/' file '/mcorr_' mcorr_method '/' ];
            if ~exist( filedir, 'dir' ), mkdir( filedir ); end
        fname_tif_gr = [filedir file '_2P_XYT_green_mcorr.tif'];
        fname_tif_red = [filedir file '_2P_XYT_red_mcorr.tif'];
        fname_mat = [filedir file '_mcorr_output.mat'];
    else
        filedir = [ data_locn 'Data/' file(1:8) '/Processed/' file '/mcorr_' mcorr_method '_ref' reffile '/' ];
            if ~exist( filedir, 'dir' ), mkdir( filedir ); end
        fname_tif_gr = [filedir file '_2P_XYT_green_imreg_ref' reffile '.tif'];
        fname_tif_red = [filedir file '_2P_XYT_red_imreg_ref' reffile '.tif'];
        fname_mat = [filedir file '_imreg_ref' reffile '_output.mat'];
    end

    if all( [exist(fname_tif_gr,'file'),...
            exist(fname_tif_red,'file'),...
            exist(fname_mat,'file')] )
        check = 1;
    else
        check = 0;
    end
