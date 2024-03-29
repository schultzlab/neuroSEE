% Written by Ann Go
% 
% This function checks whether motion correction OR image registration to a
% template has been done for a file. The check only returns true if the tif
% images and the output data matrix exist.

function check = checkfor_mcorrIm( data_locn, file, mcorr_method, reffile, requireRed )
    if nargin<5, requireRed = true; end
    if nargin<4, reffile = []; end
    if nargin<3, mcorr_method = 'normcorre'; end
    
    if isempty(reffile)
        filedir = [ data_locn 'Data/' file(1:8) '/Processed/' file '/mcorr_' mcorr_method '/' ];
        fname_tif_gr = [filedir file '_2P_XYT_green_mcorr.tif'];
        fname_tif_red = [filedir file '_2P_XYT_red_mcorr.tif'];
        fname_mat = [filedir file '_mcorr_output.mat'];
    else
        filedir = [ data_locn 'Data/' file(1:8) '/Processed/' file '/imreg_' mcorr_method '_ref' reffile '/' ];
        fname_tif_gr = [filedir file '_2P_XYT_green_imreg_ref' reffile '.tif'];
        fname_tif_red = [filedir file '_2P_XYT_red_imreg_ref' reffile '.tif'];
        fname_mat = [filedir file '_imreg_ref' reffile '_output.mat'];
    end

    if requireRed
        if all( [exist(fname_tif_gr,'file'),...
                exist(fname_tif_red,'file'),...
                exist(fname_mat,'file')] )            
            check = true;
        else
            check = false;
        end
    else
        if all( [exist(fname_tif_gr,'file'),...
                exist(fname_mat,'file')] )            
            check = true;
        else
            check = false;
        end
    end
