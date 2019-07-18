% Written by Ann Go
%
% INPUTS
%   data_locn : directory where GCaMP6 data is
%   file      : part of filename of 2P image in the format
%               yyyymmdd_HH_MM_SS
%   force     : set to 1 to overwrite existing tif files (not motion corrected)
%   suffix    : =[]          to load original raw or tif files
%               ='_mcorr'    to load motion corrected tif files   
%   params.mcorr_method: ['normcorre' or 'fftRigid'] motion correction method
%                       * CaImAn NoRMCorre method OR fft-rigid method (Katie's)
%   params.segment_method: ['CaImAn' or 'ABLE'] roi segmentation method

% OUTPUTS
%   imG       : matrix of green channel image stack
%   imR       : matrix of red channel image stack

function [ imG, imR ] = load_imagefile( data_locn, file, force, suffix, params )
    
    if nargin < 5 
        mcorr_method  = 'normcorre'; 
        segment_method  = 'CaImAn'; 
    else
        mcorr_method  = params.methods.mcorr_method;
        segment_method  = params.methods.segment_method;
    end
    if nargin < 4, suffix  = [];    end
    if nargin < 3, force  = 0;    end
    if strcmpi(segment_method,'CaImAn'), load_imR = false; else, load_imR = true; end
    err = 1;

    str = sprintf( '%s: Loading %s images...\n', file, suffix );
    cprintf( str )
    
    % Define file names
    if isempty( suffix )
        % original tif files are in the 2P directory
        dir_2P = [ data_locn 'Data/' file(1:8) '/2P/' file '_2P/' ];
        fname_tif_gr = [ dir_2P file '_2P_XYT_green.tif' ];
        fname_tif_red = [ dir_2P file '_2P_XYT_red.tif' ];
        
        fname_raw = [ dir_2P file '_2P_XYT.raw' ];
        if ~exist(fname_raw,'file')
            fname_raw = [ dir_2P file '_2P_XYTZ.raw' ];
        end
        
        if any([ force, ~exist(fname_tif_gr,'file'), ~exist(fname_tif_red,'file') ]) 
            % load raw and create tif stacks
            if exist( fname_raw, 'file' )
                options.skipN = 2; %skip every other frame
                if force || ~exist(fname_tif_gr,'file')
                    imG = read_file( fname_raw, 1, Inf, options );
                    writeTifStack( imG, fname_tif_gr );
                end
                if force || ~exist(fname_tif_red,'file')
                    imR = read_file( fname_raw, 2, Inf, options );
                    writeTifStack( imR, fname_tif_red );
                end
                % Note: This is somehow faster than reading entire raw file, saving half of
                % it as imG, the other half as imR and then creating tif stacks.
                err = 0;
            else
                str = sprintf( '%s: No raw file found\n', file );
                cprintf( 'Error', str )
                imG = []; imR = [];
            end
        else
            imG = read_file( fname_tif_gr );
            imR = read_file( fname_tif_red ); 
            err = 0;
        end
    else % if there's a suffix '_mcorr'
        % motion corrected tif stacks are in the Processed directory
        if strcmpi(mcorr_method,'normcorre')
            dir_processed = fullfile( data_locn, 'Data/', file(1:8), '/Processed/', file, '/mcorr_normcorre/' );
        else
            dir_processed = fullfile( data_locn, 'Data/', file(1:8), '/Processed/', file, '/mcorr_fftRigid/' );
        end
        fname_tif_gr = [ dir_processed file '_2P_XYT_green' suffix '.tif' ];
        
        imG = read_file( fname_tif_gr );
       
        if load_imR
            fname_tif_red = [ dir_processed file '_2P_XYT_red' suffix '.tif' ];
            imR = read_file( fname_tif_red );
        else
            imR = [];
        end
        err = 0;
    end

    if ~err
        newstr = sprintf('%s: %s Images loaded\n',file,suffix);
        refreshdisp(  newstr,str );
    end
