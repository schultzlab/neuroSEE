% Written by Ann Go
%
% INPUTS
%   data_locn : directory where GCaMP6 data is
%   file      : part of filename of 2P image in the format
%               yyyymmdd_HH_MM_SS
%   forceRaw  : flag to overwrite existing tif files (not motion corrected)
%               and read raw file
%   suffix    : = []         to load original raw or tif files
%               ='_mcorr'    to load motion corrected tif files   
%   mcorr_method : ['normcorre' or 'fftRigid'] motion correction method
%                       * CaImAn NoRMCorre method OR fft-rigid method (Katie's)
%   load_imR  : flag to load imR

% OUTPUTS
%   imG       : matrix of green channel image stack
%   imR       : matrix of red channel image stack

function [ imG, imR ] = load_imagefile( data_locn, file, forceRaw, suffix, mcorr_method, load_imR, reffile, imreg_method )
    if nargin < 3, forceRaw = false; end
    if nargin < 4, suffix  = []; end
    if nargin < 5, mcorr_method  = 'normcorre'; end
    if nargin < 6, load_imR = true; end
    if nargin < 7, reffile = []; imreg_method = []; end
    if nargin < 8 || isempty(imreg_method)
        if ~isempty(mcorr_method)
            imreg_method = mcorr_method; 
        else
            imreg_method = [];
        end
    end

    if strcmpi(suffix, '_imreg') && isempty(reffile)
        beep
        cprintf('Errors','Reference file for registered image required!\n');
        return
    end
    str = sprintf( '%s: Loading %s images...\n', file, suffix );
    cprintf( str )
    
    % Define file names
    % NON-MOTION-CORRECTED
    if isempty( suffix )
        % original tif files are in the 2P directory
        dir_2P = [ data_locn 'Data/' file(1:8) '/2P/' file '_2P/' ];
        fname_tif_gr = [ dir_2P file '_2P_XYT_green.tif' ];
        fname_tif_red = [ dir_2P file '_2P_XYT_red.tif' ];
        
        fname_raw = [ dir_2P file '_2P_XYT.raw' ];
        if ~exist(fname_raw,'file')
            fname_raw = [ dir_2P file '_2P_XYTZ.raw' ];
        end
        
        if ~any([ forceRaw, ~exist(fname_tif_gr,'file'), ~exist(fname_tif_red,'file') ])
            imG = read_file( fname_tif_gr );
            imR = read_file( fname_tif_red ); 
            if any( size(imG) ~= size(imR) )
                forceRaw = true;
            else
                err = 0;
            end
        end
        
        if any([ forceRaw, ~exist(fname_tif_gr,'file'), ~exist(fname_tif_red,'file') ]) 
            % load raw and create tif stacks
            if exist( fname_raw, 'file' )
                options.skipN = 2; %skip every other frame
                if forceRaw 
                    imG = read_file( fname_raw, 1, Inf, options );
                    writeTifStack( imG, fname_tif_gr );
                    imR = read_file( fname_raw, 2, Inf, options );
                    writeTifStack( imR, fname_tif_red );
                end 
                if ~exist(fname_tif_gr,'file')
                    imG = read_file( fname_raw, 1, Inf, options );
                    writeTifStack( imG, fname_tif_gr );
                    imR = read_file( fname_tif_red ); 
                    if any( size(imG) ~= size(imR) )
                        imR = read_file( fname_raw, 2, Inf, options );
                        writeTifStack( imR, fname_tif_red );
                    end
                end
                if  ~exist(fname_tif_red,'file')
                    imG = read_file( fname_tif_gr );
                    imR = read_file( fname_raw, 2, Inf, options );
                    writeTifStack( imR, fname_tif_red );
                    if any( size(imG) ~= size(imR) )
                        imG = read_file( fname_raw, 1, Inf, options );
                        writeTifStack( imG, fname_tif_gr );
                    end

                end
                % Note: This is somehow faster than reading entire raw file, saving half of
                % it as imG, the other half as imR and then creating tif stacks.
                err = 0;
            else
                str = sprintf( '%s: No raw file found\n', file );
                cprintf( 'Error', str )
                imG = []; imR = [];
            end
        end
    elseif strcmpi(suffix, '_mcorr')  
        % MOTION-CORRECTED
        % motion corrected tif stacks are in the Processed directory
        if strcmpi(mcorr_method,'normcorre')
            dir_processed = fullfile( data_locn, 'Data/', file(1:8), '/Processed/', file, '/mcorr_normcorre/' );
        elseif strcmpi(mcorr_method,'normcorre-r')
            dir_processed = fullfile( data_locn, 'Data/', file(1:8), '/Processed/', file, '/mcorr_normcorre-r/' );
        elseif strcmpi(mcorr_method,'normcorre-nr')
            dir_processed = fullfile( data_locn, 'Data/', file(1:8), '/Processed/', file, '/mcorr_normcorre-nr/' );
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
    elseif strcmpi(suffix, '_imreg')
        % REGISTERED TO A REFERENCE FILE
        if strcmpi(imreg_method,'normcorre')
            dir_processed = fullfile( data_locn, 'Data/', file(1:8), '/Processed/', file, '/imreg_normcorre_', reffile, '/' );
        elseif strcmpi(imreg_method,'normcorre-r')
            dir_processed = fullfile( data_locn, 'Data/', file(1:8), '/Processed/', file, '/imreg_normcorre-r_', reffile, '/' );
        elseif strcmpi(imreg_method,'normcorre-nr')
            dir_processed = fullfile( data_locn, 'Data/', file(1:8), '/Processed/', file, '/imreg_normcorre-nr_', reffile, '/' );
        else
            dir_processed = fullfile( data_locn, 'Data/', file(1:8), '/Processed/', file, '/imreg_fftRigid_', reffile, '/' );
        end
        fname_tif_gr = [ dir_processed file '_2P_XYT_green' suffix '_' reffile '.tif' ];

        imG = read_file( fname_tif_gr );

        if load_imR
            fname_tif_red = [ dir_processed file '_2P_XYT_red' suffix '_' reffile '.tif' ];
            imR = read_file( fname_tif_red );
        else
            imR = [];
        end
        err = 0;
    end

    if ~err
        newstr = sprintf('%s: %s Images loaded\n', file, suffix);
        refreshdisp(  newstr,str );
    end
