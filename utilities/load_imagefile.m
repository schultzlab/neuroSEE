% Written by Ann Go
%
% INPUTS
%   data_locn : directory where GCaMP6 data is
%   file      : part of filename of 2P image in the format
%               yyyymmdd_HH_MM_SS
%   force     : set to 1 to overwrite existing tif files (not motion corrected)
%   suffix    : =[]         to load original raw or tif files
%               ='_mcorr'    to load motion corrected tif files   
%   mcorr_method: motion correction method
%                   1- CaImAn NoRMCorre method
%                   2- fft-rigid method (Katie's)

% OUTPUTS
%   imG       : matrix of green channel image stack
%   imR       : matrix of red channel image stack

function [ imG, imR ] = load_imagefile( data_locn, file, force, suffix, mcorr_method )
    if nargin < 5, mcorr_method  = 'normcorre';    end
    if nargin < 4, suffix  = [];    end
    if nargin < 3, force  = 0;    end
    if strcmpi(mcorr_method,'normcorre'), load_imR = false; end
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
        
        if force % load raw and create tif stacks
            if exist( fname_raw, 'file' )
                options.skipN = 2; %skip every other frame
                imG = read_file( fname_raw, 1, Inf, options );
                    writeTifStack( imG, fname_tif_gr );
                imR = read_file( fname_raw, 2, Inf, options );
                    writeTifStack( imR, fname_tif_red );
                % Note: This is somehow faster than reading entire raw file, saving half of
                % it as imG, the other half as imR and then creating tif stacks.
                err = 0;
            else
                str = sprintf( '%s: No raw file found\n', file );
                cprintf( 'Error', str )
                imG = []; imR = [];
            end
        else
            % Find out if either green or red tif stack exists
            yn_tif_gr = exist(fname_tif_gr,'file');
            yn_tif_red = exist(fname_tif_red,'file');

            % If both exist, load tif stacks
            if yn_tif_gr && yn_tif_red
                imG = read_file( fname_tif_gr );
                imR = read_file( fname_tif_red );
                err = 0;

                % If the sizes of green and red stacks don't match, load raw
                if ~all( size(imG) == size(imR) )
                    if exist( fname_raw, 'file' )
                        options.skipN = 2; %skip every other frame
                        imG = read_file( fname_raw, 1, Inf, options );
                            writeTifStack( imG, fname_tif_gr );
                        imR = read_file( fname_raw, 2, Inf, options );
                            writeTifStack( imR, fname_tif_red );
                    else
                        str = sprintf('%s: tif file sizes dont match. No raw file found\n',file);
                        cprintf('Error',str)
                        imG = []; imR = [];
                        err = 1;
                    end
                end

            % If only green tif exists, load red from raw and create red tif 
            elseif yn_tif_gr && ~yn_tif_red
                imG = read_file( fname_tif_gr );
                if exist(fname_raw,'file')
                    options.skipN = 2; %skip every other frame
                    imR = read_file( fname_raw, 2, Inf, options );
                        writeTifStack( imR, fname_tif_red );
                        err = 0;
                    % If the sizes of green and red stacks don't match, load green from raw
                    if ~all( size(imG) == size(imR) )
                        imG = read_file( fname_raw, 1, Inf, options );
                            writeTifStack( imG, fname_tif_gr );
                    end
                else
                    str = sprintf( '%s: No red tif nor raw files found\n', file );
                    cprintf('Error',str)
                    imG = []; imR = [];
                end

            % If only red tif exists, load green from raw and create green tif   
            elseif ~yn_tif_gr && yn_tif_red
                imR = read_file( fname_tif_red );
                if exist(fname_raw,'file')
                    options.skipN = 2; %skip every other frame
                    imG = read_file( fname_raw, 1, Inf, options );
                        writeTifStack( imG, fname_tif_gr );
                        err = 0;
                    % If the sizes of green and red stacks don't match, load raw
                    if ~all( size(imG) == size(imR) )
                        imR = read_file( fname_raw, 2, Inf, options );
                            writeTifStack( imR, fname_tif_red );
                    end
                else
                    str = sprintf('%s: No green tif nor raw files found\n',file);
                    cprintf('Error',str)
                    imG = []; imR = [];
                end
            % If no tif stacks exists, load all from raw and create tif stacks
            else
                if exist(fname_raw,'file')
                    options.skipN = 2; %skip every other frame
                    imG = read_file( fname_raw, 1, Inf, options );
                        writeTifStack( imG, fname_tif_gr );
                    imR = read_file( fname_raw, 2, Inf, options );
                        writeTifStack( imR, fname_tif_red );
                    % This is somehow faster than reading entire raw file, saving half of
                    % it as imG, the other half as imR and then creating tif stacks.
                    err = 0;
                else
                    str = sprintf('%s: No raw file found\n',file);
                    cprintf('Error',str)
                    imG = []; imR = [];
                end
            end
        end

    else % if there's a suffix '_mcorr'
        % motion corrected tif stacks are in the Processed directory
        if strcmpi(mcorr_method,'normcorre')
            dir_processed = fullfile( data_locn, 'Data/', file(1:8), '/Processed/', file, '/mcorr_normcorre/' );
        else
            dir_processed = fullfile( data_locn, 'Data/', file(1:8), '/Processed/', file, '/mcorr_fftRigid/' );
        end
        fname_tif_gr = [ dir_processed file '_2P_XYT_green' suffix '.tif' ];
        fname_tif_red = [ dir_processed file '_2P_XYT_red' suffix '.tif' ];
        
        imG = read_file( fname_tif_gr );
        if load_imR
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
