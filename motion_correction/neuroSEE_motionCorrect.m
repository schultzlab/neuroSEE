% Written by Ann Go
%
% This function implements motion correction on red and green channel
% images
%
% INPUTS
%   imG         : matrix of green image stack
%   imR         : matrix of red image stack
%   data_locn   : GCaMP data repository
%   file        : part of file name of image stacks in the format
%                   yyyymmdd_HH_MM_SS
%   mcorr_method: ['normcorre', 'normcorre-r', 'normcorre-nr' or 'fftRigid'] motion correction method
%                   * CaImAn NoRMCorre method OR fft-rigid method (Katie's)
%   force       : (optional, default: 0) if =1, motion correction will be done even though motion
%                   corrected images already exist
%   reffile     : if registering image stack to a file, reffile is the
%                   template. This file must have already been motion corrected.
% OUTPUTS
%   imG         : matrix of motion corrected green image stack
%   imR         : matrix of motion corrected red image stack
%   mcorr_output: cell array containing
%                   green.[ meanframe, meanregframe ]
%                   red.[ meanframe, meanregframe ]
%                   shifts.[ zipper_shift, shifts ]
%                   template
%   params_mcorr: parameters for specific motion correction method

function [ imG, imR, mcorr_output, params_mcorr ] = neuroSEE_motionCorrect(...
                                                imG, imR, data_locn, file, mcorr_method, params_mcorr, reffile, force )
    
    if nargin<8, force = false; end
    if nargin<7, reffile = []; end
    refChannel = params_mcorr.refChannel;
    
    if isempty(reffile)
        filedir = [ data_locn 'Data/' file(1:8) '/Processed/' file '/mcorr_' mcorr_method '/' ];
            if ~exist( filedir, 'dir' ), mkdir( filedir ); end
        fname_tif_gr_mcorr = [filedir file '_2P_XYT_green_mcorr.tif'];
        fname_tif_red_mcorr = [filedir file '_2P_XYT_red_mcorr.tif'];
        fname_mat_mcorr = [filedir file '_mcorr_output.mat'];
        fname_fig = [filedir file '_mcorr_summary.fig'];
    else
        filedir = [ data_locn 'Data/' file(1:8) '/Processed/' file '/mcorr_' mcorr_method '_ref' reffile '/' ];
            if ~exist( filedir, 'dir' ), mkdir( filedir ); end
        fname_tif_gr_mcorr = [filedir file '_2P_XYT_green_imreg_ref' reffile '.tif'];
        fname_tif_red_mcorr = [filedir file '_2P_XYT_red_imreg_ref' reffile '.tif'];
        fname_mat_mcorr = [filedir file '_imreg_ref' reffile '_output.mat'];
        fname_fig = [filedir file '_imreg_ref' reffile '_summary.fig'];
    end

    if any([ force, ~exist(fname_tif_gr_mcorr,'file'), ~exist(fname_tif_red_mcorr,'file'), ~exist(fname_mat_mcorr,'file') ])
        if isempty(reffile)
            str = sprintf( '%s: Starting motion correction\n', file );
            template = [];
        else
            str = sprintf( '%s: Starting image registration to %s\n', file, reffile );
            refdir = [data_locn 'Data/' reffile(1:8) '/Processed/' reffile '/mcorr_' mcorr_method '/'];
            c = load([refdir reffile '_mcorr_output.mat']);
            if strcmpi(refChannel,'green')
                template = c.green.meanregframe;
            else
                template = c.red.meanregframe;
            end
        end
        cprintf( 'Text', str );
        
        if strcmpi(mcorr_method,'normcorre')
            mcorr_output.params.normcorre_r = params_mcorr.normcorre_r;
            mcorr_output.params.normcorre_r = params_mcorr.normcorre_r;
            
            % rename variables to check if normcorre-r has been done
            if isempty(reffile)
                filedir = [ data_locn 'Data/' file(1:8) '/Processed/' file '/mcorr_normcorre-r/' ];
                    if ~exist( filedir, 'dir' ), mkdir( filedir ); end
                fname_tif_gr_mcorr = [filedir file '_2P_XYT_green_mcorr.tif'];
                fname_tif_red_mcorr = [filedir file '_2P_XYT_red_mcorr.tif'];
                fname_mat_mcorr = [filedir file '_mcorr_output.mat'];
                fname_fig = [filedir file '_mcorr_summary.fig'];
            else
                filedir = [ data_locn 'Data/' file(1:8) '/Processed/' file '/mcorr_normcorre-r_ref' reffile '/' ];
                    if ~exist( filedir, 'dir' ), mkdir( filedir ); end
                fname_tif_gr_mcorr = [filedir file '_2P_XYT_green_imreg_ref' reffile '.tif'];
                fname_tif_red_mcorr = [filedir file '_2P_XYT_red_imreg_ref' reffile '.tif'];
                fname_mat_mcorr = [filedir file '_imreg_ref' reffile '_output.mat'];
                fname_fig = [filedir file '_imreg_ref' reffile '_summary.fig'];
            end
            
            if any([ force, ~exist(fname_tif_gr_mcorr,'file'), ~exist(fname_tif_red_mcorr,'file'), ~exist(fname_mat_mcorr,'file') ])
                fprintf( '\tFirst doing rigid correction\n' )
                % first do rigid correction
                if strcmpi(refChannel,'green')
                    [ imG, imR, out_g, out_r, col_shift, shifts, template, ~ ] = normcorre_2ch( imG, imR, params_mcorr.normcorre_r, template );
                else
                    [ imR, imG, out_r, out_g, col_shift, shifts, template, ~ ] = normcorre_2ch( imR, imG, params_mcorr.normcorre_r, template );
                end
                % Save summary figure, tif images, motion correction/registration output matrix
                if force || ~exist(fname_fig,'file'), makeplot(out_g,out_r); end
                saveTifOutputM(out_g, out_r, shifts, col_shift, template, imG, imR);
            else
                % read motion corrected tif files for normcorre-r
                imG = read_file( fname_tif_gr_mcorr );
                imR = read_file( fname_tif_red_mcorr ); 
            end

            % then do non-rigid correction
            fprintf( '\tNow doing non-rigid correction\n' )
            if strcmpi(refChannel,'green')
                [ imG, imR, out_g, out_r, col_shift, shifts, template, ~ ] = normcorre_2ch( imG, imR, params_mcorr.normcorre_nr, template );
            else
                [ imR, imG, out_r, out_g, col_shift, shifts, template, ~ ] = normcorre_2ch( imR, imG, params_mcorr.normcorre_nr, template );
            end
            % Save summary figure, tif images, motion correction/registration output matrix
            if isempty(reffile)
                filedir = [ data_locn 'Data/' file(1:8) '/Processed/' file '/mcorr_normcorre-nr/' ];
                    if ~exist( filedir, 'dir' ), mkdir( filedir ); end
                fname_tif_gr_mcorr = [filedir file '_2P_XYT_green_mcorr.tif'];
                fname_tif_red_mcorr = [filedir file '_2P_XYT_red_mcorr.tif'];
                fname_mat_mcorr = [filedir file '_mcorr_output.mat'];
                fname_fig = [filedir file '_mcorr_summary.fig'];
            else
                filedir = [ data_locn 'Data/' file(1:8) '/Processed/' file '/mcorr_normcorre-nr_ref' reffile '/' ];
                    if ~exist( filedir, 'dir' ), mkdir( filedir ); end
                fname_tif_gr_mcorr = [filedir file '_2P_XYT_green_imreg_ref' reffile '.tif'];
                fname_tif_red_mcorr = [filedir file '_2P_XYT_red_imreg_ref' reffile '.tif'];
                fname_mat_mcorr = [filedir file '_imreg_ref' reffile '_output.mat'];
                fname_fig = [filedir file '_imreg_ref' reffile '_summary.fig'];
            end
            makeplot(out_g,out_r);
            saveTifOutputM(out_g, out_r, shifts, col_shift, template, imG, imR)
        
        elseif strcmpi(mcorr_method,'normcorre-r')
            mcorr_output.params.normcorre_r = params_mcorr.normcorre_r;
            if strcmpi(refChannel,'green')
                [ imG, imR, out_g, out_r, col_shift, shifts, template, ~ ] = normcorre_2ch( imG, imR, params_mcorr.normcorre_r, template );
            else
                [ imR, imG, out_r, out_g, col_shift, shifts, template, ~ ] = normcorre_2ch( imR, imG, params_mcorr.normcorre_r, template );
            end
            % Save summary figure, tif images, motion corrections/registration output matrix
            makeplot(out_g,out_r);
            saveTifOutputM(out_g, out_r, shifts, col_shift, template, imG, imR)
        
        elseif strcmpi(mcorr_method,'normcorre-nr')
            mcorr_output.params.normcorre_nr = params_mcorr.normcorre_nr;
            if strcmpi(refChannel,'green')
                [ imG, imR, out_g, out_r, col_shift, shifts, template, ~ ] = normcorre_2ch( imG, imR, params_mcorr.normcorre_nr, template );
            else
                [ imR, imG, out_g, out_r, col_shift, shifts, template, ~ ] = normcorre_2ch( imR, imG, params_mcorr.normcorre_nr, template );
            end
            % Save summary figure, tif images, motion corrections/registration output matrix
            makeplot(out_g,out_r);
            saveTifOutputM(out_g, out_r, shifts, col_shift, template, imG, imR)
        
        else
            if ~isempty(template)
                str = sprintf('%s: Image registration cannot be done with fftRigid motion correction method.', file);
                cprintf('Errors',str);    
            else
                mcorr_output.params.fftRigid = params_mcorr.fftRigid;
                [ imG, imR, out_g, out_r, col_shift, shifts, template, ~ ] = motionCorrectToNearestPixel( double(imG), double(imR), ...
                                                                                params_mcorr.fftRigid.imscale, params_mcorr.fftRigid.Nimg_ave, ...
                                                                                refChannel, params_mcorr.fftRigid.redoT );
                % Save summary figure, tif images, motion corrections/registration output matrix
                makeplot(out_g,out_r);
                saveTifOutputM(out_g, out_r, shifts, col_shift, template, imG, imR)
            end
        end
    else
        [imG, imR] = load_imagefile( data_locn, file, 0, '_mcorr', params_mcorr );
        mcorr_output = load(fname_mat_mcorr);
        if isfield(mcorr_output,'params') 
            params_mcorr = mcorr_output.params;  % applies to proc data from July 2019
        end
        if isfield(mcorr_output,'options')
            params_mcorr = mcorr_output.options; % applies to proc data prior to July 2019 
        end

        if ~exist_fig
            % If summary fig doesn't exist, create it   
            out_g = mcorr_output.green;
            out_r = mcorr_output.red;
            makeplot(out_g,out_r);
        end
    end
    
    function makeplot(out_g,out_r)
        fh = figure; 
        subplot(221), 
            imagesc( out_g.meanframe ); 
            axis image; colorbar; axis off;
            title( 'Mean frame for raw green' );
        subplot(222), 
            imagesc( out_g.meanregframe ); 
            axis image; colorbar; axis off; 
            title( 'Mean frame for corrected green' );
        subplot(223), 
            imagesc( out_r.meanframe ); 
            axis image; colorbar; axis off; 
            title( 'Mean frame for raw red' );
        subplot(224), 
            imagesc( out_r.meanregframe ); 
            axis image; colorbar; axis off;
            title( 'Mean frame for corrected red' );
        axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0  1],'Box','off',...
            'Visible','off','Units','normalized', 'clipping' , 'off');
            titletext = [file(1:8) '-' file(10:11) '.' file(13:14) '.' file(16:17)];
            text(0.5, 0.98,titletext);
        fname_fig = [filedir file '_mcorr_summary'];
            savefig( fh, fname_fig );
            saveas( fh, fname_fig, 'png' );
        close( fh );
    end

    function saveTifOutputM(out_g, out_r, shifts, col_shift, template, imG, imR)
        % Save output
        mcorr_output.green = out_g;
        mcorr_output.red = out_r;
        mcorr_output.shifts = shifts;
        mcorr_output.col_shift = col_shift;
        mcorr_output.template = template;
        save(fname_mat_mcorr,'-struct','mcorr_output');

        % Save motion corrected tif images
        prevstr = sprintf( '%s: Saving motion corrected tif images...\n', file );
        cprintf('Text',prevstr);
            writeTifStack( imG,fname_tif_gr_mcorr );
            writeTifStack( imR,fname_tif_red_mcorr );
        str = sprintf( '%s: Motion corrected tif images saved\n', file );
        refreshdisp( str, prevstr );
    end
end