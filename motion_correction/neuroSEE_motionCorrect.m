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
%   params.methods.mcorr_method: ['normcorre' or 'fftRigid'] motion correction method
%                   * CaImAn NoRMCorre method OR fft-rigid method (Katie's)
%   force       : (optional, default: 0) if =1, motion correction will be done even though motion
%                   corrected images already exist
% OUTPUTS
%   imG         : matrix of motion corrected green image stack
%   imR         : matrix of motion corrected red image stack
%   mcorr_output: cell array containing
%                   green.[ meanframe, meanregframe ]
%                   red.[ meanframe, meanregframe ]
%                   shifts.[ zipper_shift, shifts ]
%                   template
%   params      : parameters for specific motion correction method

function [ imG, imR, mcorr_output, params ] = neuroSEE_motionCorrect(...
                                                imG, imR, data_locn, file, params, force )
    
    if nargin<6, force = 0;      end
    mcorr_method = params.methods.mcorr_method;
    
    if strcmpi(mcorr_method,'normcorre')
        filedir = fullfile( data_locn, 'Data/', file(1:8), '/Processed/', file, '/mcorr_normcorre/' );
        if ~isempty(imG)
            params.nonrigid = NoRMCorreSetParms(...
                'd1',size(imG,1),...         
                'd2',size(imG,2));    
        end
    else
        filedir = fullfile( data_locn, 'Data/', file(1:8), '/Processed/', file, '/mcorr_fftRigid/' );
    end
    if ~exist( filedir, 'dir' ), mkdir( filedir ); end
 
    fname_tif_gr_mcorr = [filedir file '_2P_XYT_green_mcorr.tif'];
    fname_tif_red_mcorr = [filedir file '_2P_XYT_red_mcorr.tif'];
    fname_mat_mcorr = [filedir file '_mcorr_output.mat'];
    fname_fig = [filedir file '_mcorr_summary.fig'];
    
    exist_gr = exist(fname_tif_gr_mcorr,'file');
    exist_red = exist(fname_tif_red_mcorr,'file');
    exist_mat = exist(fname_mat_mcorr,'file');
    exist_fig = exist(fname_fig,'file');
        
    if any([ force, ~exist_gr, ~exist_red, ~exist_mat ])
        str = sprintf( '%s: Starting motion correction\n', file );
        cprintf( 'Text', str );
        if strcmpi(mcorr_method,'normcorre')
            mcorr_output.params = params.nonrigid;
            [ imG, imR, out_g, out_r, col_shift, shifts, template, ~ ] = normcorre_2ch( imG, imR, params.nonrigid );
        else
            mcorr_output.params = params.fftRigid;
            imscale = params.fftRigid.imscale;
            Nimg_ave = params.fftRigid.Nimg_ave;
            refChannel = params.fftRigid.refChannel;
            redoT = params.fftRigid.redoT;
            [ imG, imR, out_g, out_r, col_shift, shifts, template, ~ ] = motionCorrectToNearestPixel( double(imG), double(imR), file, ...
                                                                            imscale, Nimg_ave, refChannel, redoT );
        end
        
        % Save summary figure
        makeplot(out_g,out_r);
        
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
    else
        [imG, imR] = load_imagefile( data_locn, file, 0, '_mcorr', params );
        mcorr_output = load(fname_mat_mcorr);
        if strcmpi(mcorr_method,'normcorre')
            if isfield(mcorr_output,'params') 
                params.nonrigid = mcorr_output.params;  % applies to proc data from July 2019
            end
            if isfield(mcorr_output,'options')
                params.nonrigid = mcorr_output.options; % applies to proc data prior to July 2019 
            end
        else
            params.fftRigid = mcorr_output.params;            
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
            saveas( fh, fname_fig, 'jpg' );
        close( fh );
    end
end