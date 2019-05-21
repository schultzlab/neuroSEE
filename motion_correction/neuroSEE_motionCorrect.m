% Written by Ann Go
%
% This function implements Katie's motionCorrectToNearestPixel fxn if 
% force = 1 or if either green or red motion corrected tif files don't yet
% exist in filedir.
%
% INPUTS
%   imG         : matrix of green image stack
%   imR         : matrix of red image stack
%   data_locn   : GCaMP data repository
%   file        : part of file name of image stacks in the format
%                   yyyymmdd_HH_MM_SS
%   params.Nimg_ave    : number of images to be averaged when correcting for pixel
%                           shift (zippering)
%   params.imscale     : downsampling factor for motion correction
%   force       : if =1, motion correction will be done even though motion
%                   corrected images already exist
% OUTPUTS
%   imG         : matrix of motion corrected green image stack
%   imR         : matrix of motion corrected red image stack
%   mcorr_output: cell array containing
%                   output.green & output.red
%                   each contains
%                       output.<colour>.template : template used for registration
%                       output.<colour>.meanframe : mean frame for original stack 
%                       output.<colour>.meanregframe : mean frame for registered stack
%                       output.<colour>.shift : matrix of x&y shift for each frame
%   params    : params.imscale & params.Nimg_ave

function [imG, imR, mcorr_output, params] = neuroSEE_motionCorrect(imG, imR, data_locn, file, params, force)
    if nargin<6, force = 0;      end
    
    filedir = fullfile(data_locn,'Data/',file(1:8),'/Processed/',file,'/');
        if ~exist(filedir,'dir'), mkdir(filedir); end
    fname_tif_gr_mcorr = [filedir file '_2P_XYT_green_mcorr.tif'];
    fname_tif_red_mcorr = [filedir file '_2P_XYT_red_mcorr.tif'];
%     fname_tif_gr_mcorr = [file '_2P_XYT_green_mcorr.tif'];
%     fname_tif_red_mcorr = [file '_2P_XYT_red_mcorr.tif'];
    fname_mat_mcorr = [filedir file '_2P_mcorr_output.mat'];
    fname_fig = [filedir file '_2P_mcorr_summary.fig'];
        
    % If asked to force overwrite, run motion correction right away
    if force
        imscale = params.imscale;
        Nimg_ave = params.Nimg_ave;
        refChannel = params.refChannel;
        redoT = params.redoT;
        [imG, imR, mcorr_output, fh] = motionCorrectToNearestPixel(double(imG), double(imR), file, imscale, Nimg_ave, refChannel, redoT);
        
        % Save figure
        savefig( fh, fname_fig );
        saveas( fh, fname_fig(1:end-4), 'pdf' );
        close( fh );

        % Save output
        save(fname_mat_mcorr,'-struct','mcorr_output');

        prevstr = sprintf( '%s: Saving motion corrected tif images...\n', file );
        cprintf('Text',prevstr);
            writeTifStack( imG,fname_tif_gr_mcorr );
            writeTifStack( imR,fname_tif_red_mcorr );
        str = sprintf( '%s: Motion corrected tif images saved\n', file );
        refreshdisp( str, prevstr );
    else
        yn_gr_mcorr = exist(fname_tif_gr_mcorr,'file');
        yn_red_mcorr = exist(fname_tif_red_mcorr,'file');
        yn_mat_mcorr = exist(fname_mat_mcorr,'file');
        yn_fig_mcorr = exist(fname_fig,'file');

        % If any of motion corrected tif stacks or motion correction output
        % mat doesn't exist, run motion correction
        if any([~yn_gr_mcorr,~yn_red_mcorr,~yn_mat_mcorr])
            imscale = params.imscale;
            Nimg_ave = params.Nimg_ave;
            refChannel = params.refChannel;
            redoT = params.redoT;
            [imG, imR, mcorr_output, fh] = motionCorrectToNearestPixel(double(imG), double(imR), file, imscale, Nimg_ave, refChannel, redoT);
            
            % Save figure
            savefig( fh, fname_fig );
            saveas( fh, fname_fig(1:end-4), 'pdf' );
            close( fh );

            % Save output
            save(fname_mat_mcorr,'-struct','mcorr_output');

            prevstr = sprintf( '%s: Saving motion corrected tif images...\n', file );
            cprintf('Text',prevstr);
                writeTifStack( imG,fname_tif_gr_mcorr );
                writeTifStack( imR,fname_tif_red_mcorr );
            str = sprintf( '%s: Motion corrected tif images saved\n', file );
            refreshdisp( str, prevstr );
        else
            % If they do exist, load motion corrected tif stacks
            [imG, imR] = load_imagefile(data_locn,file,force,'_mcorr');
            mcorr_output = load(fname_mat_mcorr);
            params.imscale = mcorr_output.params.imscale;
            params.Nimg_ave = mcorr_output.params.Nimg_ave;
            if ~yn_fig_mcorr
                % If summary fig doesn't exist, create it   
                out_g = mcorr_output.green;
                out_r = mcorr_output.red;
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
                fname_fig = [filedir file '_2P_mcorr_summary'];
                    savefig( fh, fname_fig );
                    saveas( fh, fname_fig(1:end-4), 'pdf' );
                close( fh );
            end
        end
    end
