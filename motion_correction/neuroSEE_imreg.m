% Written by Ann Go
%
% This function registers imG to the template using normcorre.
%
% INPUTS
%   imG         : matrix of green image stack
%   data_locn   : GCaMP data repository
%   file        : part of file name of image stacks in the format
%                   yyyymmdd_HH_MM_SS
%   templatefile : filename of reference template
%   template    : matrix of reference template
%   params
%   force       : (optional, default: 0) if =1, motion correction will be done even though motion
%                   corrected images already exist
% OUTPUTS
%   imG         : matrix of motion corrected green image stack
%   mcorr_output: cell array containing
%                   green.[ meanframe, meanregframe ]
%                   shifts
%                   template
%   params      : parameters for specific motion correction method

function [ imG, imR, mcorr_output, params ] = neuroSEE_imreg(...
                                                imG, imR, data_locn, file, reffile, params, force )
                                            
    if nargin<7, force = false; end
    if isfield(params.mcorr,'refChannel')
        refChannel = params.mcorr.refChannel;
    else
        refChannel = 'green';
    end
    mcorr_method = params.methods.mcorr_method;
   
    % filenames
    filedir = [data_locn 'Data/' file(1:8) '/Processed/' file '/mcorr_' mcorr_method '_ref' reffile '/'];
    if ~exist( filedir, 'dir' ), mkdir( filedir ); end

    fname_tif_gr = [filedir file '_2P_XYT_green_imreg_ref' reffile '.tif'];
    fname_tif_red = [filedir file '_2P_XYT_red_imreg_ref' reffile '.tif'];
    fname_mat_mcorr = [filedir file '_imreg_ref' reffile '_output.mat'];
    
    if any([ force, ~exist(fname_tif_gr,'file'), ~exist(fname_tif_red,'file'), ~exist(fname_mat_mcorr,'file') ])
        % template file
        refdir = [data_locn 'Data/' reffile(1:8) '/Processed/' reffile '/mcorr_' mcorr_method '/'];
        c = load([refdir reffile '_mcorr_output.mat']);
        template_g = c.green.meanregframe;
        template_r = c.red.meanregframe;

        % image registration
        fprintf( '%s: Starting image registration to %s\n', file, reffile );
        if strcmpi(refChannel,'green')
            if strcmpi(mcorr_method,'normcorre-r')
                [ imG, imR, out_g, out_r, col_shift, shifts, ~, ~ ] = normcorre_2ch( imG, imR, params.mcorr.normcorre_r, template_g );
            elseif strcmpi(mcorr_method,'normcorre-nr')
                [ imG, imR, out_g, out_r, col_shift, shifts, ~, ~ ] = normcorre_2ch( imG, imR, params.mcorr.normcorre_nr, template_g );
            elseif strcmpi(mcorr_method,'normcorre')
                [ imG, imR, ~, ~, ~, ~, ~, ~ ] = normcorre_2ch( imG, imR, params.mcorr.normcorre_r, template_g );
                [ imG, imR, out_g, out_r, col_shift, shifts, ~, ~ ] = normcorre_2ch( imG, imR, params.mcorr.normcorre_nr, template_g );
            end
        else
            if strcmpi(mcorr_method,'normcorre-r')
                [ imR, imG, out_r, out_g, col_shift, shifts, ~, ~ ] = normcorre_2ch( imR, imG, params.mcorr.normcorre_r, template_r );
            elseif strcmpi(mcorr_method,'normcorre-nr')
                [ imR, imG, out_r, out_g, col_shift, shifts, ~, ~ ] = normcorre_2ch( imR, imG, params.mcorr.normcorre_nr, template_r );
            elseif strcmpi(mcorr_method,'normcorre')
                [ imR, imG, ~, ~, ~, ~, ~, ~ ] = normcorre_2ch( imR, imG, params.mcorr.normcorre_r, template_r );
                [ imR, imG, out_r, out_g, col_shift, shifts, ~, ~ ] = normcorre_2ch( imR, imG, params.mcorr.normcorre_nr, template_r );
            end
        end

        % Save summary figure
        makeplot(out_g,out_r);

        % Save output
        mcorr_output.green = out_g;
        mcorr_output.red = out_r;
        mcorr_output.shifts = shifts;
        mcorr_output.col_shift = col_shift;
        mcorr_output.template_g = template_g;
        mcorr_output.template_r = template_r;
        mcorr_output.params = params.mcorr;
        save(fname_mat_mcorr,'-struct','mcorr_output');

        % Save motion corrected tif images
        prevstr = sprintf( '%s: Saving registered tif images...\n', file );
        cprintf('Text',prevstr);
        writeTifStack( imG,fname_tif_gr );
        writeTifStack( imR,fname_tif_red );
        str = sprintf( '%s: Registered tif images saved\n', file );
        refreshdisp( str, prevstr );
    else
        fprintf( '%s: Registered image found. Skipping registration\n', file, reffile );
    end

    function makeplot(out_g,out_r)
        fh = figure; 
        subplot(221), 
            imagesc( out_g.meanframe ); 
            axis image; axis off; colormap(gray);
            titletext = ['Non-motion-corrected ' file(1:8) '-' file(10:11) '-' file(13:14)];
            title(titletext);
        subplot(222), 
            imagesc( template_g ); 
            axis image; axis off; colormap(gray); 
            titletext = ['Template: ' reffile(1:8) '-' reffile(10:11) '-' reffile(13:14)];
            title(titletext);
        subplot(223), 
            C1 = imfuse( out_g.meanframe, template_g, 'falsecolor', 'Scaling', 'joint', 'ColorChannels', [1 2 0]);
            imshow(C1);  
            title( 'Before registration' );
        subplot(224), 
            C2 = imfuse(out_g.meanregframe,template_g,'falsecolor','Scaling','joint','ColorChannels',[1 2 0]);
            imshow(C2); 
            title( 'After registration' );
        figname = [filedir file '_green_imreg_ref' reffile '_summary'];
            savefig( fh, figname );
            saveas( fh, figname, 'png' );
        close( fh );

        fh = figure; 
        subplot(221), 
            imagesc( out_r.meanframe ); 
            axis image; axis off; colormap(gray);
            titletext = ['Non-motion-corrected ' file(1:8) '-' file(10:11) '-' file(13:14)];
            title(titletext);
        subplot(222), 
            imagesc( template_r ); 
            axis image; axis off; colormap(gray); 
            titletext = ['Template: ' reffile(1:8) '-' reffile(10:11) '-' reffile(13:14)];
            title(titletext);
        subplot(223), 
            C1 = imfuse( out_r.meanframe, template_r, 'falsecolor', 'Scaling', 'joint', 'ColorChannels', [1 2 0]);
            imshow(C1);  
            title( 'Before registration' );
        subplot(224), 
            C2 = imfuse(out_r.meanregframe,template_r,'falsecolor','Scaling','joint','ColorChannels',[1 2 0]);
            imshow(C2); 
            title( 'After registration' );
        figname = [filedir file '_red_imreg_ref' reffile '_summary'];
            savefig( fh, figname );
            saveas( fh, figname, 'png' );
        close( fh );

    end
end
