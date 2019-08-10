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
    
    if nargin<7, force = 0; end
    
    filedir = [data_locn 'Data/' file(1:8) '/Processed/' file '/mcorr_normcorre_ref' reffile '/'];
    if ~exist( filedir, 'dir' ), mkdir( filedir ); end

    fname_tif_gr = [filedir file '_2P_XYT_green_imreg_ref' reffile '.tif'];
    fname_tif_red = [filedir file '_2P_XYT_red_imreg_ref' reffile '.tif'];
    fname_mat_mcorr = [filedir file '_imreg_ref' reffile '_output.mat'];
    fname_fig_gr = [filedir file '_green_imreg_ref' reffile '_summary.fig'];
    fname_fig_red = [filedir file '_red_imreg_ref' reffile '_summary.fig'];
    
    % template file
    refdir = [data_locn 'Data/' reffile(1:8) '/Processed/' reffile '/mcorr_normcorre/'];
    c = load([refdir reffile '_mcorr_output.mat']);
    template = c.template;
    template_r = c.red.meanregframe;
    clear c        
    
    if any([ force, ~exist(fname_tif_gr,'file'), ~exist(fname_tif_red,'file'), ~exist(fname_mat_mcorr,'file') ])
        fprintf( '%s: Starting image registration to %s\n', file, reffile );
        
        params.nonrigid = NoRMCorreSetParms(...
                            'd1', size(imG,1),...         
                            'd2', size(imG,2),...
                            'grid_size', params.nonrigid.grid_size,...
                            'overlap_pre', params.nonrigid.overlap_pre,...
                            'overlap_post', params.nonrigid.overlap_post,...
                            'iter', params.nonrigid.iter,...
                            'use_parallel', params.nonrigid.use_parallel,...
                            'max_shift', params.nonrigid.max_shift,...
                            'mot_uf', params.nonrigid.mot_uf,...
                            'bin_width', params.nonrigid.bin_width,...
                            'max_dev', params.nonrigid.max_dev,...
                            'us_fac', params.nonrigid.us_fac,...
                            'init_batch', params.nonrigid.init_batch,...
                            'correct_bidir', params.nonrigid.correct_bidir);        

        [ imG, imR, out_g, out_r, col_shift, shifts, ~, ~ ] = normcorre_2ch( imG, imR, params.nonrigid, template );
        
        % Save summary figure
        makeplot(out_g,out_r);
        
        % Save output
        mcorr_output.green = out_g;
        mcorr_output.red = out_r;
        mcorr_output.shifts = shifts;
        mcorr_output.col_shift = col_shift;
        mcorr_output.template = template;
        mcorr_output.params = params.nonrigid;
        save(fname_mat_mcorr,'-struct','mcorr_output');

        % Save motion corrected tif images
        prevstr = sprintf( '%s: Saving registered tif images...\n', file );
        cprintf('Text',prevstr);
        writeTifStack( imG,fname_tif_gr );
        writeTifStack( imR,fname_tif_red );
        str = sprintf( '%s: Registered tif images saved\n', file );
        refreshdisp( str, prevstr );
    else
        mcorr_output = load(fname_mat_mcorr);
        if ~exist(fname_fig_gr,'file') || ~exist(fname_fig_red,'file')
            out_g = mcorr_output.green;
            out_r = mcorr_output.red;
            makeplot(out_g,out_r);
        end
    end
    
    function makeplot(out_g,out_r)
        fh = figure; 
        subplot(221), 
            imagesc( out_g.meanframe ); 
            axis image; axis off; colormap(gray);
            titletext = ['Non-motion-corrected ' file(1:8) '-' file(10:11) '-' file(13:14)];
            title(titletext);
        subplot(222), 
            imagesc( template ); 
            axis image; axis off; colormap(gray); 
            titletext = ['Template: ' reffile(1:8) '-' reffile(10:11) '-' reffile(13:14)];
            title(titletext);
        subplot(223), 
            C1 = imfuse( out_g.meanframe, template, 'falsecolor', 'Scaling', 'joint', 'ColorChannels', [1 2 0]);
            imshow(C1);  
            title( 'Before registration' );
        subplot(224), 
            C2 = imfuse(out_g.meanregframe,template,'falsecolor','Scaling','joint','ColorChannels',[1 2 0]);
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