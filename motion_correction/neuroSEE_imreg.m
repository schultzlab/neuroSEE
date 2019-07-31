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
%   params.methods.mcorr_method: ['normcorre' or 'fftRigid'] motion correction method
%                   * CaImAn NoRMCorre method OR fft-rigid method (Katie's)
%   force       : (optional, default: 0) if =1, motion correction will be done even though motion
%                   corrected images already exist
% OUTPUTS
%   imG         : matrix of motion corrected green image stack
%   mcorr_output: cell array containing
%                   green.[ meanframe, meanregframe ]
%                   shifts
%                   template
%   params      : parameters for specific motion correction method

function [ imG, mcorr_output, params ] = neuroSEE_imreg(...
                                                imG, data_locn, file, templatefile, template, params, force )
    
    if nargin<7, force = 0; end
    
    filedir = [data_locn 'Data/' file(1:8) '/Processed/' file '/mcorr_normcorre_ref' templatefile '/'];
    if ~exist( filedir, 'dir' ), mkdir( filedir ); end

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
    
 
    fname_tif_gr = [filedir file '_2P_XYT_green_imreg_ref' templatefile '.tif'];
    fname_mat_mcorr = [filedir file '_imreg_output_ref' templatefile '.mat'];
    fname_fig = [filedir file '_imreg_summary_ref' templatefile '.fig'];
    
    if force || ~exist(fname_tif_gr,'file')
        fprintf( '%s: Starting image registration to %s\n', file, templatefile );
        out_g.meanframe = mean(imG,3);
        [imG, shifts, ~, ~, ~] = normcorre( imG, params.nonrigid, template );
        out_g.meanregframe = mean(imG,3);
        
        % Save summary figure
        makeplot(out_g, template);
        
        % Save output
        mcorr_output.green = out_g;
        mcorr_output.shifts = shifts;
        mcorr_output.template = template;
        mcorr_output.params = params.nonrigid;
        save(fname_mat_mcorr,'-struct','mcorr_output');

        % Save motion corrected tif images
        prevstr = sprintf( '%s: Saving registered tif images...\n', file );
        cprintf('Text',prevstr);
        writeTifStack( imG,fname_tif_gr );
        str = sprintf( '%s: Registered tif images saved\n', file );
        refreshdisp( str, prevstr );
    else
        if ~exist(fname_fig,'file')
            mcorr_output = load(fname_tif_gr);
            out_g = mcorr_output.green;
            stemplate = mcorr_output.template;
            makeplot(out_g, stemplate);
        end
    end
    
    function makeplot(out_g)
        fh = figure; 
        subplot(221), 
            imagesc( out_g.meanframe ); 
            axis image; colorbar; axis off;
            titletext = ['Motion-corrected ' file];
            title(titletext);
        subplot(222), 
            imagesc( template ); 
            axis image; colorbar; axis off; 
            titletext = [templatefile ' (template)'];
            title(titletext);
        subplot(223), 
            C1 = imfuse( out_g.meanframe, template, 'falsecolor', 'Scaling', 'joint', 'ColorChannels', [1 2 0]);
            imshow(C1);  
            axis image; colorbar; axis off; 
            title( 'Before registration' );
        subplot(224), 
            C2 = imfuse(out_g.meanregframe,template,'falsecolor','Scaling','joint','ColorChannels',[1 2 0]);
            imshow(C2); 
            axis image; colorbar; axis off; 
            title( 'After registration' );
        figname = [filedir file '_imreg_summary_ref' templatefile];
            savefig( fh, figname );
            saveas( fh, figname, 'jpg' );
        close( fh );
    end
end