% Written by Ann Go
% This script takes as input an image file optionally registered to a local template
% templateloc (template for day/session) and outputs the same file registered 
% to a global template templateglob (template for all sessions).
% INPUTS
%   file            : in the format 'yyyymmdd_HH_MM_SS'
%   templateglob    : in the same format as file. Must be of the same
%       environment as file
%   imregr_params   : struct of fields shifts, options, col_shift. For
%       rigid registration
%   imregnr_params  : struct of same fields as imregr_params. For non-rigid
%       registration
%   templateloc     : (optional) in the same format as file. Must be an image from the same
%       session as file (same day)
%
% OUTPUT
%   file registered to global template templateglob

function imG_globalreg = imreg_global( file, templateglob, imregr_params, imregnr_params, templateloc, force )
    if nargin<5, templateloc = []; end
    if nargin<6, force = false; end
    
    if strcmpi(file, templateglob)
        beep
        cprintf('Text','File to be registered is the same as template. Skipping registration.');    
        return
    end
        
    tic
    
    [data_locn,~,err] = load_neuroSEEmodules;
    if ~isempty(err)
        beep
        cprintf('Errors',err);    
        return
    end
    
    % filenames to save outputs to
    filedir = [ data_locn 'Data/' file(1:8) '/Processed/' file '/imreg_' mcorr_method '_ref' templateglob '/' ];
    fname_tif_gr_mcorr = [filedir file '_2P_XYT_green_imreg_ref' templateglob '.tif'];
    fname_mat_mcorr = [filedir file '_imreg_ref' templateglob '_output.mat'];
    fname_fig = [filedir file '_imreg_summary.fig'];
    
    if any([ force, ~exist(fname_tif_gr_mcorr,'file'), ~exist(fname_mat_mcorr,'file') ])    mcorr_method = 'normcorre';
        % load image file registered to templateloc
        if ~isempty(templateloc)
            [ imG, ~ ] = load_imagefile( data_locn, file, false, '_imreg', mcorr_method, false, templateloc, mcorr_method );
        else
            [ imG, ~ ] = load_imagefile( data_locn, file, false, '_mcorr', mcorr_method, false );
        end
        
        % load templateglob 
        fprintf( '%s: Starting image registration to %s\n', file, templateglob );
        refdir = [data_locn 'Data/' templateglob(1:8) '/Processed/' templateglob '/mcorr_' mcorr_method '/'];
        c = load([refdir templateglob '_mcorr_output.mat']);
        template_g = c.green.meanregframe;
    
        % register to templateglob by applying pre-determined shifts
        % rigid registration
        imG_imregr = apply_shifts( imG, imregr_params.shifts, imregr_params.options, 0, 0, 0, imregr_params.col_shift );

        % non-rigid registration
        imG_globalreg = apply_shifts( imG_imregr, imregnr_params.shifts, imregnr_params.options, 0, 0, 0, imregnr_params.col_shift );

        % summary figure
        out_g.meanframe = mean(imG,3);
        out_g.meanregframe = mean(imG_globalreg,3);
        fh = figure;
        subplot(221), 
            C1 = imfuse( out_g.meanframe, template_g, 'falsecolor', 'Scaling', 'joint', 'ColorChannels', [1 2 0]);
            imshow(C1);  
            title( 'Green: Before registration' );
        subplot(222), 
            C2 = imfuse( out_g.meanregframe,template_g,'falsecolor','Scaling','joint','ColorChannels',[1 2 0] );
            imshow(C2); 
            title( 'Green: After registration' );
        subplot(223), 
            C1 = imfuse( out_g.meanframe, template_g, 'falsecolor', 'Scaling', 'joint', 'ColorChannels', [1 2 0]);
            imshow(C1);  
            title( 'Reversed colours' );
        subplot(224), 
            C2 = imfuse( out_g.meanregframe,template_g,'falsecolor','Scaling','joint','ColorChannels',[1 2 0] );
            imshow(C2); 
            title( 'Reversed colours' );
        if ~exist( filedir, 'dir' ), mkdir( filedir ); end
        savefig( fh, fname_fig );
        saveas( fh, fname_fig, 'png' );
        close( fh );
        
        % save globally registered image
        shifts.r = imregr_params.shifts;
        shifts.nr = imregnr_params.shifts;
        col_shift.r = imregr_params.col_shift;
        col_shift.nr = imregnr_params.col_shift;
        params_mcorr.r = imregr_params.options;
        params_mcorr.nr = imregnr_params.options;
        saveTifOutput( out_g, [], shifts, col_shift, template_g, imG_globalreg, [], [], [], params_mcorr, ...
                    file, fname_mat_mcorr, fname_tif_gr_mcorr, [], templateglob )
    else
        cprintf( '%s already registered to %s. To overwrite existing file, specify force argument as true.\n', file, templateglob );
    end
    
    t = toc;
    str = sprintf('%s: Processing done in %g hrs\n', file, round(t/3600,2));
    cprintf(str)
end