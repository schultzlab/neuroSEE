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

function [ imG, mcorr_output, params_mcorr, imR ] = neuroSEE_motionCorrect2( imG, imR, data_locn, file, mcorr_method, params_mcorr, reffile, force, list, requireRed, mode )    

    if nargin<6, mcorr_method = 'normcorre'; end
    if nargin<7, reffile = []; end
    if nargin<8, force = false; end
    if nargin<9, list = []; end
    if nargin<10, requireRed = true; end
    if nargin<11, mode = 1; end % 1: for motion correction, 2: for image registration
    
    
    refChannel = params_mcorr.refChannel;
    
    if isempty(reffile)
        filedir = [ data_locn 'Data/' file(1:8) '/Processed/' file '/mcorr_' mcorr_method '/' ];
        fname_tif_gr_mcorr = [filedir file '_2P_XYT_green_mcorr.tif'];
        fname_tif_red_mcorr = [filedir file '_2P_XYT_red_mcorr.tif'];
        fname_mat_mcorr = [filedir file '_mcorr_output.mat'];
        fname_fig = [filedir file '_mcorr_summary.fig'];
    else
        if strcmpi(file,reffile)
            beep
            cprintf('Text','File to be registered is the same as template. Skipping registration.');    
            return
        end
        filedir = [ data_locn 'Data/' file(1:8) '/Processed/' file '/imreg_' mcorr_method '_ref' reffile '/' ];
        fname_tif_gr_mcorr = [filedir file '_2P_XYT_green_imreg_ref' reffile '.tif'];
        fname_tif_red_mcorr = [filedir file '_2P_XYT_red_imreg_ref' reffile '.tif'];
        fname_mat_mcorr = [filedir file '_imreg_ref' reffile '_output.mat'];
        fname_fig = [filedir file '_imreg_ref' reffile '_summary.fig'];
    end

    if any([ force,...
             ~exist(fname_tif_gr_mcorr,'file'),...
             and(requireRed, ~exist(fname_tif_red_mcorr,'file')),...
             ~exist(fname_mat_mcorr,'file') ])
        if isempty(reffile)
            str = sprintf( '%s: Doing motion correction\n', file );
            template = [];
            template_g = [];
            template_r = [];
        else
            str = sprintf( '%s: Registering image to %s\n', file, reffile );
            refdir = [data_locn 'Data/' reffile(1:8) '/Processed/' reffile '/mcorr_' mcorr_method '/'];
            c = load([refdir reffile '_mcorr_output.mat']);
            if strcmpi(refChannel,'green')
                template = c.green.meanregframe;
            else
                template = c.red.meanregframe;
            end
            template_g = c.green.meanregframe;
            template_r = c.red.meanregframe;
        end
        cprintf( 'Text', str );
        
        if strcmpi(mcorr_method,'normcorre')     
            % rename variables to check if normcorre-r has been done
            if isempty(reffile)
                filedir = [ data_locn 'Data/' file(1:8) '/Processed/' file '/mcorr_normcorre-r/' ];
                fname_tif_gr_mcorr = [filedir file '_2P_XYT_green_mcorr.tif'];
                fname_tif_red_mcorr = [filedir file '_2P_XYT_red_mcorr.tif'];
                fname_mat_mcorr = [filedir file '_mcorr_output.mat'];
                fname_fig = [filedir file '_mcorr_summary.fig'];
            else
                filedir = [ data_locn 'Data/' file(1:8) '/Processed/' file '/imreg_normcorre-r_ref' reffile '/' ];
                fname_tif_gr_mcorr = [filedir file '_2P_XYT_green_imreg_ref' reffile '.tif'];
                fname_tif_red_mcorr = [filedir file '_2P_XYT_red_imreg_ref' reffile '.tif'];
                fname_mat_mcorr = [filedir file '_imreg_ref' reffile '_output.mat'];
                fname_fig = [filedir file '_imreg_ref' reffile '_summary.fig'];
            end
            
            if any([ force,...
                     ~exist(fname_tif_gr_mcorr,'file'),...
                     and(requireRed, ~exist(fname_tif_red_mcorr,'file')),...
                     ~exist(fname_mat_mcorr,'file') ])
                if isempty(reffile)
                    fprintf( '\tFirst doing rigid correction\n' );
                else
                    fprintf( '\tFirst doing rigid registration\n' );
                end
                % first do rigid correction
                if strcmpi(refChannel,'green')
                    if mode == 1
                        [ imG, imR, out_g, out_r, col_shift, shifts, template, ~ ] = normcorre_2ch( imG, imR, params_mcorr.normcorre_r, template );
                    else
                        out_g.meanframe = mean(imG,3); out_r.meanframe = mean(imR,3);
                        [~, shifts, ~, params_mcorr.normcorre_r, col_shift] = normcorre(mean(imG,3), params_mcorr.normcorre_r, template);
                        A = struct('shifts', repmat(shifts, size(imG,3), 1));
                        imG = apply_shifts( imG, A.shifts, params_mcorr.normcorre_r, 0, 0, 0, col_shift );
                        imR = apply_shifts( imR, A.shifts, params_mcorr.normcorre_r, 0, 0, 0, col_shift );
                        out_g.meanregframe = mean(imG,3); out_r.meanregframe = mean(imR,3);
                    end
                else
                    if mode == 1
                        [ imR, imG, out_r, out_g, col_shift, shifts, template, ~ ] = normcorre_2ch( imR, imG, params_mcorr.normcorre_r, template );
                    else
                        out_g.meanframe = mean(imG,3); out_r.meanframe = mean(imR,3);
                        [~, shifts, ~, params_mcorr.normcorre_r, col_shift] = normcorre(mean(imR,3), params_mcorr.normcorre_r, template);
                        A = struct('shifts', repmat(shifts, size(imR,3), 1));
                        imG = apply_shifts( imG, A.shifts, params_mcorr.normcorre_r, 0, 0, 0, col_shift );
                        imR = apply_shifts( imR, A.shifts, params_mcorr.normcorre_r, 0, 0, 0, col_shift );
                        out_g.meanregframe = mean(imG,3); out_r.meanregframe = mean(imR,3);
                    end
                end
                % Save summary figure, tif images, motion correction/registration output matrix
                if force || ~exist(fname_fig,'file'), makeplot(out_g,out_r); end
                saveTifOutputM(out_g, out_r, shifts, col_shift, template, imG, imR, template_g, template_r, params_mcorr.normcorre_r, file);
            else
                % read motion corrected tif files for normcorre-r
                imG = read_file( fname_tif_gr_mcorr );
                imR = read_file( fname_tif_red_mcorr ); 
            end

            % then do non-rigid correction
            if isempty(reffile)
                fprintf( '\tNow doing non-rigid correction\n' );
            else
                fprintf( '\tNow doing non-rigid registration\n' );
            end
            if strcmpi(refChannel,'green')
                if mode == 1
                    [ imG, imR, out_g, out_r, col_shift, shifts, template, ~ ] = normcorre_2ch( imG, imR, params_mcorr.normcorre_nr, template );
                else
                    out_g.meanframe = mean(imG,3); out_r.meanframe = mean(imR,3);
                    [~, shifts, ~, params_mcorr.normcorre_nr, col_shift] = normcorre(mean(imG,3), params_mcorr.normcorre_nr, template);
                    A = struct('shifts', repmat(shifts, size(imG,3), 1));
                    imG = apply_shifts( imG, A.shifts, params_mcorr.normcorre_nr, 0, 0, 0, col_shift );
                    imR = apply_shifts( imR, A.shifts, params_mcorr.normcorre_nr, 0, 0, 0, col_shift );
                    out_g.meanregframe = mean(imG,3); out_r.meanregframe = mean(imR,3);
                end
            else
                if mode == 1
                    [ imR, imG, out_r, out_g, col_shift, shifts, template, ~ ] = normcorre_2ch( imR, imG, params_mcorr.normcorre_nr, template );
                else
                    out_g.meanframe = mean(imG,3); out_r.meanframe = mean(imR,3);
                    [~, shifts, ~, params_mcorr.normcorre_nr, col_shift] = normcorre(mean(imR,3), params_mcorr.normcorre_nr, template);
                    A = struct('shifts', repmat(shifts, size(imR,3), 1));
                    imG = apply_shifts( imG, A.shifts, params_mcorr.normcorre_nr, 0, 0, 0, col_shift );
                    imR = apply_shifts( imR, A.shifts, params_mcorr.normcorre_nr, 0, 0, 0, col_shift );
                    out_g.meanregframe = mean(imG,3); out_r.meanregframe = mean(imR,3);
                end
            end
            % Save summary figure, tif images, motion correction/registration output matrix
            if isempty(reffile)
                filedir = [ data_locn 'Data/' file(1:8) '/Processed/' file '/mcorr_normcorre/' ];
                fname_tif_gr_mcorr = [filedir file '_2P_XYT_green_mcorr.tif'];
                fname_tif_red_mcorr = [filedir file '_2P_XYT_red_mcorr.tif'];
                fname_mat_mcorr = [filedir file '_mcorr_output.mat'];
                fname_fig = [filedir file '_mcorr_summary.fig'];
            else
                filedir = [ data_locn 'Data/' file(1:8) '/Processed/' file '/imreg_normcorre_ref' reffile '/' ];
                fname_tif_gr_mcorr = [filedir file '_2P_XYT_green_imreg_ref' reffile '.tif'];
                fname_tif_red_mcorr = [filedir file '_2P_XYT_red_imreg_ref' reffile '.tif'];
                fname_mat_mcorr = [filedir file '_imreg_ref' reffile '_output.mat'];
                fname_fig = [filedir file '_imreg_ref' reffile '_summary.fig'];
            end
            makeplot(out_g,out_r);
            saveTifOutputM(out_g, out_r, shifts, col_shift, template, imG, imR, template_g, template_r, params_mcorr, file);
        
        elseif  strcmpi(mcorr_method,'normcorre-r')     
            mcorr_output.params.normcorre_r = params_mcorr.normcorre_r;
            if strcmpi(refChannel,'green')
                if mode == 1
                    [ imG, imR, out_g, out_r, col_shift, shifts, template, ~ ] = normcorre_2ch( imG, imR, params_mcorr.normcorre_r, template );
                else
                    out_g.meanframe = mean(imG,3); out_r.meanframe = mean(imR,3);
                    [~, shifts, ~, params_mcorr.normcorre_r, col_shift] = normcorre(mean(imG,3), params_mcorr.normcorre_r, template);
                    A = struct('shifts', repmat(shifts, size(imG,3), 1));
                    imG = apply_shifts( imG, A.shifts, params_mcorr.normcorre_r, 0, 0, 0, col_shift );
                    imR = apply_shifts( imR, A.shifts, params_mcorr.normcorre_r, 0, 0, 0, col_shift );
                    out_g.meanregframe = mean(imG,3); out_r.meanregframe = mean(imR,3);
                end
            else
                if mode == 1
                    [ imR, imG, out_r, out_g, col_shift, shifts, template, ~ ] = normcorre_2ch( imR, imG, params_mcorr.normcorre_r, template );
                else
                    out_g.meanframe = mean(imG,3); out_r.meanframe = mean(imR,3);
                    [~, shifts, ~, params_mcorr.normcorre_r, col_shift] = normcorre(mean(imR,3), params_mcorr.normcorre_r, template);
                    A = struct('shifts', repmat(shifts, size(imR,3), 1));
                    imG = apply_shifts( imG, A.shifts, params_mcorr.normcorre_r, 0, 0, 0, col_shift );
                    imR = apply_shifts( imR, A.shifts, params_mcorr.normcorre_r, 0, 0, 0, col_shift );
                    out_g.meanregframe = mean(imG,3); out_r.meanregframe = mean(imR,3);
                end
            end
            % Save summary figure, tif images, motion corrections/registration output matrix
            makeplot(out_g,out_r);
            saveTifOutputM(out_g, out_r, shifts, col_shift, template, imG, imR, template_g, template_r, params_mcorr.normcorre_r, file);
        
        elseif  strcmpi(mcorr_method,'normcorre-nr')
            mcorr_output.params.normcorre_nr = params_mcorr.normcorre_nr;
            if strcmpi(refChannel,'green')
                if mode == 1
                    [ imG, imR, out_g, out_r, col_shift, shifts, template, ~ ] = normcorre_2ch( imG, imR, params_mcorr.normcorre_nr, template );
                else
                    out_g.meanframe = mean(imG,3); out_r.meanframe = mean(imR,3);
                    [~, shifts, ~, params_mcorr.normcorre_nr, col_shift] = normcorre(mean(imG,3), params_mcorr.normcorre_nr, template);
                    A = struct('shifts', repmat(shifts, size(imG,3), 1));
                    imG = apply_shifts( imG, A.shifts, params_mcorr.normcorre_nr, 0, 0, 0, col_shift );
                    imR = apply_shifts( imR, A.shifts, params_mcorr.normcorre_nr, 0, 0, 0, col_shift );
                    out_g.meanregframe = mean(imG,3); out_r.meanregframe = mean(imR,3);
                end
            else
                if mode == 1
                    [ imR, imG, out_r, out_g, col_shift, shifts, template, ~ ] = normcorre_2ch( imR, imG, params_mcorr.normcorre_nr, template );
                else
                    out_g.meanframe = mean(imG,3); out_r.meanframe = mean(imR,3);
                    [~, shifts, ~, params_mcorr.normcorre_nr, col_shift] = normcorre(mean(imR,3), params_mcorr.normcorre_nr, template);
                    A = struct('shifts', repmat(shifts, size(imR,3), 1));
                    imG = apply_shifts( imG, A.shifts, params_mcorr.normcorre_nr, 0, 0, 0, col_shift );
                    imR = apply_shifts( imR, A.shifts, params_mcorr.normcorre_nr, 0, 0, 0, col_shift );
                    out_g.meanregframe = mean(imG,3); out_r.meanregframe = mean(imR,3);
                end
            end
            % Save summary figure, tif images, motion corrections/registration output matrix
            makeplot(out_g,out_r);
            saveTifOutputM(out_g, out_r, shifts, col_shift, template, imG, imR, template_g, template_r, params_mcorr.normcorre_nr, file);
        
        else
            if ~isempty(template)
                str = sprintf('%s: Cannot do image registration with fftRigid method. Choose another method.', file);
                cprintf('Errors',str);    
            else
                mcorr_output.params.fftRigid = params_mcorr.fftRigid;
                [ imG, imR, out_g, out_r, col_shift, shifts, template, ~ ] = motionCorrectToNearestPixel( double(imG), double(imR), ...
                                                                                params_mcorr.fftRigid.imscale, params_mcorr.fftRigid.Nimg_ave, ...
                                                                                refChannel, params_mcorr.fftRigid.redoT );
                % Save summary figure, tif images, motion corrections/registration output matrix
                makeplot(out_g,out_r);
                saveTifOutputM(out_g, out_r, shifts, col_shift, template, imG, imR, template_g, template_r, params_mcorr.fftRigid, file);
            end
        end
        
        if isempty(reffile)
            fprintf( '%s: Motion correction done\n', file );
        else
            fprintf( '%s: Image registration to %s done\n', file, reffile );
        end
    else
        if nargout>3, load_imR = true; else, load_imR = false; end
        if isempty(reffile)
            [imG, imR] = load_imagefile( data_locn, file, false, '_mcorr', mcorr_method, load_imR );
        else
            if ~isempty(list)
                [mouseid,expname] = find_mouseIDexpname( list );
                str = sprintf('%s: Loading registered images...\n', [mouseid '_' expname '_' file]);
            else
                str = sprintf('%s: Loading registered images...\n', [file '_ref' reffile]);
            end
            cprintf( 'Text', str );
            
            imG = read_file(fname_tif_gr_mcorr);
            if load_imR
                imR = read_file(fname_tif_red_mcorr);
            else
                imR = [];
            end
            
            if ~isempty(list)
                newstr = sprintf('%s: Registered images loaded\n', [mouseid '_' expname '_' file]);
            else
                newstr = sprintf('%s: Registered images loaded\n', [file '_ref' reffile]);
            end
            refreshdisp( newstr, str )
        end
        
        mcorr_output = load(fname_mat_mcorr);
        if isfield(mcorr_output,'params') 
            params_mcorr = mcorr_output.params;  % applies to proc data from July 2019
        end
        if isfield(mcorr_output,'options')
            params_mcorr = mcorr_output.options; % applies to proc data prior to July 2019 
        end

        if ~exist(fname_fig,'file')
            % If summary fig doesn't exist, create it
            if ~isempty(reffile)
                refdir = [data_locn 'Data/' reffile(1:8) '/Processed/' reffile '/mcorr_' mcorr_method '/'];
                c = load([refdir reffile '_mcorr_output.mat']);
                template_g = c.green.meanregframe;
                template_r = c.red.meanregframe;
            end
            out_g = mcorr_output.green;
            out_r = mcorr_output.red;
            makeplot(out_g,out_r);
        end
    end
    
    function makeplot(out_g,out_r)
        fh = figure;
        if isempty(reffile)
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
        else
            subplot(221), 
                C1 = imfuse( out_g.meanframe, template_g, 'falsecolor', 'Scaling', 'joint', 'ColorChannels', [1 2 0]);
                imshow(C1);  
                title( 'Green: Before registration' );
            subplot(222), 
                C2 = imfuse(out_g.meanregframe,template_g,'falsecolor','Scaling','joint','ColorChannels',[1 2 0]);
                imshow(C2); 
                title( 'Green: After registration' );
            subplot(223), 
                C1 = imfuse( out_r.meanframe, template_r, 'falsecolor', 'Scaling', 'joint', 'ColorChannels', [1 2 0]);
                imshow(C1);  
                title( 'Red: Before registration' );
            subplot(224), 
                C2 = imfuse(out_r.meanregframe,template_r,'falsecolor','Scaling','joint','ColorChannels',[1 2 0]);
                imshow(C2); 
                title( 'Red: After registration' );
        end
        if isempty(reffile)
            fname_fig = [filedir file '_mcorr_summary'];
        else
            fname_fig = [filedir file '_imreg_ref' reffile '_summary'];
        end
            if ~exist( filedir, 'dir' ), mkdir( filedir ); fileattrib filedir +w '' s; end
            savefig( fh, fname_fig );
            saveas( fh, fname_fig, 'png' );
        close( fh );
    end

    function saveTifOutputM(out_g, out_r, shifts, col_shift, template, imG, imR, template_g, template_r, params_mcorr, file)
        % Save output
        mcorr_output.green = out_g;
        if ~isempty(out_r)
            mcorr_output.red = out_r;
        end
        mcorr_output.shifts = shifts;
        mcorr_output.col_shift = col_shift;
        mcorr_output.template = template;
        mcorr_output.params = params_mcorr;
        if ~isempty(template_g), mcorr_output.template_g = template_g; end
        if ~isempty(template_r), mcorr_output.template_r = template_r; end
        save(fname_mat_mcorr,'-struct','mcorr_output');

        % Save motion corrected tif images
        if isempty(reffile)
            prevstr = sprintf( '%s: Saving motion corrected tif images...\n', file );
        else
            prevstr = sprintf( '%s: Saving registered tif images...\n', file );
        end
        cprintf('Text',prevstr);
            writeTifStack( imG,fname_tif_gr_mcorr );
            if ~isempty(imR), writeTifStack( imR,fname_tif_red_mcorr ); end
        if isempty(reffile)
            str = sprintf( '%s: Motion corrected tif images saved\n', file );
        else
            str = sprintf( '%s: Registered tif images saved\n', file );
        end
        refreshdisp( str, prevstr );
    end
end