% Written by Ann Go
%
% When force = 1 or if figure with ROIs don't yet exist in filedir, this function
% implements either
% (1) the ABLE fxn (renamed version of neuroSEE_segment
% which Katie adapted from Seig's code) OR
% (2) CaImAn
%
% INPUTS
%   imG         : matrix of green image stack
%   imR         : matrix of red image stack
%   data_locn   : GCaMP data repository
%   file        : part of file name of image stacks in the format
%               yyyymmdd_hh_mm_ss
%   params      : parameters for specific roi segmentation method
%   force       : if =1, roi segmentation will be done even though roi
%               segmentation output already exists
% OUTPUTS
%   tsG         : raw time series from green channel
%   df_f
%   masks       : ROI masks
%   corr_image  : correlation image from green channel
%   params      : parameters for specific roi segmentation method

function [tsG, df_f, masks, corr_image, params] = neuroSEE_segment( imG, data_locn, file, params, force, mean_imR, list, reffile )
    
    if nargin<8, reffile = []; end
    if nargin<7, list = []; end
    if nargin<6, mean_imR = []; end
    if nargin<5, force = false; end

    mcorr_method = params.methods.mcorr_method;
    segment_method = params.methods.segment_method;
    runpatches = params.methods.runpatches;
    dofissa = params.methods.dofissa; 
        if dofissa, str_fissa = 'FISSA'; else, str_fissa = 'noFISSA'; end
    roiarea_thr = params.ROIsegment.roiarea_thr;      % area of roi to be considered a cell
        

    if isempty(list)
        filedir = [ data_locn, 'Data/', file(1:8), '/Processed/', file, '/mcorr_', mcorr_method , '/', segment_method, '/' ];
        fname_mat = [filedir file '_segment_output.mat'];
        fname_fig1 = [filedir file '_ROIs.fig'];
        fname_fig2 = [filedir file '_raw_timeseries.fig'];
        fname_fig3 = [filedir file '_df_f.fig'];
        fname_pref = [filedir file];
    else
        [ mouseid, expname ] = find_mouseIDexpname(list);
        groupreg_method = params.methods.groupreg_method;
        imreg_method = params.methods.imreg_method;
        if strcmpi(imreg_method, mcorr_method)
            filedir = [ data_locn 'Analysis/' mouseid '/' mouseid '_' expname '/group_proc/' groupreg_method '_' imreg_method '_' segment_method '_'...
                        str_fissa '/' mouseid '_' expname '_imreg_ref' reffile '/'];
        else
            filedir = [ data_locn 'Analysis/' mouseid '/' mouseid '_' expname '/group_proc/' groupreg_method '_' imreg_method '_' segment_method '_'...
                        str_fissa '/' mouseid '_' expname '_imreg_ref' reffile '_' mcorr_method '/'];
        end
        fname_mat = [filedir mouseid '_' expname '_ref' reffile '_segment_output.mat'];
        fname_fig1 = [filedir mouseid '_' expname '_ref' reffile '_ROIs.fig'];
        fname_fig2 = [filedir mouseid '_' expname '_ref' reffile '_raw_timeseries.fig'];
        fname_fig3 = [filedir mouseid '_' expname '_ref' reffile '_df_f.fig'];
        fname_pref = [filedir mouseid '_' expname '_ref' reffile ];
    end

    if force || ~exist(fname_mat,'file')
        if isempty(reffile)
            fprintf( '%s: Starting ROI segmentation\n', file );
        else
            fprintf( '%s: Starting ROI segmentation\n', [mouseid '_' expname] );
        end
        if strcmpi(segment_method,'ABLE')
            maxcells = params.ROIsegment.ABLE.maxcells;
            cellrad = params.ROIsegment.ABLE.cellrad;
            df_prctile = params.ROIsegment.ABLE.df_prctile;
            df_medfilt1 = params.ROIsegment.ABLE.df_medfilt1;
            fr = params.ROIsegment.ABLE.fr;

            [tsG_all, masks_all, corr_image] = ABLE_manfredi( imG, mean_imR, maxcells, cellrad );

            % Calculate df_f
            df_f_all = zeros(size(tsG_all));
            for i = 1:size(tsG_all,1)
                x = lowpass( medfilt1(tsG_all(i,:),df_medfilt1), 1, fr );
                fo = ones(size(x)) * prctile(x,df_prctile);
                while fo == 0
                    fo = ones(size(x)) * prctile(x,df_prctile+5);
                    df_prctile = df_prctile+5;
                end
                df_f_all(i,:) = (x - fo) ./ fo;
            end

        else
            if runpatches
                [tsG_all, df_f_all, masks_all, corr_image, F0, GUIdata] = CaImAn_patches( imG, params.ROIsegment.CaImAn );
            else
                [tsG_all, df_f_all, masks_all, corr_image, F0, GUIdata] = CaImAn( imG, params.ROIsegment.CaImAn );
            end
            df_f_all = full(df_f_all);
        end
        
        % Eliminate very small rois and rois touching image border
        area = zeros(size(masks_all,3),1);
        borderpix = 4;
        for j = 1:size(masks_all,3)
            mask = masks_all(borderpix:size(masks_all,1)-borderpix,borderpix:size(masks_all,2)-borderpix,j);
            im = imclearborder(mask);
            c = regionprops(im,'area');
            if ~isempty(c)
                area(j) = c.Area;                    % area of each ROI
            end
        end
        masks = masks_all(:,:,area>roiarea_thr);
        elim_masks = masks_all(:,:,area<roiarea_thr);
        tsG = tsG_all(area>roiarea_thr,:);
        df_f = df_f_all(area>roiarea_thr,:);
        
        % Save output
        output.tsG = tsG;
        output.df_f = df_f;
        output.masks = masks;
        output.elim_masks = elim_masks;
        output.corr_image = corr_image;
        output.F0 = F0;
        output.GUIdata = GUIdata;
        output.params = params.ROIsegment;
        if ~exist( filedir, 'dir' ), mkdir( filedir ); end
        save(fname_mat,'-struct','output');

        % Plot masks on correlation image and save plot
        plotROIsegmentdata(corr_image, masks, elim_masks, tsG, df_f, fname_pref);

        fprintf('%s: ROI segmentation done\n',file);

    else
        if isempty(reffile)
            prevstr = sprintf( '%s: Segmentation output found. Loading...\n', file );
        else
            prevstr = sprintf( '%s: Segmentation output found. Loading...\n', [mouseid '_' expname] );
        end
        cprintf(prevstr)
        
        % If it exists, load it
        segmentOutput = load(fname_mat);
        masks = segmentOutput.masks;
        if isfield(segmentOutput,'cell_tsG')
            tsG = segmentOutput.cell_tsG;
        elseif isfield(segmentOutput,'tsG')
            tsG = segmentOutput.tsG;
        end
        if isfield(segmentOutput,'df_f')
            df_f = segmentOutput.df_f;
        else
            df_f = zeros(size(tsG));
        end
        if isfield(segmentOutput,'corr_image')
            corr_image = segmentOutput.corr_image;
        else
            corr_image = zeros(size(mean_imR));
        end
        if isfield(segmentOutput,'elim_masks')
            elim_masks = segmentOutput.elim_masks;
        else
            elim_masks = [];
        end
        params.ROIsegment = segmentOutput.params;

        % If ROI image doesn't exist, create & save figure
        if any([~exist(fname_fig1,'file'),~exist(fname_fig2,'file'),~exist(fname_fig3,'file')])
           plotROIsegmentdata(corr_image, masks, elim_masks, tsG, df_f, fname_pref);
        end
        if isempty(reffile)
            newstr = sprintf( '%s: Segmentation output found and loaded\n', file );
        else
            newstr = sprintf( '%s: Segmentation output found and loaded\n', [mouseid '_' expname] );
        end
        refreshdisp(newstr, prevstr)
    end
end
