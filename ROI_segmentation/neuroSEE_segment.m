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
%   data_locn   : GCaMP data repository
%   params      : parameters for specific roi segmentation method
%   fileorlist  : file - part of file name of image stacks in the format
%               yyyymmdd_hh_mm_ss
%               * Can also be a list - name of text file containing filenames of files 
%               to be concatenated for ROI segmentation. Typically in the format 
%               'list_m##_expname.txt'. When a list is given, a reffile
%               should be provided.
%   reffile     : (optional) file to be used as registration template. Required if 
%               fileorlist above is a list.This file is usually part of 'list'
%               but does not have to be. 
%   conc_runs    : (optional) flag if rois were segmented from concatenated files from
%               different environments e.g. fam1fam2fam1-fam1 but rois are
%               for fam1fam2fam1. DO NOT flag for fam1fam2fam1 files since in
%               this case it is understood that the rois are from the
%               concatenated environments.
%   force       : (optional) if =1, roi segmentation will be done even though roi
%               segmentation output already exists
%   mean_imR    : (optional) average of red image stack (only required for ABLE)
%
% OUTPUTS
%   tsG         : raw time series from green channel
%   df_f        : deltaf/f time series from green channel
%   masks       : ROI masks
%   corr_image  : correlation image from green channel
%   params      : parameters for specific roi segmentation method

function [tsG, df_f, masks, corr_image, params] = neuroSEE_segment( imG, data_locn, params, fileorlist, reffile, savedir, force, mean_imR )
                                                                    
    if nargin<8, mean_imR = []; end
    if nargin<7, force = false; end
    if nargin<6, list = []; end
    if nargin<5, reffile = []; end
    
    mcorr_method = params.methods.mcorr_method;
    segment_method = params.methods.segment_method;
    runpatches = params.methods.runpatches;

    roiarea_min = params.ROIsegment.roiarea_min;
    roiarea_max = params.ROIsegment.roiarea_max;
    invcirc_max = params.ROIsegment.invcirc_max;
    overlap_thr = params.ROIsegment.overlap_thr;
        
    % determine whether fileorlist is a file or list
    if ~strncmp(fileorlist,'list',4)
        file = fileorlist; list = [];
    else
        list = fileorlist; file = [];
    end
    if isempty(list)
        savedir = [ data_locn, 'Data/', file(1:8), '/Processed/', file, '/mcorr_', mcorr_method , '/', segment_method, '/' ];
        fname_mat = [savedir file '_segment_output.mat'];
        fname_fig1 = [savedir file '_ROIs.fig'];
        fname_fig2 = [savedir file '_raw_timeseries.fig'];
        fname_fig3 = [savedir file '_df_f.fig'];
        fname_pref = [savedir file];
    else
        [ mouseid, expname ] = find_mouseIDexpname(list);
        fname_mat = [savedir mouseid '_' expname '_ref' reffile '_segment_output.mat'];
        fname_fig1 = [savedir mouseid '_' expname '_ref' reffile '_ROIs.fig'];
        fname_fig2 = [savedir mouseid '_' expname '_ref' reffile '_raw_timeseries.fig'];
        fname_fig3 = [savedir mouseid '_' expname '_ref' reffile '_df_f.fig'];
        fname_pref = [savedir mouseid '_' expname '_ref' reffile ];
    end
    
    prevstr = [];
    if force || ~exist(fname_mat,'file')
        if isempty(list)
            str = sprintf( '%s: Starting ROI segmentation\n', file );
        else
            str = sprintf( '%s: Starting ROI segmentation\n', [mouseid '_' expname] );
        end
        refreshdisp(str, prevstr);
        
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
                [tsG_all, df_f_all, masks_all, corr_image, F0, A] = CaImAn_patches( imG, params.ROIsegment.CaImAn );
            else
                [tsG_all, df_f_all, masks_all, corr_image, F0, A] = CaImAn( imG, params.ROIsegment.CaImAn );
            end
            df_f_all = full(df_f_all);
        end
        
        % Eliminate 
        %   1. rois touching image border 
        %   2. very small and very large rois
        %   3. rois with very high inverse circularity (not round)
        %   4. highly overlapping rois
        
        borderpix = 4;
        area = zeros(size(masks_all,3),1);
        invcirc = zeros(size(masks_all,3),1);
        inc = []; exc = []; 
        try % ROI elimination
            for j = 1:size(masks_all,3)
                mask = masks_all(borderpix:size(masks_all,1)-borderpix,borderpix:size(masks_all,2)-borderpix,j);
                im = imclearborder(mask);
                c = regionprops(im,'area','perimeter');
                try
                    if ~isempty(c)
                        area(j) = c.Area;  % area of each ROI
                        invcirc(j) = (c.Perimeter.^2)/(4*pi*c.Area);
                        if all([area(j)>roiarea_min,...
                                area(j)<roiarea_max,...
                                invcirc(j)<invcirc_max])
                            inc = [inc; j];
                        else
                            exc = [exc; j];
                        end
                    else
                        exc = [exc; j];
                    end
                catch % happens when the ROI is composed of 2 rois
                    exc = [exc; j];
                end
            end
            % eliminate overlapping rois
            exc2 = [];
            for j = 1:length(inc)
                for k = 1:length(inc)
                    if j~=k
                        [~, overlap1, overlap2] = calcROIoverlap(masks_all(:,:,j), masks_all(:,:,k));
                        if overlap1>overlap_thr || overlap2>overlap_thr
                            exc2 = [exc2; j; k];
                            inc(inc == j) = [];
                            inc(inc == k) = [];
                        end
                    end
                end
            end
        
            exc = unique([exc; exc2]);
            inc = setdiff(1:size(masks_all,3),exc);
            masks = masks_all(:,:,inc);     elim_masks = masks_all(:,:,exc);
            tsG = tsG_all(inc,:);           elim_tsG = tsG_all(exc,:);
            df_f = df_f_all(inc,:);         elim_df_f = df_f_all(exc,:);

            % Save output
            if isempty(list)
                str = sprintf( '%s: Saving ROI segmentation output\n', file );
            else
                str = sprintf( '%s: Saving ROI segmentation output\n', [mouseid '_' expname] );
            end
            refreshdisp(str, prevstr);
            prevstr = str;
            
            output.incmasks = inc;  output.excmasks = exc;
            output.tsG = tsG;       output.elim_tsG = elim_tsG;
            output.df_f = df_f;     output.elim_df_f = elim_df_f;
            output.masks = masks;   output.elim_masks = elim_masks;
            output.corr_image = corr_image;
            output.F0 = F0;
            output.A = A;
            output.params = params.ROIsegment;
        catch % Save ROI segmentation output before elimination criteria
            masks = masks_all;     
            tsG = tsG_all;           
            df_f = df_f_all;         

            % Save output
            if isempty(list)
                str = sprintf( '%s: Error in ROI elimination step, skipping step and saving all data\n', file );
            else
                str = sprintf( '%s: Error in ROI elimination step, skipping step and saving all data\n', [mouseid '_' expname] );
            end
            refreshdisp(str, prevstr);
            prevstr = str;
            
            output.tsG = tsG_all;       
            output.df_f = df_f_all;  
            output.masks = masks_all;
            output.corr_image = corr_image;
            output.F0 = F0;
            output.A = A;
            output.params = params.ROIsegment;
        end
        
        if ~exist( savedir, 'dir' ), mkdir( savedir ); fileattrib(savedir,'+w','g','s'); end
        save(fname_mat,'-struct','output');

        % Plot masks on correlation image and save plot
        plotROIsegmentdata(corr_image, masks, elim_masks, tsG, df_f, fname_pref);

        if isempty(list)
            str = sprintf( '%s: ROI segmentation done\n', file );
        else
            str = sprintf( '%s: ROI segmentation done\n', [mouseid '_' expname] );
        end
        refreshdisp(str, prevstr);
    else
        % If it exists, load it
        if isempty(list)
            str = sprintf( '%s: Segmentation output found. Loading...\n', file );
        else
            str = sprintf( '%s: Segmentation output found. Loading...\n', [mouseid '_' expname] );
        end
        refreshdisp(str, prevstr);
        prevstr = str;
        
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
        if isempty(list)
            str = sprintf( '%s: Segmentation output found and loaded\n', file );
        else
            str = sprintf( '%s: Segmentation output found and loaded\n', [mouseid '_' expname] );
        end
        refreshdisp(str, prevstr)
    end
end
