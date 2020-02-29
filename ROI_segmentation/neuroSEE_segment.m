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

    if nargin<7, list = []; end
    if nargin<6, mean_imR = []; end
    if nargin<5, force = false; end

    mcorr_method = params.methods.mcorr_method;
    segment_method = params.methods.segment_method;
    dofissa = params.methods.dofissa; 
        if dofissa, str_fissa = 'FISSA'; else, str_fissa = 'noFISSA'; end

    if isempty(list)
        filedir = [ data_locn, 'Data/', file(1:8), '/Processed/', file, '/mcorr_', mcorr_method , '/', segment_method, '/' ];
        fname_mat = [filedir file '_segment_output.mat'];
        fname_fig1 = [filedir file '_ROIs.fig'];
        fname_fig2 = [filedir file '_raw_timeseries.fig'];
        fname_fig3 = [filedir file '_df_f.fig'];
    else
        [ mouseid, expname ] = find_mouseIDexpname(list);
        filedir = [ data_locn 'Analysis/' mouseid '/' mouseid '_' expname '/group_proc/' mcorr_method '_' segment_method '_'...
                    str_fissa '/' mouseid '_' expname '_imreg_ref' reffile '/'];
        fname_mat = [filedir mouseid '_' expname '_ref' reffile '_segment_output.mat'];
        fname_fig1 = [filedir mouseid '_' expname '_ref' reffile '_ROIs.fig'];
        fname_fig2 = [filedir mouseid '_' expname '_ref' reffile '_raw_timeseries.fig'];
        fname_fig3 = [filedir mouseid '_' expname '_ref' reffile '_df_f.fig'];
    end

    if force || ~exist(fname_mat,'file')
        maxcells = params.ROIsegment.maxcells;
        cellrad = params.ROIsegment.cellrad;

        if strcmpi(segment_method,'ABLE')
            df_prctile = params.ROIsegment.df_prctile;
            df_medfilt1 = params.ROIsegment.df_medfilt1;
            fr = params.fr;

            [tsG, masks, corr_image] = ABLE_manfredi( imG, mean_imR, maxcells, cellrad );

            % Calculate df_f
            df_f = zeros(size(tsG));
            for i = 1:size(tsG,1)
                x = lowpass( medfilt1(tsG(i,:),df_medfilt1), 1, fr );
                fo = ones(size(x)) * prctile(x,df_prctile);
                while fo == 0
                    fo = ones(size(x)) * prctile(x,df_prctile+5);
                    df_prctile = df_prctile+5;
                end
                df_f(i,:) = (x - fo) ./ fo;
            end

        else
            [df_f, masks, corr_image] = CaImAn( imG, maxcells, cellrad );
            df_f = full(df_f);

            % Extract raw timeseries
            [d1,d2,T] = size(imG);
            tsG = zeros(size(df_f));
            for n = 1:size(masks,3)
                maskind = masks(:,:,n);
                for j = 1:T
                    imG_reshaped = reshape( imG(:,:,j), d1*d2, 1);
                    tsG(n,j) = mean( imG_reshaped(maskind) );
                end
            end
        end

        % Save output
        output.tsG = tsG;
        output.df_f = df_f;
        output.masks = masks;
        output.corr_image = corr_image;
        output.params = params.ROIsegment;
        if ~exist( filedir, 'dir' ), mkdir( filedir ); end
        save(fname_mat,'-struct','output');

        % Plot masks on correlation image and save plot
        makeplot(corr_image, masks);

        fprintf('%s: ROI segmentation done\n',file);

    else
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
        params.ROIsegment = segmentOutput.params;

        % If ROI image doesn't exist, create & save figure
        if any([~exist(fname_fig1,'file'),~exist(fname_fig2,'file'),~exist(fname_fig3,'file')])
           makeplot(corr_image, masks);
        end
        fprintf('%s: Segmentation output found and loaded\n',file);
    end

    function makeplot(corr_image, masks)
        % ROIs overlayed on correlation image
        plotopts.plot_ids = 1; % set to 1 to view the ID number of the ROIs on the plot
        fig = plotContoursOnSummaryImage(corr_image, masks, plotopts);
        savefig(fig, fname_fig1(1:end-4));
        saveas(fig, fname_fig1(1:end-4), 'png');
        close(fig);
        
        % raw timeseries
        fig = figure;
        iosr.figures.multiwaveplot(1:size(tsG,2),1:size(tsG,1),tsG,'gain',5); yticks([]); xticks([]); 
        title('Raw timeseries','Fontweight','normal','Fontsize',12); 
        savefig(fig, fname_fig2(1:end-4));
        saveas(fig, fname_fig2(1:end-4), 'png');
        close(fig);
        
        % dF/F
        fig = figure;
        iosr.figures.multiwaveplot(1:size(df_f,2),1:size(df_f,1),df_f,'gain',5); yticks([]); xticks([]); 
        title('dF/F','Fontweight','normal','Fontsize',12); 
        savefig(fig, fname_fig3(1:end-4));
        saveas(fig, fname_fig3(1:end-4), 'png');
        close(fig);
    end
end
