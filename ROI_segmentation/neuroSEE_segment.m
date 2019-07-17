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


function [tsG, df_f, masks, corr_image, params] = neuroSEE_segment(imG, mean_imR, data_locn, file, mcorr_method, segment_method, params, force)

if nargin<6, force = 0;      end


filedir = [data_locn,'Data/',file(1:8),'/Processed/',file,'/mcorr_',mcorr_method,'/',segment_method,'/'];
if ~exist( filedir, 'dir' ), mkdir( filedir ); end
fname_mat = [filedir file '_segment_output.mat'];
fname_fig = [filedir file '_ROIs.fig'];

% If asked to force overwrite, do roi segmentation right away
if force
    if strcmpi(segment_method,'ABLE')
        [tsG, masks, ~, corr_image] = ABLE( imG, mean_imR, file, params.maxcells, params.cellrad );

        % Calculate df_f
        df_f = zeros(size(tsG));
        for i = 1:size(tsG,1)
            x = lowpass( tsG(i,:), 1, params.imrate );
            fo = ones(size(x)) * prctile(x,5);
            df_f(i,:) = (tsG(i,:) - fo) ./ fo;
            df_f(i,:) = lowpass( medfilt1(df_f(i,:),23), 1, params.imrate );
        end

    else
        [df_f, masks, corr_image] = CaImAn( imG, file, params.maxcells, params.cellrad );

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
    output.params.maxcells = params.maxcells;
    output.params.cellrad = params.cellrad;
    save(fname_mat,'-struct','output');

    % Plot masks on correlation image and save plot
    plotopts.plot_ids = 1; % set to 1 to view the ID number of the ROIs on the plot
    fig = plotContoursOnSummaryImage(corr_image, masks, plotopts);
    savefig(fig,fname_fig);
    saveas(fig,fname_fig,'pdf');
    close(fig);
    

    fprintf('%s: ROI segmentation done\n',file);

else
    yn_mat = exist(fname_mat,'file');
    yn_fig = exist(fname_fig,'file');

    % If timeseries mat file doesn't exist, do roi segmentation
    if ~yn_mat
        if strcmpi(segment_method,'ABLE')
            [tsG, masks, ~, corr_image] = ABLE_manfredi( imG, mean_imR, file, params.maxcells, params.cellrad );

            % Calculate df_f
            df_f = zeros(size(tsG));
            for i = 1:size(tsG,1)
                x = lowpass( tsG(i,:), 1, params.imrate );
                fo = ones(size(x)) * prctile(x,5);
                df_f(i,:) = (tsG(i,:) - fo) ./ fo;
                df_f(i,:) = lowpass( medfilt1(df_f(i,:),23), 1, params.imrate );

            end

        else
            [df_f, masks, corr_image] = CaImAn( imG, file, params.maxcells, params.cellrad );

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
        output.params.maxcells = params.maxcells;
        output.params.cellrad = params.cellrad;
        save(fname_mat,'-struct','output');

        % Plot masks on correlation image and save plot
        plotopts.plot_ids = 1; % set to 1 to view the ID number of the ROIs on the plot
        fig = plotContoursOnSummaryImage(corr_image, masks, plotopts);
        savefig(fig,fname_fig);
        saveas(fig,fname_fig,'pdf');
        close(fig);

        fprintf('%s: ROI segmentation done\n',file);
        

    else
        % If it exists, load it 
        segmentOutput = load(fname_mat);
        if isfield(segmentOutput,'cell_tsG')
            tsG = segmentOutput.cell_tsG;

        else
            tsG = segmentOutput.tsG;
        end
        if isfield(segmentOutput,'df_f')
            df_f = segmentOutput.df_f;
        else
            df_f = zeros(size(tsG));
        end
        masks = segmentOutput.masks;
        if isfield(segmentOutput,'corr_image')
            corr_image = segmentOutput.corr_image;
        else
            corr_image = zeros(size(mean_imR));
        end
        params.cellrad = segmentOutput.params.cellrad;
        params.maxcells = segmentOutput.params.maxcells;

        % If ROI image doesn't exist, create & save figure
        if ~yn_fig
           plotopts.plot_ids = 1; % set to 1 to view the ID number of the ROIs on the plot
           fig = plotContoursOnSummaryImage(corr_image, masks, plotopts);
           savefig(fig,fname_fig);
           saveas(fig,fname_fig,'pdf');
           close(fig);
        end
        str = sprintf('%s: Segmentation output found and loaded\n',file);
        cprintf(str)
    end
end
end

