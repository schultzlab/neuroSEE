% Written by Ann Go
%
% This function implements the ABLE fxn (renamed version of neuroSEE_segment
% which Katie adapted from Seig's code) for ROI segmentation when force = 1 or if figure with 
% ROIs don't yet exist in filedir. ROIs are then accepted/rejected on the
% basis of 
%   (1) saturation
%   (2) noise 
%   (3) area
%
% INPUTS
%   imG     : matrix of green image stack
%   imR     : matrix of red image stack
%   filedir : directory of image stacks
%   file    : part of file name of image stacks in the format
%               yyyymmdd_hh_mm_ss
%   params.
%       cellrad     : expected cell radius in pixels
%       maxcells    : expected max number of cells
%       satThresh   : ROI is accepted if its fluorescence is below satThresh
%       satTime     :   satTime fraction of the time
%       areaThresh  : min acceptable ROI area (pixels)
%       noiseThresh : max acceptable ROI noise level
%   force   : if =1, motion correction will be done even though motion
%               corrected images already exist
% OUTPUTS
%   tsG             : time series from green channel
%   tsR             : time series from red channel
%   masks           : ROI masks
%   mean_imratio    : mean ratio image of green and red channels

function [cell_tsG, cell_tsR, masks, mean_imratio, params] = neuroSEE_segment(stack_g, stack_r, mean_r, data_locn, file, params, force)

tic;
    if nargin<7, force = 0;      end

    filedir = fullfile(data_locn,'Data/',file(1:8),'/Processed/',file,'/');
    fname_out = [filedir file '_2P_segment_output.mat'];
    fname_fig = [filedir file '_2P_ROIs.fig'];
    
    % If asked to force overwrite, run ABLE right away
    if force
        cellrad = params.cellrad;
        maxcells = params.maxcells;
        cell_tsR = zeros(1,size(stack_r,3));
        [cell_tsG, masks, mean_imratio] = ABLE_manfredi( stack_g, mean_r, file, cellrad, maxcells );
        
        % [cell_tsG, masks] = refineSegmentation(cell_tsG, masks, params);
        for i = 1:size(masks,3)
            ind = find(masks(:,:,i));
            for j = 1:size(stack_r,3)
                imR_reshaped = reshape(stack_r(:,:,j),size(stack_r,1)*size(stack_r,2),1);
                cell_tsR(i,j) = mean(imR_reshaped(ind));
            end
        end
        
        % Save output
        % Plot masks on summary image and save plot
        plotopts.plot_ids = 1; % set to 1 to view the ID number of the ROIs on the plot
        fig = plotContoursOnSummaryImage(mean_imratio, masks, plotopts);
        savefig(fig,fname_fig);
        saveas(fig,fname_fig,'pdf');
        close(fig);

        % Save masks in a mat file
        save(fname_out,'cell_tsG','cell_tsR','masks','mean_imratio','params')

    else
        yn_ts = exist(fname_out,'file');
        yn_fig = exist(fname_fig,'file');
        
        % If timeseries mat file doesn't exist, run ABLE
        if ~yn_ts
            cellrad = params.cellrad;
            maxcells = params.maxcells;
            cell_tsR = zeros(1,size(stack_r,3));
            save('\Users\mc6017\neuroSEE\08-Jul-2019_video_data.mat');
             [cell_tsG, masks, mean_imratio] = ABLE_manfredi( stack_g, mean_r, file, cellrad, maxcells ); 
             
            % [cell_tsG, masks] = refineSegmentation(cell_tsG, masks, params);
            for i = 1:size(masks,3)
                ind = find(masks(:,:,i));
                for j = 1:size(stack_r,3)
                    imR_reshaped = reshape(stack_r(:,:,j),size(stack_r,1)*size(stack_r,2),1);
                    cell_tsR(i,j) = mean(imR_reshaped(ind));
               
                end
            end
            
    
            
            % Save output
            % Plot masks on summary image and save plot
            plotopts.plot_ids = 1; % set to 1 to view the ID number of the ROIs on the plot
            fig = plotContoursOnSummaryImage(mean_imratio, masks, plotopts);
%             fname_fig = [filedir file '_2P_ROIs'];
%             savefig(fig,fname_fig);
%             saveas(fig,fname_fig,'pdf');
%             close(fig);
     elapsed = toc;
    str = sprintf('MAC -> Neuro see segment finished processing in %2.2f minutes. No errors.',elapsed/60);
    SendSlackNotification('https://hooks.slack.com/services/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',str, '#general','@manfredi.castelli17', [], [], []);
            % Save masks in a mat file
%             save(fname_out,'cell_tsG','cell_tsR','masks','mean_imratio','params')
        else
            % If it exists, load it 
            segmentOutput = load(fname_out);
            cell_tsG = segmentOutput.cell_tsG;
            cell_tsR = segmentOutput.cell_tsR;
            masks = segmentOutput.masks;
            mean_imratio = segmentOutput.mean_imratio;
            params.cellrad = segmentOutput.params.cellrad;
            params.maxcells = segmentOutput.params.maxcells;
            
            % If ROI image doesn't exist, create & save figure
            if ~yn_fig
               plotopts.plot_ids = 1; % set to 1 to view the ID number of the ROIs on the plot
               fig = plotContoursOnSummaryImage(mean_imratio, masks, plotopts);
%                savefig(fig,fname_fig);
%                saveas(fig,fname_fig,'pdf');
               close(fig);
            end
            str = sprintf('%s: Segmentation output loaded\n',file);
            cprintf(str)
        end
    end
end
