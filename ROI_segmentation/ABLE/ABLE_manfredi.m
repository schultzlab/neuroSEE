%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   NeuroSEE: An automated Neuronal Source Extractionv
%             and Exploration toolbox
%
%   Author: Manfredi Castelli - Seigfred Prado
%   Supervisor: Simon Schultz
%   Acknowledgment: Stephanie Reynolds, Pier Luigi Dragotti
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Hybrid Level Set Segmentation
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   INPUTS:
%       stack_g:    green channel matrix
%       mean_r:     mean image of red channel matrix
%       cellRadius: expected radius of the cell in pixels
%       maxCells:   estimated maximum number of cells in a frame
%                       yyyymmdd_hh_mm_ss
%   OUTPUTS:
%       cell_tsG:   green channel time series
%       cell_tsR:   red channel time series
%       masks:      masks of the segmented ROIs
%       fig:        fig handle of mean image with ROIs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




function [tsG, masks, overl_corr] = ABLE_manfredi( stack_g, mean_r, file, cellrad, maxcells )
   tic;
   str = sprintf('%s: Extracting ROIs with ABLE\n', file);
   cprintf(str);

%% Initialise the summary images
   szsmooth = [3 3];   % higher values could end up into merging differnet cells into once

   % Compute the summary images, i.e. mean and correlation images
   mean_g = mean(stack_g,3);
   corr_g = crossCorr(stack_g(:,:,:)); %1:2:end));
   mean_imratio = mean_g ./ mean_r;

   % Set the initial images
   init_g = medfilt2( corr_g, szsmooth );
   red = mat2gray(mean_r);
   red = red - 0.3;
   red = medfilt2(red,szsmooth);

   corr = imadjust(init_g);

  % Setting Overlay to show differences between red and green channel images
   overl_corr = imfuse(imadjust(red),imadjust(corr),'falsecolor','Scaling','joint','ColorChannels',[1 2 0]);

   R0     = mean_imratio - min(mean_imratio(:));
   R0     = R0 / max(R0(:));

   % Set initial parameters
   if nargin<4
       cellrad  = 10; % cell radius
   end
   if nargin<5
       maxcells = 200;
   end

   alpha  = 0.05; % starts with a low value then increases during tuning
   init_opt.blur_radius = 4; % default is 1; radius of the blurring applied
                             % to the input summary image

   % Secondary Metric
   % - divide the red over the green average image so that the only remaining
   %   signal is from neuronal nuclei, thus eliminating field illumination
   %   inhomogeneity and non-cellular structures like blood vessels
   init_opt.secondary_metric    = R0;
   retuned_opt.secondary_metric = init_opt.secondary_metric;
   init_opt.second_alpha        = 0.5; % same as alpha

   %% Initialise by getting the candidate ROI seeds using the green channel
   % - this is activity-based initialisation
    phi_0    = initialiseROISegmentation_Manfredi(corr,red,R0, cellrad, alpha, init_opt,1,3);
    exp_ROIs = maxcells;

   %% Retune the alpha parameter
   retuned_alpha     = alpha;
   initial_masks_num = size(phi_0, 3);
   retuned_masks_num = initial_masks_num;
   prevstr = [];
   while (retuned_masks_num > exp_ROIs)
       retuned_alpha = retuned_alpha + 0.05;
       retuned_opt.second_alpha = retuned_alpha;
       str = sprintf( '\tAlpha Tuning: Current Alpha is %g', retuned_alpha );
       refreshdisp( str, prevstr );
       prevstr = str;
       % disp(['Alpha Tuning: Current Alpha is ', num2str(retuned_alpha)]);
       phi    = initialiseROISegmentation_Manfredi(corr,red,R0, cellrad, retuned_alpha, retuned_opt,1,3);
       retuned_masks_num = size(phi,3);
   end
   % disp(['Final alpha value is ', num2str(retuned_alpha)]);
   str = sprintf( '\tFinal alpha value is %g\n', retuned_alpha );
   refreshdisp( str, prevstr );

   if ~exist('phi','var')

       phi= phi_0;
   end

   %% Retune the lambda parameter
   % Initialise the level set algorithm parameters
   %seg_opt.lambda = 60; % this depends on the data
   seg_opt.mergeCorr = 0.95;  % correlation coefficient threshold above
                              % which two neigbouring ROIs will be merged

   seg_opt.mergeDuring = 1;
   retune_lambda = 0; %make this 1 if you want to retune lambda

   if retune_lambda
       % Tune lambda based on the data
       seg_opt.lambda = 60;
       [retuned_lambda, phi] = lambda_tune(phi_0, stack_g, cellrad, seg_opt);
       seg_opt.lambda = retuned_lambda;
       if seg_opt.lambda < 30
          seg_opt.lambda = 60;
       end

   else
       phi = phi_0;
       seg_opt.lambda = 60; %default value of lambda
   end


   % tic;

   %% Segmentation proper
   seg_opt.corrIm = init_g;
   seg_opt.meanIm = R0;
   %seg_opt = rmfield(seg_opt,'plot_progress');


  % [masks, tsG] = segment_cells(phi, stack_g(:,:,1:2:end), cellrad, seg_opt); 

   [masks, tsG] = segment_cells(phi, stack_g, cellrad, seg_opt);

%    mask_num = size(masks,3); % Number of detected ROIs

   % Plot masks on summary image and save plot
  % plotopts.plot_ids = 1; % set to 1 to view the ID number of the ROIs on the plot
   %fig = plotContoursOnSummaryImage(mean_imratio, masks, plotopts);
   %fname_fig = [filedir file '_2P_ROIs'];
   %savefig(fig,fname_fig,'compact');
%    saveas(fig,fname_fig,'pdf');
%  close(fig);
%
%    % Save masks in a mat file
%%fname_masks = [filedir file '_2P_segment_output.mat'];
%
% fname_masks = [file '_2P_segment_output.mat'];
% save(fname_masks,'tsG','masks','corr');

end
