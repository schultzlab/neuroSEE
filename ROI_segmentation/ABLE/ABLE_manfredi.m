%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   NeuroSEE: An automated Neuronal Source Extraction
%             and Exploration toolbox
%   
%   Author: Manfredi Castelli 
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


%function [tsG, tsR, masks, mean_imratio] = ABLE( stack_g, mean_r, cellrad, maxcells )

function [tsG, masks, overl_corr] = ABLE( stack_g, mean_r, file, cellrad, maxcells ) % by Ann
   %tic;
   str = sprintf('%s: Extracting ROIs with ABLE\n', file);
   cprintf(str);

%% Initialise the summary images
   szsmooth = [3 3];
   % Compute the summary images, i.e. mean and correlation images
   mean_g = mean(stack_g,3); 
   %tic;
   corr_g = crossCorr(stack_g(:,:,:)); %1:2:end));
   %elapsed = toc; % slack notification for manfredi's account
   %str = sprintf('MAC -> Cross Corr  for %d*%d image,finished processing in %2.2f minutes. No errors.',dim(1),dim(2),elapsed/60);
   %SendSlackNotification('https://hooks.slack.com/services/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',str, '#general','@manfredi.castelli17', [], [], []);
   
   mean_imratio = mean_g ./ mean_r;  
   % Set the initial images 
   init_g = medfilt2( corr_g, szsmooth );
   
   red = mat2gray(mean_r);
   red = medfilt2(red,[3 3]);
   red = imadjust(red);
   
   corr = mat2gray(corr_g);
   corr = medfilt2(corr,szsmooth);
   corr = imadjust(corr);
   
   overl_corr = imfuse(red,corr,'falsecolor','Scaling','joint','ColorChannels',[1 2 0]);
    
   clear corr; clear red;
  
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

   %% Initialise by getting the candidate ROI seeds using the green channel and red channel
   % - this is activity-based and morphological initialisation
   phi_0    = initialiseROISegmentation_Manfredi(init_g,mean_r, cellrad, alpha, init_opt,1,3);
   exp_ROIs = maxcells; 

   %% Retune the alpha parameter
   
   retuned_alpha     = alpha;
   initial_masks_num = size(phi_0, 3);
   retuned_masks_num = initial_masks_num;
   prevstr = [];
   if size(phi_0,3) < 20
   
   
       while (retuned_masks_num > exp_ROIs)
           retuned_alpha = retuned_alpha + 0.05;
           retuned_opt.second_alpha = retuned_alpha;
           str = sprintf( '\tAlpha Tuning: Current Alpha is %g', retuned_alpha );
           refreshdisp( str, prevstr );
           prevstr = str; 
           % disp(['Alpha Tuning: Current Alpha is ', num2str(retuned_alpha)]);
           phi = initialiseROISegmentation(init_g, cellrad, retuned_alpha, retuned_opt);
           retuned_masks_num = size(phi,3);
       end
   end
   % disp(['Final alpha value is ', num2str(retuned_alpha)]);
   str = sprintf( '\tFinal alpha value is %g\n', retuned_alpha );
   refreshdisp( str, prevstr );

%    if exist('phi','var')
%        phi_0_g = phi;
%    else 
% %        phi_0_g = phi_0; manfredi
%    end

   %% Initialise by getting the candidate ROI seeds using the red channel
   % - this is morphology-based initialisation
%    I_norm_ori = R0;
%    min_sigma  = 4;
%    max_sigma  = 8;
%    Ntheta     = 8;
%    C  = 0;
%    tm = 0;
%    p2 = 2*max_sigma*3 + 1;
%    bw = adaptivethreshold( I_norm_ori, p2, C, tm );
%    
%    [ new_coord_x, new_coord_y, absigma_leave] = SeedDetection( I_norm_ori, min_sigma, max_sigma, bw, Ntheta );
%    
%    MinArea = 6;
% 
%    outparams = fastSegmentation( I_norm_ori, new_coord_x, new_coord_y, absigma_leave, MinArea, Ntheta, min_sigma );
% 
%    [szX, szY] = size( R0 );
%    phi_0_r    = ones( szX, szY, length(outparams.xypos) );
%    for k=1:length( outparams.Pixels )
%        for j=1:length( outparams.Pixels{k} )
%            phi_0_r( outparams.Pixels{k}(2,j), outparams.Pixels{k}(1,j), k ) = -1;
%        end 
%    end




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
   seg_opt.corrIm = corr_g;
   seg_opt.meanIm = R0;
   %seg_opt = rmfield(seg_opt,'plot_progress');

   % [masks, tsG, ~, tsR, ~] ...
   %                     = segment_cells(phi, stack_g, mean_r, cellrad, seg_opt);
   [masks, tsG] = segment_cells(phi, stack_g, cellrad, seg_opt); % edit by Ann

%    runtime  = toc;
%    mask_num = size(masks,3); % Number of detected ROIs

   % Plot masks on summary image and save plot
%    plotopts.plot_ids = 1; % set to 1 to view the ID number of the ROIs on the plot
%    fig = plotContoursOnSummaryImage(mean_imratio, masks, plotopts);
%    fname_fig = [filedir file '_2P_ROIs'];
%    savefig(fig,fname_fig,'compact');
%    saveas(fig,fname_fig,'pdf');
%    close(fig);
%    
%    % Save masks in a mat file
%    fname_masks = [filedir file '_2P_segment_output.mat'];
%    save(fname_masks,'cell_tsG','cell_tsR','masks','mean_imratio');
    %elapsed = toc;
    %str = sprintf('MAC -> ABLE finished processing in %2.2f minutes. No errors.',elapsed/60);
    %SendSlackNotification('https://hooks.slack.com/services/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',str, '#general','@manfredi.castelli17', [], [], []);
end 
