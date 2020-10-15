% Written by Ann Go
% Motion correction (NoRmCorre) segment adapted from NoRMCorreSetParms.m
% ROI segmentation (CaImAn) segment adapted from CNMFSetParms.m
% Struct for setting the neuroSEE pipeline parameters. Any parameter that is
% not set gets a default value

function params = neuroSEE_setparams(varargin)
Names = [
    % methods
        'mcorr_method       ' % motion correction method (default: 'normcorre')
        'segment_method     ' % ROI segmentation method: ABLE or CaImAn (default: CaImAn)
        'runpatches         ' % flag to run patches for CaImAn processing (default: false)
        'dofissa            ' % flag to do fissa correction (default: true)
        'doasd              ' % flag to calculate place fields by asd estimation  (default: false)
        'groupreg_method    ' % method for concatenating file data (either imreg or roireg) (default: [])
    % dataset info
        'd1                 ' % number of rows
        'd2                 ' % number of cols
        'd3                 ' % number of planes (for 3d imaging, default: 1)
        'FOV                ' % size of field of view in um
        'fr                 ' % imaging frame rate in Hz (defaut: 30)
        'decay_time         ' % length of a typical transient in seconds (default: 1.0)
    % motion correction (general)
        'refChannel         ' % reference channel for motion correction (default: 'green')
    % motion correction: fftRigid    
        'imscale            ' % image downsampling factor (default: 1)
        'Nimg_ave           ' % no. of images to be averaged for calculating pixel shift (zippering) (default: 10)
        'redoT              ' % no. of frames at start of file to redo motion correction for after 1st iteration (default: 300)
    % motion correction: NoRmCorre
        % patches
        'grid_size_r        ' % size of non-overlapping regions - rigid correction (default: [d1,d2,d3])
        'grid_size_nr       ' % size of non-overlapping regions - nonrigid correction (default: [d1/4,d2/4,d3])
        'max_shift_r        ' % maximum rigid shift in each direction - rigid correction (default: [30,30,1])
        'max_shift_nr       ' % maximum rigid shift in each direction - nonrigid correction (default: [30,30,1])
        'overlap_pre        ' % size of overlapping region (default: [64,64,1])
        'min_patch_size     ' % minimum size of patch (default: [64,64,1])    
        'min_diff           ' % minimum difference between patches (default: [32,32,1])
        'us_fac             ' % upsampling factor for subpixel registration (default: 50)
        'mot_uf_r           ' % degree of patches upsampling - rigid correction (default: [1,1,1])
        'mot_uf_nr          ' % degree of patches upsampling - nonrigid correction (default: [4,4,1])
        
        'max_dev            ' % maximum deviation of patch shift from rigid shift (default: [3,3,1])
        'overlap_post       ' % size of overlapping region after upsampling (default: [64,64,1])
        'phase_flag         ' % flag for using phase correlation (default: false)
        'shifts_method      ' % method to apply shifts ('FFT','cubic','linear')
        % template updating
        'upd_template       ' % flag for online template updating (default: true)
        'init_batch         ' % length of initial batch (default: 100)
        'bin_width          ' % width of each bin (default: 10)
        'buffer_width       ' % number of local means to keep in memory (default: 50)
        'method_avetemplate ' % method for averaging the template (default: {'median';'mean})
        'iter               ' % number of data passes (default: 1)
        'boundary           ' % method of boundary treatment 'NaN','copy','zero','template' (default: 'copy')
        % misc
        'add_value          ' % add dc value to data (default: 0)
        'use_parallel       ' % for each frame, update patches in parallel (default: false)
        'memmap             ' % flag for saving memory mapped motion corrected file (default: false)
        'mem_filename       ' % name for memory mapped file (default: 'motion_corrected.mat')
        'mem_batch_size     ' % batch size during memory mapping for speed (default: 5000)
        % plotting
        'plot_flag          ' % flag for plotting results in real time (default: false)
        'make_avi_mcorr     ' % flag for making movie (default: false)
        'name_mcorravi      ' % name for movie (default: 'motion_corrected.avi')
        % output type
        'output_type        ' % 'mat' (load in memory), 'memmap', 'tiff', 'hdf5', 'bin' (default:mat)
        'h5_groupname       ' % name for hdf5 dataset (default: 'mov')
        'h5_filename        ' % name for hdf5 saved file (default: 'motion_corrected.h5')
        'tiff_filename      ' % name for saved tiff stack (default: 'motion_corrected.tif')
        'output_filename    ' % name for saved file will be used if `h5_,tiff_filename` are not specified
        % use windowing
        'use_windowing      ' % flag for windowing data before fft (default: false)
        'window_length      ' % length of window on each side of the signal as a fraction of signal length
                               %    total length = length(signal)(1 + 2*window_length). (default: 0.5)
        % bitsize for reading .raw files
        'bitsize            ' % (default: 2 (uint16). other choices 1 (uint8), 4 (single), 8 (double))
        % offset from bidirectional sampling
        'correct_bidir      ' % check for offset due to bidirectional scanning (default: true)
        'nFrames            ' % number of frames to average (default: 50)
        'bidir_us           ' % upsampling factor for bidirectional sampling (default: 10)
        'col_shift          ' % known bi-directional offset provided by the user (default: [])
        'print_msg          ' % flag to print messages (default: false)
    % ROI segmentation (general)
        'cellrad_FOV490     ' % expected radius of a cell in 490x490 um FOV (pixels) (default: 6)  
        'maxcells_FOV490    ' % estimated number of cells in 490x490 um FOV (default: 500) 
        'cellrad_FOV330     ' % expected radius of a cell in 330x330 um FOV (pixels) (default: 9)  
        'maxcells_FOV330    ' % estimated number of cells in 330x330 um FOV (default: 400)
        'df_prctile         ' % percentile to be defined as baseline (default: 20)
        'roiarea_thr        ' % area of roi to be considered a cell (default: 70)
    % ROI segmentation (ABLE)    
        'df_medfilt1        ' % degree of smoothing for df_f (default: 13)
    % ROI segmentation: CaImAn
        % patches (optional)
        'patch_size         ' % size of each patch along each dimension (optional, default: [128,128])
        'overlap            ' % amount of overlap in each dimension (optional, default: [16,16])
        % INITIALIZATION  (initialize_components.m)
        'ssub               ' % spatial downsampling factor (default: 1)
        'tsub               ' % temporal downsampling factor (default: r)
        'init_method        ' % initialization method ('greedy','greedy_corr','sparse_NMF','HALS') (default: 'greedy')
        'rem_prct           ' % percentile to be removed before initialization (default: 20)
        'noise_norm         ' % normalization by noise estimate prior to initialization (default: true)
        'noise_norm_prctile ' % minimum noise level (as percentile of P.sn) used in the normalization prior to initialization (default: 2)
        % greedy_corr parameters (greedyROI_corr.m)
        'min_corr           ' % minimum local correlation for initializing a neuron (default: 0.3)
        % greedyROI parameters (greedyROI.m)
        'gSig               ' % half size of neurons to be found (default: [9,9])
        'gSiz               ' % half size of bounding box for each neuron (default: 2*gSig+1)
        'nb                 ' % number of background components (default: 1)
        'nIter              ' % maximum number of rank-1 NMF iterations during refining
        'med_app            ' % number of timesteps to be interleaved for fast (approximate) median calculation (default: 1)
        'save_memory        ' % process data sequentially to save memory (default: 0)
        'chunkSiz           ' % filter this number of timesteps each time (default: 100)
        'windowSiz          ' % size of window over which is computed sequentially (default: 32 x 32)
        'rolling_sum        ' % flag for using rolling sum to detect new components (default: True)
        'rolling_length     ' % length of rolling window (default: 100)
        % sparse_NMF parameters (sparse_NMF_initialization.m)
        'snmf_max_iter      ' % max # of sparse NMF iterations
        'err_thr            ' % relative change threshold for stopping sparse_NMF
        'eta                ' % frobenious norm factor *max(Y(:))^2
        'beta               ' % sparsity factor
        % HALS initialization parameters (HALS_initialization.m)
        'max_iter_hals_in   ' % maximum number of HALS iterations
        % HALS parameters (HALS_2d.m)
        'bSiz               ' % expand kernel for HALS growing (default: 3)
        'maxIter            ' % maximum number of HALS iterations (default: 5)
        % Noise and AR coefficients calculation (preprocess_data.m)
        'noise_range        ' % frequency range over which to estimate the noise (default: [0.25,0.5])
        'noise_method       ' % method for which to estimate the noise level (default: 'logmexp')
        'max_timesteps      ' % maximum number of timesteps over which to estimate noise (default: 3000)
        'flag_g             ' % compute global AR coefficients (default: false)
        'lags               ' % number of extra lags when computing the AR coefficients (default: 5)
        'include_noise      ' % include early lags when computing AR coefs (default: 0)
        'pixels             ' % pixels to include when computing the AR coefs (default: 1:numel(Y)/size(Y,ndims(Y)))
        'split_data         ' % split data into patches for memory reasons (default: 0)
        'block_size         ' % block size for estimating noise std in patches (default: [64,64])
        'cluster_pixels     ' % cluster pixels to active/inactive based on the PSD density (default: false)
        'extract_max        ' % extract the maximum activity intervals for each pixel (default: false)
        'max_nlocs          ' % number of local maxima to be extracted (default: 10)
        'max_width          ' % length of each interval (default: 11)
        % UPDATING SPATIAL COMPONENTS (unpdate_spatial_components.m)
        'spatial_method     ' % method for updating spatial components 'constrained' or 'regularized' (default: 'regularized')
        'search_method      ' % method for determining footprint of spatial components 'ellipse' or 'dilate' (default: 'dilate')
        'spatial_parallel   ' % update pixels in parallel (default: 1 if present)
        % determine_search_location.m
        'min_size           ' % minimum size of ellipse axis (default: 3)
        'max_size           ' % maximum size of ellipse axis (default: 8)
        'dist               ' % expansion factor of ellipse (default: 3)
        'se                 ' % morphological element for dilation (default: strel('disk',1,0))
        % threshold_components.m
        'thr_method         ' % method to threshold ('max' or 'nrg', default 'max')
        'maxthr             ' % threshold of max value below which values are discarded (default: 0.25)
        'nrgthr             ' % energy threshold (default: 0.995)
        'clos_op            ' % morphological element for closing (default: strel('square',3))
        'medw               ' % size of median filter (default: [3,3])
        'conn_comp          ' % extract largest connected component (binary, default: true)
        % UPDATING TEMPORAL COMPONENTS (update_temporal_components.m)
        'p                  ' % order of AR model dynamics (default: 2)
        'deconv_method      ' % method for spike deconvolution (default: 'constrained_foopsi')
        'restimate_g        ' % flag for updating the time constants for each component (default: 1)
        'temporal_iter      ' % number of block-coordinate descent iterations (default: 2)
        'temporal_parallel  ' % flag for parallel updating of temporal components (default: true if present)
        'full_A             ' % if true turn A into full matrix. If false turn Y into double precision (default: false)
        % CONSTRAINED DECONVOLUTION (constrained_foopsi.m)
        'spkinfer_method    ' % methods for performing spike inference ('dual','cvx','spgl1','lars') (default:'cvx')
        'bas_nonneg         ' % flag for setting the baseline lower bound. if 1, then b >= 0 else b >= min(y) (default 1)
        'fudge_factor       ' % scaling constant to reduce bias in the time constant estimation (default 1 - no scaling)
        'resparse           ' % number of times that the solution is resparsened (default: 0)
        % MERGING (merge_ROIs.m)
        'merge_thr          ' % merging threshold (default: 0.8)
        'fast_merge         ' % flag for using fast merging (default 1)
        % DF/F (extract_DF_F.m)
        'df_window          ' % length of running window (default [], no window)
        % CONTOUR PLOTS (plot_contours.m)
        'plot_bck_image     ' % plot background image or overlay on existing one (deafult: true)
        'cont_threshold     ' 
        % VIDEO (make_patch_video.m)
        'ind                ' % indices of components to be shown (deafult: 1:4)
        'skip_frame         ' % skip frames when showing the video (default: 1 (no skipping))
        'sx                 ' % half size of representative patches (default: 16)
        'make_avi_rois      ' % flag for saving avi video (default: 0)
        'show_background    ' % flag for displaying the background in the denoised panel (default: 1)
        'show_contours      ' % flag for showing the contour plots of the patches in the FoV (default: 0)
        'cmap               ' % colormap for plotting (default: 'default')
        'name_roiavi        ' % name of saved video file (default: based on current date)
        % PLOT COMPONENTS (view_patches.m)
        'plot_df            ' % flag for displaying DF/F estimates (default: 1)
        'make_gif           ' % save animation (default: 0)
        'save_avi           ' % save video (default: 0)
        'pause_time         ' % time to pause between each component (default: Inf, user has to click)
        % CLASSIFY COMPONENTS PIXELS (classify_components_pixels.m)
        'cl_thr             ' % overlap threshold for energy for a component to be classified as true (default: 0.8)
        % CLASSIFY COMPONENTS with CORRELATION (classify_comp_corr.m)
        'space_thresh       ' % threshold for r-value in space (default: 0.5)
        'time_thresh        ' % threshold for r-value in time (default: 0.4)
        'A_thresh           ' % threshold for determining overlap (default: 0.1)
        'Npeaks             ' % # of peaks to be considered (default: 5)
        'peak_int           ' % interval around the peak (default: -2:5)
        'MinPeakDist        ' % minimum peak distance for finding points of high activity  (default: 10)
        % ORDER COMPONENTS (order_components.m)
        'nsd                ' % number of standard deviations (default: 3)
        'nfr                ' % number of consecutive frames (default: 3)
        % PATCHES          (run_CNMF_patches.m)
        'gnb                ' % number of global background components (default: 3)
        'create_memmap      ' % create a memory mapped file if it is not provided in the input (default: false)    
        'classify_comp      ' % classify components based on correlation values (default: true)
        'refine_flag        ' % refine components within patch processing after merging (default: true)    
        'patch_space_thresh ' % space correlation threshold within patch (default: 0.25)
        'patch_time_thresh  ' % time correlation threshold within patch (default: 0.25)
        'patch_min_SNR      ' % minimum SNR for accepting exceptional events within patch (default: 0.5)
        'patch_min_fitness  ' % maximum fitness threshold within patch (default: log(normcdf(-patch_min_SNR))*N_samples_exc)
        'patch_min_fit_delta' % maximum fitness_delta threshold within patch (default: -2)
        'patch_cnn_thr      ' % threshold for CNN classifier within a patch (default: 0.05)
        % parameters for max probability test (trace_fit_extreme.m)
        'max_pr_thr         ' % threshold for keeping components (default: 0.9)
        't_int              ' % length of each trial in sec (default: 0.25)
        'sn_fac             ' % multiplicative factor for estimated noise level (default: 1)
        % parameters for thresholding based on size (classify_components.m)
        'max_size_thr       ' % maximum size of each component in pixels (default: 300)
        'min_size_thr       ' % minimum size of each component in pixels (default: 9)
        'size_thr           ' % fraction of max value for thresholding each component before determining its size (default 0.2)
        % parameters for registering components across different sessions (register_ROIs.m)
        'dist_exp           ' % exponent for calculating the distance between different ROIs (default: 1)
        'dist_thr           ' % distance threshold above which dist = Inf (default: 0.5)
        'dist_maxthr        ' % max thresholding for components before turing into binary masks (default: 0.15)
        'dist_overlap_thr   ' % threshold for detecting if one ROI is a subset of another (deafult: 0.8)
        'plot_reg           ' % plot registered ROIs (default: true)
        % parameters for computing event exceptionality (compute_event_exceptionality.m)
        'min_SNR            ' % minimum SNR for accepting exceptional events (default: 3)
        'robust_std         ' % use robust std for computing noise in traces (false)
        'N_samples_exc      ' % number of samples over which to compute (default: ceil(decay_time*fr))
        'min_fitness        ' % threshold on time variability  (default: log(normcdf(-min_SNR))*N_samples_exc)    
        'min_fitness_delta  ' % threshold on the derivative of time variability
        % parameters for CNN classifier (cnn_classifier.m)
        'cnn_thr            ' % threshold for CNN classifier (default: 0.2)
    % neuropil correction
        'ddf_prctile        ' % percentile to be used for estimating baseline (default: 5)
        'ddf_medfilt1       ' % degree of smoothing for ddf_f (default: 17)
    % spike extraction
        'bl_prctile         ' % percentile to be used for estimating baseline (default: 85)
        'spk_SNR            ' % spike SNR for min spike value (default: 1)
        'lam_pr             ' % false positive probability for determing lambda penalty (default: 0.99)
    % PF mapping
        'Nepochs            ' % number of epochs into which data will be divided for analysis (default: 1)
        'histsmoothWin      ' % smoothing window for histogram method (default: 7)
        'gaussfiltSigma     ' % stdev for 2d gaussian filter (default: 1.5)    
        'Vthr               ' % speed threshold (mm/s) Note: David Dupret uses 20 mm/s (default: 20)
                              %                              Neurotar uses 8 mm/s
        'prctile_thr        ' % percentile threshold for filtering nonPCs (default: 99)
        'pfactivet_thr      ' % fraction of dwell time in place field cell is required to be active (default: 0.05)
        'activetrials_thr   ' % fraction of trials cell is required to be active (default: 0.5)
        'Nrand              ' % number of shuffles for bootstrap test (default: 1000)
        'Nbins_1D           ' % no. of position bins in 103-cm linear track (default: 50)
        'Nbins_2D           ' % position bins in 325-mm diameter open field arena (default: [16,16])
   ]; 
   
[m,~] = size(Names);
names = lower(Names);

% Combine all leading options structures o1, o2, ... in l1Set(o1,o2,...).
options = [];
for j = 1:m
    eval(['options.' Names(j,:) '= [];']);
end
i = 1;
while i <= nargin
    arg = varargin{i};
    if ischar(arg), break; end
    if ~isempty(arg)                      % [] is a valid options argument
        if ~isa(arg,'struct')
            error(sprintf(['Expected argument %d to be a string parameter name ' ...
                'or an options structure\ncreated with OPTIMSET.'], i));
        end
        for j = 1:m
            if any(strcmp(fieldnames(arg),deblank(Names(j,:))))
                eval(['val = arg.' Names(j,:) ';']);
            else
                val = [];
            end
            if ~isempty(val)
                eval(['options.' Names(j,:) '= val;']);
            end
        end
    end
    i = i + 1;
end

% A finite state machine to parse name-value pairs.
if rem(nargin-i+1,2) ~= 0
    error('Arguments must occur in name-value pairs.');
end
expectval = 0;                          % start expecting a name, not a value
while i <= nargin
    arg = varargin{i};
    
    if ~expectval
        if ~ischar(arg)
            error(sprintf('Expected argument %d to be a string parameter name.', i));
        end
        
        lowArg = lower(arg);
        j = strmatch(lowArg,names);
        if isempty(j)                       % if no matches
            error(sprintf('Unrecognized parameter name ''%s''.', arg));
        elseif length(j) > 1                % if more than one match
            % Check for any exact matches (in case any names are subsets of others)
            k = strmatch(lowArg,names,'exact');
            if length(k) == 1
                j = k;
            else
                msg = sprintf('Ambiguous parameter name ''%s'' ', arg);
                msg = [msg '(' deblank(Names(j(1),:))];
                for k = j(2:length(j))'
                    msg = [msg ', ' deblank(Names(k,:))];
                end
                msg = sprintf('%s).', msg);
                error(msg);
            end
        end
        expectval = 1;                      % we expect a value next
        
    else
        eval(['options.' Names(j,:) '= arg;']);
        expectval = 0;
        
    end
    i = i + 1;
end

if expectval
    error(sprintf('Expected value for parameter ''%s''.', arg));
end

Values = [
    % methods
        {'normcorre'}         % motion correction method (default: 'normcorre')
        {'CaImAn'}            % ROI segmentation method: ABLE or CaImAn (default: CaImAn)
        {false}               % flag to run patches for CaImAn processing (default: false)
        {true}                % flag to do fissa correction (default: true)
        {false}               % flag to calculate place fields by asd estimation  (default: false)
        {[]}                  % method for concatenating file data (either imreg or roireg) (default: 'imreg')
    % dataset info
        {512}                 % number of rows
        {512}                 % number of columns
        {1}
        {330}                 % size of field of view in um
        {30.9}                % imaging frame rate in Hz (defaut: 30)
        {1.0}                 % length of a typical transient in seconds (default: 1.0)
    % motion correction (general)
        {'green'}             % reference channel for motion correction (default: 'green')
    % motion correction: fftRigid    
        {1}                   % image downsampling factor (default: 1)
        {10}                  % no. of images to be averaged for calculating pixel shift (zippering) (default: 10)
        {300}                 % no. of frames at start of file to redo motion correction for after 1st iteration (default: 300)
    % NoRmCorre
        % patches
        {[512,512,1]}         % size of non-overlapping regions - rigid correction (default: [d1,d2,d3])
        {[64,64,1]}         % size of non-overlapping regions - nonrigid correction (default: [d1/4,d2/4,d3])
        {[30,30,1]}           % maximum rigid shift in each direction - rigid correction (default: [30,30,1])
        {[30,30,1]}           % maximum rigid shift in each direction - nonrigid correction (default: [20,20,1])
        {[16,16,1]}           % size of overlapping region (default: [32,32,1])
        {[16,16,1]}           % minimum size of patch (default: [32,32,1])    
        {[8,8,1]}           % minimum difference between patches (default: [16,16,1])
        {50}                  % upsampling factor for subpixel registration (default: 50)
        {[1,1,1]}             % degree of patches upsampling - rigid correction (default: [1,1,1])
        {[4,4,1]}             % degree of patches upsampling - nonrigid correction (default: [4,4,1])
        {[5,5,1]}             % maximum deviation of patch shift from rigid shift (default: [3,3,1])
        {[64,64,1]}           % size of overlapping region after upsampling (default: [32,32,1])
        {false}               % use phase correlation (good for high SNR)
        {'FFT'}               % method for applying shifts ('FFT', 'linear', 'cubic')
        % template updating
        {true}                % flag for online template updating (default: true)
        {100}                 % length of initial batch (default: 100)
        {50}                  % width of each bin (default: 10)
        {50}                  % number of local means to keep in memory (default: 50)
        {{'median';'mean'}}   % method for averaging the template (default: {'median';'mean'}
        {1}                   % number of data passes (default: 1)
        {'copy'}              % method of boundary treatment (default: 'copy')
        % misc
        {0}                   % add dc value to data (default: 0)
        {false}               % for each frame, update patches in parallel (default: false)
        {false}               % flag for saving memory mapped motion corrected file (default: false)
        {'motion_corrected.mat'} % name for memory mapped file (default: 'motion_corrected.mat')
        {1000}                % batch size used during memory mapping for faster mapping
        % plotting
        {false}               % flag for plotting results in real time (default: false)
        {false}               % flag for making movie (default: false)
        {'motion_corrected.avi'} % name for movie (default: 'motion_corrected.avi')
        % output_type    
        {'mat'}
        {'mov'}
        {'motion_corrected.h5'}
        {'motion_corrected.tif'}
        {''}
        % use_windowing
        {false}
        {0.5}
        % bitsize for reading .raw files
        {2}
        % offset from bidirectional sampling
        {true}
        {50}
        {10}
        {[]}
        {false}               % flag to print messages (default: false)
    % ROI segmentation
        {6}                   % expected radius of a cell in 490x490 um FOV (default: 6)
        {500}                 % estimated number of cells in 490x490 um FOV (default: 500)
        {9}                   % expected radius of a cell in 330x330 um FOV (pixels) (default: 9)  
        {400}                 % estimated number of cells in 330x330 um FOV (default: 400)
        {20}                  % percentile to be defined as baseline (default 20)
        {70}                  % area of roi to be considered a cell (default: 70)
    % ROI segmentation (ABLE)
        {13}                  % degree of smoothing for df_f (default: 13)
    % ROI segmentation: CaImAn
        % patches (optional)
        {[128,128]}           % size of each patch along each dimension (optional, default: [128,128])
        {[16,16]}             % amount of overlap in each dimension (optional, default: [16,16])
        % INITIALIZATION  (initialize_components.m)
        {1}                   % spatial downsampling factor (default: 1)
        {5}                   % temporal downsampling factor (default: 5)
        {'greedy'}            % initialization method ('greedy','greedy_corr','sparse_NMF','HALS') (default: 'greedy')
        {20}                  % percentile to be removed before initialization (default: 20)
        {true}                % normalization by noise estimate prior to initialization (default: true)
        {2}                   % minimum noise level (as percentile of P.sn) used in the normalization prior to initialization (default: 2)
        % greedy_corr parameters (greedyROI_corr.m)
        {0.3}                 % minimum local correlation for initializing a neuron (default: 0.3)
        % greedyROI parameters (greedyROI.m)
        {9}                   % half size of neurons to be found (default: [9,9])
        {[]}                  % half size of bounding box for each neuron (default: 2*gSig+1)
        {1}                   % number of background components (default: 1)
        {5}                   % maximum number of rank-1 NMF iterations during refining
        {1}                   % number of timesteps to be interleaved for fast (approximate) median calculation (default: 1)
        {0}                   % process data sequentially to save memory (default: 0)
        {100}                 % filter this number of timesteps each time (default: 100)
        {[32,32]}             % size of window over which is computed sequentially (default: 32 x 32)
        {true}                % flag for using rolling sum to detect new components (default: True)
        {100}                 % length of rolling window (default: 100)
        % sparse_NMF parameters (sparse_NMF_initialization.m)
        {100}                 % max # of sparse NMF iterations
        {1e-4}                % relative change threshold for stopping sparse_NMF
        {1}                   % frobenious norm factor *max(Y(:))^2
        {.5}                  % sparsity factor
        % Noise and AR coefficients calculation (preprocess_data.m)
        {5}                   % maximum number of HALS iterations
        % HALS parameters (HALS_2d.m)
        {3}                   % expand kernel for HALS growing (default: 3)
        {5}                   % maximum number of HALS iterations (default: 5)
        % Noise and AR coefficients calculation (preprocess_data.m)
        {[0.25,0.5]}          % frequency range over which to estimate the noise (default: [0.25,0.5])
        {'mean'}              % method for which to estimate the noise level (default: 'logmexp')
        {3000}                % maximum number of timesteps over which to estimate noise (default: 3000)
        {false}               % compute global AR coefficients (default: false)
        {5}                   % number of extra lags when computing the AR coefficients (default: 5)
        {false}               % include early lags when computing AR coefs (default: 0)
        {[]}                  % pixels to include when computing the AR coefs (default: 1:numel(Y)/size(Y,ndims(Y)))
        {false}               % split data into patches for memory reasons (default: 0)
        {[64,64]}             % block size for estimating noise std in patches (default: [64,64])
        {false}               % cluster pixels to active/inactive based on the PSD density (default: false)
        {false}               % extract the maximum activity intervals for each pixel (default: false)
        {30}                  % number of local maxima to be extracted (default: 10)
        {21}                  % length of each interval (default: 11)
        % UPDATING SPATIAL COMPONENTS (unpdate_spatial_components.m)
        {'regularized'}       % method for updating spatial components 'constrained' or 'regularized' (default: 'regularized')
        {'dilate'}            % method for determining footprint of spatial components 'ellipse' or 'dilate' (default: 'dilate')
        {~isempty(which('parpool'))}    % update pixels in parallel (default: 1 if present)
        % determine_search_location.m
        {3}                   % minimum size of ellipse axis (default: 3)
        {8}                   % maximum size of ellipse axis (default: 8)
        {3}                   % expansion factor of ellipse (default: 3)
        {strel('disk',1,0)}   % morphological element for dilation (default: strel('disk',1,0))
        % CONSTRAINED DECONVOLUTION (constrained_foopsi.m)
        {'max'}               % method to threshold ('max' or 'nrg', default 'max')
        {0.25}                % threshold of max value below which values are discarded (default: 0.25)
        {0.995}               % energy threshold (default: 0.995)
        {strel('square',3)}   % morphological element for closing (default: strel('square',3))
        {[3,3]}               % size of median filter (default: [3,3])
        {true}                % extract largest connected component (binary, default: true)
        % UPDATING TEMPORAL COMPONENTS (update_temporal_components.m)
        {2}                   % order of AR model dynamics (default: 1)
        {'constrained_foopsi'} % method for spike deconvolution (default: 'constrained_foopsi')
        {1}                   % flag for updating the time constants for each component (default: 1)
        {4}                   % number of block-coordinate descent iterations (default: 2)
        {~isempty(which('parpool'))} % flag for parallel updating of temporal components (default: true if present)
        {false}               % if true turn A into full matrix. If false turn Y into double precision (default: false)
        % CONSTRAINED DECONVOLUTION (constrained_foopsi.m)
        {'cvx'}               % methods for performing spike inference ('dual','cvx','spgl1','lars') (default:'cvx')
        {1}                   % flag for setting the baseline lower bound. if 1, then b >= 0 else b >= min(y) (default 1)
        {0.98}                % scaling constant to reduce bias in the time constant estimation (default 1 - no scaling)
        {0}                   % number of times that the solution is resparsened (default: 0)
        % MERGING (merge_ROIs.m)
        {0.8}                % merging threshold (default: 0.8)
        {1}                   % flag for using fast merging (default 1)
        % DF/F (extract_DF_F.m)
        {[]}                  % length of running window (default [], no window)
        % CONTOUR PLOTS (plot_contours.m)
        {true}                % plot background image or overlay on existing one (deafult: true)
        {0.9}    
        % VIDEO (make_patch_video.m)
        {1:4}                 % indeces of components to be shown (deafult: 1:4)
        {1}                   % skip frames when showing the video (default: 1 (no skipping))
        {16}                  % half size of representative patches (default: 16)
        {0}                   % flag for saving avi video (default: 0)
        {1}                   % flag for displaying the background in the denoised panel (default: 1)
        {1}                   % flag for showing the contour plots of the patches in the FoV (default: 0)
        {'default'}           % colormap for plotting (default: 'default')
        {['video_',datestr(now,30),'.avi']} % name of saved video file (default: based on current date)
        % PLOT COMPONENTS (view_patches.m)
        {1}                   % flag for displaying DF/F estimates (default: 1)
        {0}                   % save animation (default: 0)
        {0}                   % save video (default: 0)
        {Inf}                 % time to pause between each component (default: Inf, user has to click)
        % CLASSIFY COMPONENTS PIXELS (classify_components_pixels.m)
        {0.8}                 % overlap threshold for energy for a component to be classified as true (default: 0.8)
        % CLASSIFY COMPONENTS with CORRELATION (classify_comp_corr.m)
        {0.5}                 % threshold for r-value in space (default: 0.5)
        {0.4}                 % threshold for r-value in time (default: 0.4)
        {0.1}                 % threshold for determining overlap (default: 0.1)
        {5}                   % # of peaks to be considered (default: 5)
        {-2:6}                % interval around the peak (default: -2:5)
        {10}                  % minimum peak distance for finding points of high activity  (default: 10)
        % ORDER COMPONENTS (order_components.m)
        {3}                   % number of standard deviations (default: 3)
        {5}                   % number of consecutive frames (default: 3)
        % PATCHES          (run_CNMF_patches.m)
        {3}                   % number of global background components (default: 3)
        {false}               % create a memory mapped file if it is not provided in the input (default: false)    
        {true}                % classify components based on correlation values (default: true)
        {true}                % refine components within patch processing after merging (default: true)    
        {0.25}                % space correlation threshold within patch (default: 0.25)
        {0.25}                % time correlation threshold within patch (default: 0.25)
        {0.5}                 % minimum SNR for accepting exceptional events within patch (default: 0.5)
        {[]}                  % maximum fitness threshold within patch (default: log(normcdf(-patch_min_SNR))*N_samples_exc)
        {-2}                  % maximum fitness_delta threshold within patch (default: -2)
        {0.05}                % threshold for CNN classifier within a patch (default: 0.05)
        % parameters for max probability test (trace_fit_extreme.m)
        {0.9}                 % threshold for keeping components (default: 0.9)
        {0.25}                % length of each trial in sec (default: 0.25)
        {1}                   % multiplicative factor for estimated noise level (default: 1)
        % parameters for thresholding based on size (classify_components.m)
        {320}                 % maximum size of each component in pixels (default: 300)
        {9}                   % minimum size of each component in pixels (default: 9)
        {0.2}                 % fraction of max value for thresholding each component before determining its size (default 0.2)
        % parameters for registering components across different sessions (register_ROIs.m)
        {1}                   % exponent for calculating the distance between different ROIs (default: 1)
        {0.5}                 % distance threshold above which dist = Inf (default: 0.5)
        {0.15}                % max thresholding for components before turing into binary masks (default: 0.15)
        {0.8}                 % threshold for detecting if one ROI is a subset of another (deafult: 0.8)
        {true}                % plot registered ROIs (default: true)
        % parameters for computing event exceptionality (compute_event_exceptionality.m)
        {3}                   % minimum SNR for accepting exceptional events (default: 3)
        {false}               % use robust std for computing noise in traces (false)
        {[]}                  % number of samples over which to compute (default: ceil(decay_time*fr))
        {[]}                  % threshold on time variability  (default: log(normcdf(-min_SNR))*N_samples_exc)    
        {-5}                  % threshold on the derivative of time variability
        % parameters for CNN classifier (cnn_classifier.m)
        {0.2}                 % threshold for CNN classifier (default: 0.2)    
    % neuropil correction
        {5}                   % percentile to be used for estimating baseline (default:5)
        {17}                  % degree of smoothing for ddf_f (default: 17)
    % spike extraction
        {85}                  % percentile to be used for estimating baseline (default: 85)
        {1}                   % spike SNR for min spike value (default: 1)
        {0.99}                % false positive probability for determing lambda penalty (default: 0.99)
    % PF mapping
        {1}                   % number of epochs into which data will be divided for analysis (default: 1)
        {7}                   % smoothing window for histogram method (default: 7)
        {1.5}                 % stdev for 2d gaussian filter (default: 1.5)
        {20}                  % speed threshold (mm/s) Note: David Dupret uses 20 mm/s (default: 20)
                              %                              Neurotar uses 8 mm/s
        {99}                  % percentile threshold for filtering nonPCs (default: 99)
        {0.05}                % fraction of dwell time in place field cell is required to be active (default: 0.05)
        {0.5}                 % fraction of trials cell is required to be active (default: 0.5)
        {1000}                % number of shuffles for bootstrap test
        {50}                  % no. of position bins in 103-cm linear track (default: 50)
        {[16,16]}             % position bins in 325-mm diameter open field arena (default: [16,16])
    ];

for j = 1:m
    if eval(['isempty(options.' Names(j,:) ')'])
        eval(['options.' Names(j,:) '= Values{j};']);
    end
end

if ~isempty(options.output_filename)
    out_type = options.output_type;
    [filepath,name,~] = fileparts(options.output_filename);
    output_filename = fullfile(filepath,name);
    if strcmpi(options.h5_filename,'motion_corrected.h5') && (strcmpi(out_type,'h5') || strcmpi(out_type,'hdf5'))
        options.h5_filename = [output_filename,'.h5'];
    end
    if strcmpi(options.tiff_filename,'motion_corrected.tif') && (strcmpi(out_type,'tif') || strcmpi(out_type,'tiff'))
        options.tiff_filename = [output_filename,'.tif'];
    end    
end

% NoRmCorre quality check
if isempty(options.d1); options.d1 = input('What is the total number of rows? \n'); end
if isempty(options.d2); options.d2 = input('What is the total number of columns? \n'); end
%if options.d3 == 1; nd = 2; else nd = 3; end
if isempty(options.grid_size_r); options.grid_size_r = [options.d1,options.d2,options.d3]; end
if isempty(options.grid_size_nr); options.grid_size_nr = [options.d1/4,options.d2/4,options.d3]; end
if length(options.grid_size_r) == 1; options.grid_size_r = options.grid_size_r*ones(1,3); end
if length(options.grid_size_r) == 2; options.grid_size_r(3) = 1; end
if length(options.grid_size_nr) == 1; options.grid_size_nr = options.grid_size_nr*ones(1,3); end
if length(options.grid_size_nr) == 2; options.grid_size_nr(3) = 1; end
if length(options.min_patch_size) == 1; options.min_patch_size = options.min_patch_size*ones(1,3); end
if length(options.min_patch_size) == 2; options.min_patch_size(3) = 1; end
if length(options.min_diff) == 1; options.min_diff = options.min_diff*ones(1,3); end
if length(options.min_diff) == 2; options.min_diff(3) = 1; end
if length(options.overlap_pre) == 1; options.overlap_pre = options.overlap_pre*ones(1,3); end
if length(options.overlap_pre) == 2; options.overlap_pre(3) = 1; end
if length(options.overlap_post) == 1; options.overlap_post = options.overlap_post*ones(1,3); end
if length(options.overlap_post) == 2; options.overlap_post(3) = 1; end
if length(options.max_shift_r) == 1; options.max_shift_r = options.max_shift_r*ones(1,3); end
if length(options.max_shift_r) == 2; options.max_shift_r(3) = 1; end
if length(options.max_shift_nr) == 1; options.max_shift_nr = options.max_shift_nr*ones(1,3); end
if length(options.max_shift_nr) == 2; options.max_shift_nr(3) = 1; end
if length(options.max_dev) == 1; options.max_dev = options.max_dev*ones(1,3); end
if length(options.max_dev) == 2; options.max_dev(3) = 1; end
if length(options.mot_uf_r) == 1; options.mot_uf_r = options.mot_uf_r*ones(1,3); end
if length(options.mot_uf_r) == 2; options.mot_uf_r(3) = 1; end
if length(options.mot_uf_nr) == 1; options.mot_uf_nr = options.mot_uf_nr*ones(1,3); end
if length(options.mot_uf_nr) == 2; options.mot_uf_nr(3) = 1; end
options.mot_uf_r(options.grid_size_r >= [options.d1,options.d2,options.d3]) = 1;
options.mot_uf_nr(options.grid_size_nr >= [options.d1,options.d2,options.d3]) = 1;

% CNMF quality check
if isempty(options.N_samples_exc); options.N_samples_exc = ceil(options.fr*options.decay_time); end
if isempty(options.min_fitness); options.min_fitness = log(normcdf(-options.min_SNR))*options.N_samples_exc; end
if isempty(options.patch_min_fitness); options.patch_min_fitness = log(normcdf(-options.patch_min_SNR))*options.N_samples_exc; end

% output
% methods
f = {'mcorr_method'; 'segment_method'; 'runpatches'; 'dofissa'; 'doasd'; ...
     'groupreg_method'};
fn = fieldnames(options);
for i = 1:length(fn)
    for j = 1:length(f)
        if strcmpi(fn{i},f{j})
            params.methods.(fn{i}) = options.(fn{i});
        end
    end
end
% motion correction
params.mcorr.refChannel = options.refChannel;
if strcmpi(options.mcorr_method,'fftRigid')
    f = {'imscale'; 'Nimg_ave'; 'redoT'};
    fn = fieldnames(options);
    for i = 1:length(fn)
        for j = 1:length(f)
            if strcmpi(fn{i},f{j})
                params.fftRigid.(fn{i}) = options.(fn{i});
            end
        end
    end
end
if strcmpi(options.mcorr_method,'normcorre') || strcmpi(options.mcorr_method,'normcorre-r')
    f = {'d1'; 'd2'; 'd3'; ...
        'overlap_pre'; 'min_patch_size'; 'min_diff'; 'us_fac'; 'max_dev';'overlap_post';...
        'phase_flag'; 'shifts_method'; 'upd_template'; 'init_batch'; 'bin_width'; 'buffer_width';...
        'iter'; 'boundary'; 'add_value'; 'use_parallel'; 'memmap';...
        'mem_filename'; 'mem_batch_size'; 'plot_flag'; ...
        'output_type'; 'h5_groupname'; 'h5_filename'; 'tiff_filename'; 'output_filename';...
        'use_windowing'; 'window_length'; 'bitsize'; 'correct_bidir'; 'nFrames';...
        'bidir_us'; 'col_shift'; 'print_msg'};
    fn = fieldnames(options);
    for i = 1:length(fn)
        for j = 1:length(f)
            if strcmpi(fn{i},f{j})
                params.mcorr.normcorre_r.(fn{i}) = options.(fn{i});
            end
        end
    end
    params.mcorr.normcorre_r.grid_size = options.grid_size_r;
    params.mcorr.normcorre_r.max_shift = options.max_shift_r;
    params.mcorr.normcorre_r.mot_uf = options.mot_uf_r;
    params.mcorr.normcorre_r.method = options.method_avetemplate;
    params.mcorr.normcorre_r.make_avi = options.make_avi_mcorr;
    params.mcorr.normcorre_r.name = options.name_mcorravi;
    params.mcorr.normcorre_r.fr = round(options.fr,-1);
end
if strcmpi(options.mcorr_method,'normcorre') || strcmpi(options.mcorr_method,'normcorre-nr')
    f = {'d1'; 'd2'; 'd3'; ...
        'overlap_pre'; 'min_patch_size'; 'min_diff'; 'us_fac'; 'max_dev';'overlap_post';...
        'phase_flag'; 'shifts_method'; 'upd_template'; 'init_batch'; 'bin_width'; 'buffer_width';...
        'method'; 'iter'; 'boundary'; 'add_value'; 'use_parallel'; 'memmap';...
        'mem_filename'; 'mem_batch_size'; 'plot_flag'; ...
        'output_type'; 'h5_groupname'; 'h5_filename'; 'tiff_filename'; 'output_filename';...
        'use_windowing'; 'window_length'; 'bitsize'; 'correct_bidir'; 'nFrames';...
        'bidir_us'; 'col_shift'; 'print_msg'};
    fn = fieldnames(options);
    for i = 1:length(fn)
        for j = 1:length(f)
            if strcmpi(fn{i},f{j})
                params.mcorr.normcorre_nr.(fn{i}) = options.(fn{i});
            end
        end
    end
    params.mcorr.normcorre_nr.grid_size = options.grid_size_nr;
    params.mcorr.normcorre_nr.max_shift = options.max_shift_nr;
    params.mcorr.normcorre_nr.mot_uf = options.mot_uf_nr;
    params.mcorr.normcorre_nr.method = options.method_avetemplate;
    params.mcorr.normcorre_nr.make_avi = options.make_avi_mcorr;
    params.mcorr.normcorre_nr.name = options.name_mcorravi;
    params.mcorr.normcorre_nr.fr = round(options.fr,-1);

end
% ROI segmentation
params.ROIsegment.roiarea_thr = options.roiarea_thr;
params.ROIsegment.FOV = options.FOV;
if options.FOV == 330
    maxcells = options.maxcells_FOV330;
    cellrad = options.cellrad_FOV330;
else
    maxcells = options.maxcells_FOV490;
    cellrad = options.cellrad_FOV490;
end
if strcmpi(options.segment_method,'ABLE')
    f = {'df_prctile'; 'df_medfilt1'; 'fr'};
    fn = fieldnames(options);
    for i = 1:length(fn)
        for j = 1:length(f)
            if strcmpi(fn{i},f{j})
                params.ROIsegment.ABLE.(fn{i}) = options.(fn{i});
            end
        end
    end
    params.ROIsegment.ABLE.maxcells = maxcells;
    params.ROIsegment.ABLE.cellrad = cellrad;
end
if strcmpi(options.segment_method,'CaImAn')
    f = {'d1'; 'd2'; 'd3'; 'fr'; 'decay_time'; 'df_prctile';...
        'patch_size'; 'overlap'; 'ssub'; 'tsub'; 'init_method'; 'rem_prct'; 'noise_norm'; ...
        'noise_norm_prctile'; 'min_corr'; 'gSig'; 'gSiz'; 'nb'; 'nIter'; 'med_app'; 'save_memory';...
        'chunkSiz'; 'windowSiz'; 'rolling_sum'; 'rolling_length'; 'snmf_max_iter';...
        'err_thr'; 'eta'; 'beta'; 'max_iter_hals_in'; 'bSiz'; 'maxIter';...
        'noise_range'; 'noise_method'; 'max_timesteps'; 'flag_g'; 'lags'; 'include_noise'; 'pixels';...
        'split_data'; 'block_size'; 'cluster_pixels'; 'extract_max'; 'max_nlocs'; 'max_width';...
        'spatial_method'; 'search_method'; 'spatial_parallel'; 'min_size'; 'max_size';...
        'dist'; 'se'; 'thr_method'; 'maxthr'; 'nrgthr'; 'clos_op'; 'medw'; 'conn_comp';...
        'p'; 'deconv_method'; 'restimate_g'; 'temporal_iter'; 'temporal_parallel'; 'full_A';...
        'bas_nonneg'; 'fudge_factor'; 'resparse'; 'merge_thr'; 'fast_merge'; 'df_window';...
        'plot_bck_image'; 'cont_threshold'; 'ind'; 'skip_frame'; 'sx'; 'show_background';...
        'show_contours'; 'cmap'; 'plot_df'; 'make_gif'; 'save_avi'; 'pause_time'; 'cl_thr';...
        'space_thresh'; 'time_thresh'; 'A_thresh'; 'Npeaks'; 'peak_int'; 'MinPeakDist'; 'nsd'; 'nfr';...
        'gnb'; 'create_memmap'; 'classify_comp'; 'refine_flag'; 'patch_space_thresh';...
        'patch_time_thresh'; 'patch_min_SNR'; 'patch_min_fitness'; 'patch_min_fit_delta';...
        'patch_cnn_thr'; 'max_pr_thr'; 't_int'; 'sn_fac'; 'max_size_thr'; 'min_size_thr';...
        'size_thr'; 'dist_exp'; 'dist_thr'; 'dist_maxthr'; 'dist_overlap_thr'; 'plot_reg';...
        'min_SNR'; 'robust_std'; 'N_samples_exc'; 'min_fitness'; 'min_fitness_delta';...
        'cnn_thr'; 'spk_SNR'; 'lam_pr'};
    fn = fieldnames(options);
    for i = 1:length(fn)
        for j = 1:length(f)
            if strcmpi(fn{i},f{j})
                params.ROIsegment.CaImAn.(fn{i}) = options.(fn{i});
            end
        end
    end
    params.ROIsegment.CaImAn.method = options.spkinfer_method;
    params.ROIsegment.CaImAn.make_avi = options.make_avi_rois;      
    params.ROIsegment.CaImAn.name = options.name_roiavi;   
    params.ROIsegment.CaImAn.maxcells = maxcells;
    params.ROIsegment.CaImAn.gSig = cellrad;
end
% neuropil correction
if options.dofissa
    f = {'ddf_prctile'; 'ddf_medfilt1'; 'fr'};
    fn = fieldnames(options);
    for i = 1:length(fn)
        for j = 1:length(f)
            if strcmpi(fn{i},f{j})
                params.fissa.(fn{i}) = options.(fn{i});
            end
        end
    end
end
% spike extraction
f = {'bl_prctile'; 'spk_SNR'; 'decay_time'; 'lam_pr'; 'fr'};
fn = fieldnames(options);
for i = 1:length(fn)
    for j = 1:length(f)
        if strcmpi(fn{i},f{j})
            params.spkExtract.(fn{i}) = options.(fn{i});
        end
    end
end
% PF mapping
f = {'fr'; 'Nepochs'; 'histsmoothWin'; 'gaussfiltSigma'; 'Vthr'; 'prctile_thr'; 'pfactivet_thr'; ...
      'activetrials_thr'; 'Nrand'; 'Nlaps_thr'; 'Nbins_1D'; 'Nbins_2D'};
fn = fieldnames(options);
for i = 1:length(fn)
    for j = 1:length(f)
        if strcmpi(fn{i},f{j})
            params.PFmap.(fn{i}) = options.(fn{i});
        end
    end
end

