clear; close all;
tic

%% USER: Set basic settings
                    % Set to
test = 0;           % [0,1] 1: debug mode (this will use one of smaller files in test folder)
default = 1;        % [0,1] 1: default parameters
force = [0;...      % [0,1] 1: load raw images instead of tif, i.e. force raw to tif conversion of images
         0;...      % [0,1] 1: redo motion correction even if motion corrected images exist
         0];        % [0,1] 1: force pipeline to redo all steps after motion correction
                    %       0, processing step is skipped if output of said
                    %          step already exists in individual file directory.
mcorr_method = 1;   % [1,2] 1: CaImAn NoRMCorre method, 2: fft-rigid method (Katie's)
segment_method = 2; % [1,2] 1: ABLE, 2: CaImAn
fissa_yn = 1;       % [0,1] 1: implement FISSA, 0: skip step


%% USER: Specify file

file = '20181016_10_11_35'; 


%% Load module folders and define data directory

data_locn = load_neuroSEEmodules(test);


%% Check if file has been processed. If not, continue processing unless forced to overwrite 
% existing processed data

fname_allData = [data_locn 'Data/' file(1:8) '/Processed/' file '/' file '_allData.mat'];
if exist(fname_allData,'file')
    if ~any(force)
        str = sprintf( '%s: File has been processed. Skipping processing\n', file );
        cprintf(str)
        return
    else
        if force(2)>0
            force(3) = 1; % redoing motion correction step affects all succeeding steps
        end
    end
end


%% Load image files
% force = 1 forces raw images to be loaded
% force = 0 loads tif files if they exist, but reverts to raw images if tif
%           files don't exist

[imG,imR] = load_imagefile( data_locn, file, force(1) );


%% Set default parameters

if default
    params = load( 'default_params.mat' );
end

params.mcorr_method = mcorr_method;
if mcorr_method == 1
    fields = {'imscale','Nimg_ave','refChannel','redoT'};
    params = rmfield(params,fields);
else
    params = rmfield(params,'nonrigid');
end

params.segment_method = segment_method;
if segment_method == 1
    
else
    fields = {'cellrad','maxcells','satThresh','satTime','areaThresh','noiseThresh'};
    params = rmfield(params,fields);
end

params.fissa_yn = fissa_yn;


%% USER: Set parameters (if not using default)

if ~default
    % motion correction
        % neurosee method
        if mcorr_method == 2
            params.imscale = 1;         % image downsampling factor
            params.Nimg_ave = 14;       % no. of images to be averaged for calculating pixel shift (zippering)
            params.refChannel = 2;      % channel to be used for calculating image shift (1: red, 2: green)
            params.redoT = 300;         % no. of frames at start of file to redo motion correction for after 1st iteration
        end
        % NoRMCorre
        if mcorr_method == 1
            params.nonrigid = NoRMCorreSetParms(...
                        'd1',size(imG,1),...
                        'd2',size(imG,2),...
                        'grid_size',[32,32],...
                        'overlap_pre',[32,32],...
                        'overlap_post',[32,32],...
                        'iter',1,...
                        'use_parallel',false,...
                        'max_shift',30,...
                        'mot_uf',4,...
                        'bin_width',200,...
                        'max_dev',3,...
                        'us_fac',50,...
                        'init_batch',200);
        end
    % ROI segmentation
        % ABLE
        if segment_method == 1
            params.cellrad = 8;             % expected radius of a cell (pixels)
            params.maxcells = 200;          % estimated number of cells in FOV
            params.satThresh = 2900;        % ROI is accepted if its fluorescence is below satThresh
            params.satTime = 0.3;           %   satTime fraction of the time
            params.areaThresh = 70;         % min acceptable ROI area (pixels)
            params.noiseThresh = 490;       % max acceptable ROI noise level
        end
        % CaImAn
        if segment_method == 2
        end
    % spike extraction
        params.RsmoothFac = 23;         % smoothing window for R
        params.g = 0.91;                % fluorescence impulse factor (OASIS)
        params.lambda = 2.4;            % sparsity penalty (OASIS)
    % PF mapping
        params.mode_dim = 2;            % 1: 1D, 2: 2D 
        params.FOV = 490;               % FOV area = FOV x FOV, FOV in um
        params.mode_method = 2;         % 1: ASD, 2: histogram estimation
        params.imrate = 30.9;           % image scanning frame rate
        params.Nbins = 30;             % number of location bins
        params.Nepochs = 1;             % number of epochs for each 4 min video
        params.smoothFac = 5;           % Gaussian smoothing window for histogram estimation
        params.Vthr = 20;               %s speed threshold (mm/s) Note: David Dupret uses 20 mm/s
                                        %                               Neurotar uses 8 mm/s
end

if test
    params.maxcells = 60;  
end

%% Do motion correction
% Saved in file folder: 
% (1) motion corrected tif files 
% (2) summary fig & pdf, 
% (3) mat with fields 
%       green.[ meanframe, meanregframe ]
%       red.[ meanframe, meanregframe ]
%       template
%       shifts
%       params.nonrigid OR params.[ imscale, Nimg_ave, refChannel, redoT ]

[imG, imR, mcorr_output, params] = neuroSEE_motionCorrect( imG, imR, data_locn, file, mcorr_method, params, force(2) );


%% Use ABLE to extract ROIs and raw time series
% Saved: image with ROIs (fig, pdf), mat with fields {tsG, tsR, masks, mean_imratio, params}

% [tsG, tsR, masks, mean_imratio, params] = neuroSEE_segment( imG, imR, mcorr_output.red.meanregframe, ...
%                                                             data_locn, file, params, force );

%% Run FISSA to extract neuropil-corrected time-series


%% Calculate ratiometric Ca time series (R) and extract spikes
% Saved: mat with fields {R, spikes, params}

% [R, spikes, params] = neuroSEE_extractSpikes( tsG, tsR, data_locn, file, params, force );

%% Find tracking file then load it
% Saved: fig & pdf of trajectory
%        mat with fields {time, r, phi, x, y , speed, w, alpha, TTLout, filename}

% fname_track = findMatchingTrackingFile(data_locn, file, force);
% trackdata = load_trackfile(data_locn,file, fname_track, force);

%% Generate place field maps
% Saved: fig & pdf of summary consisting of occMap, infoMap, spikeMap and placeMap
%        mat file with fields {occMap, spikeMap, infoMap, placeMap, downData, activeData,...
%                               placeMap_smooth, sorted_placeMap, sortIdx, params}

% [occMap, spikeMap, infoMap, placeMap, downData, activeData, placeMap_smooth, sorted_placeMap, normsorted_placeMap, sortIdx,...
%     params] = neuroSEE_mapPF( spikes, trackdata, data_locn, file, params, force);

%% Save output. These are all the variables needed for viewing data with GUI

% save(fname_allData,'file','mcorr_output','tsG','tsR','masks','mean_imratio','R','spikes',...
%                     'fname_track','occMap','spikeMap','infoMap','placeMap','downData','activeData',...
%                     'placeMap_smooth','sorted_placeMap','sortIdx','params');

t = toc;
str = sprintf('%s: Processing done in %g hrs', file, round(t/3600,2));
cprintf(str)