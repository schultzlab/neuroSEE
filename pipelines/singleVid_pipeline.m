clear; close all;
tic

%% Basic setup
test = 1;       % set to 1 if testing, this will use one of smaller files in ../test
default = 0;    % set to 1 to use default parameters
display = 0;    % set to 1 to display results (all results are saved in 
                %    individual file directories with summary pdfs regardless)
force = 0;      % set to 1 to overwrite saved processed files. This will 
                %    force pipeline to redo all steps incl. raw to tif
                %    conversion. If force = 0, processing step is skipped
                %    if output of said step already exists in individual
                %    file directory.

% gcp;           % start parallel pool
addpath(genpath('../behaviour'));
addpath(genpath('../intervideo_processing'));
addpath(genpath('../motion_correction'));
addpath(genpath('../PF_mapping'));
addpath(genpath('../ROI_segmentation'));
addpath(genpath('../spike_extraction'));
addpath(genpath('../utilities'));
addpath(genpath('../pipelines'));

% Data location
data_locn = '/Volumes/thefarm2/live/CrazyEights/AD_2PCa/';
if ~exist(data_locn,'dir')
    data_locn = '/Volumes/RDS/project/thefarm2/live/CrazyEights/AD_2PCa/';
end
if ~exist(data_locn,'dir')
    data_locn = '/rds/general/user/mgo/projects/thefarm2/live/CrazyEights/AD_2PCa/';
end

if test
    data_locn = [data_locn 'Data_forTesting/'];
end

    
%% Read files

file = '20181016_10_11_35'; 
% file = '20181016_09_09_43';

% file = '20190406_20_27_07'; 

%% Set parameters

if default
    params = load('default_params.mat');
end

if ~default
    % Specify own parameters
    % motion correction
        params.imscale = 1;             % image downsampling factor
        params.Nimg_ave = 14;           % number of images to be averaged for calculating pixel shift (zippering)
        params.refChannel = 2;
        params.redoT = 300;
    % ROI segmentation
        params.cellrad = 8;            % expected radius of a cell (pixels)
        params.maxcells = 200;          % estimated number of cells in FOV
        params.satThresh = 2900;        % ROI is accepted if its fluorescence is below satThresh
        params.satTime = 0.3;           %   satTime fraction of the time
        params.areaThresh = 70;         % min acceptable ROI area (pixels)
        params.noiseThresh = 490;       % max acceptable ROI noise level
    % spike extraction
        params.RsmoothFac = 23;         % smoothing window for R
        params.g = 0.91;                % fluorescence impulse factor (OASIS)
        params.lambda = 2.4;            % sparsity penalty (OASIS)
    % PF mapping
        params.mode_dim = 2;            % 1: 1D, 2: 2D 
        params.FOV = 490;               % FOV area = FOV x FOV, FOV in um
        params.mode_method = 2;         % 1: ASD, 2: histogram estimation
        params.imrate = 30.9;           % image scanning frame rate
        params.Nbins = 180;             % number of location bins
        params.Nepochs = 1;             % number of epochs for each 4 min video
        params.smoothFac = 5;           % Gaussian smoothing window for histogram estimation
        params.Vthr = 10;               %s speed threshold (mm/s) Note: David Dupret uses 20 mm/s
                                        %                               Neurotar uses 8 mm/s
end

if test
    params.maxcells = 60;  
end

%% Check if file has been processed. If not, continue processing unless forced to overwrite 
% existing processed data

fname_allData = [data_locn 'Data/' file(1:8) '/Processed/' file '/' file '_allData.mat'];
if exist(fname_allData,'file')
    if ~force
        str = sprintf( '%s: File has been processed. Skipping processing\n', file );
        cprintf(str)
        return
    end
end

%% Load raw and save tif files or load tif files if they exist

[imG,imR] = load_imagefile( data_locn, file, force );

%% Dezipper and do motion correction
% Saved: tif files, summary fig & pdf, 
%        mat with fields 
%           green.template, meanframe, meanregframe, shift
%           red.template, meanframe, meanregframe, shift
%           params.imscale, Nimg_ave

% [imG, imR, mcorr_output, params] = neuroSEE_motionCorrect(imG, imR, data_locn, file, params, force );
% params.nonrigid = NoRMCorreSetParms('d1',size(imG,1),'d2',size(imG,2),'grid_size',[32,32],'mot_uf',4,'bin_width',200,'max_shift',15,'max_dev',3,'us_fac',50,'init_batch',200);
% [imG, imR, shifts, template, options, col_shift] = normcorre_2ch( imG, imR, options);

%% Use ABLE to extract ROIs and raw time series
% Saved: image with ROIs (fig, pdf), mat with fields {tsG, tsR, masks, mean_imratio, params}

[tsG, tsR, masks, mean_imratio, params] = neuroSEE_segment( imG, imR, mcorr_output.red.meanregframe, ...
                                                            data_locn, file, params, force );

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
   params] = neuroSEE_mapPF( spikes, trackdata, data_locn, file, params, force);

%% Save output. These are all the variables needed for viewing data with GUI

% save(fname_allData,'file','mcorr_output','tsG','tsR','masks','mean_imratio','R','spikes',...
%                     'fname_track','occMap','spikeMap','infoMap','placeMap','downData','activeData',...
%                     'placeMap_smooth','sorted_placeMap','sortIdx','params');

t = toc;
str = sprintf('%s: Processing done in %g hrs', file, round(t/3600,2));
cprintf(str)