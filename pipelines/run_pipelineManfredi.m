% Written by Ann Go
%
% This script runs the complete data processing pipeline for a single
% image file. Processing steps include:
% (1) motion correction (and dezippering)
% (2) roi segmentation
% (3) neuropil decontamination and timeseries extraction
% (4) spike extraction
% (5) tracking data extraction
% (6) place field mapping
%
% The sections labeled "USER:..." require user input
%
% Matlab version requirement: neuroSEE_neuropilDecon requires at least Matlab R2018

clear; % close all;
tic

%% USER: Set basic methods

test = false;                    % flag to use one of smaller files in test folder)
default = true;                 % flag to use default parameters
                                % flag to force
force = [false;...              % (1) motion correction even if motion corrected images exist
         true;...              % (2) roi segmentation
         false;...              % (3) force neuropil decontamination
         false;...              % (4) force spike extraction
         false;...              % (5) force tracking data extraction
         false];                % (6) force place field mapping

mcorr_method = 'normcorre';     % [normcorre,fftRigid] CaImAn NoRMCorre method, fft-rigid method (Katie's)
segment_method = 'ABLE';      % [ABLE,CaImAn]
dofissa = true;                 % flag to implement FISSA (when false, overrides force(3) setting)
manually_refine_spikes = false; % flag to manually refine spike estimates
manuall_refine_PFmap = false;   % flag to tweak PFmap including changing number of epochs or velocity threshold


%% Load module folders and define data directory

[data_locn,comp,err] = load_neuroSEEmodules(test);
if ~isempty(err)
%     beep
%     cprintf('Errors',err);
    %return mANFREDI
end
if ~strcmpi(comp,'mac')
    manually_refine_spikes = false;     % no gui when running in linuxbox or hpc
end
if manually_refine_spikes
    global spikes params
end

%% USER: Specify file

file = '20181016_09_44_06';


%% USER: Set parameters (if not using default)

if ~default
    params.fr = 30.9;                                       % imaging frame rate [default: 30.9]
    % motion correction
        % neurosee method
        if strcmpi(mcorr_method,'fftRigid')
            params.fftRigid.imscale = 1;             % image downsampling factor                                             [default: 1]
            params.fftRigid.Nimg_ave = 10;           % no. of images to be averaged for calculating pixel shift (zippering)  [default: 10]
            params.fftRigid.refChannel = 'green';    % channel to be used for calculating image shift (green,red)            [default: 'green']
            params.fftRigid.redoT = 300;             % no. of frames at start of file to redo motion correction for after 1st iteration [default: 300]
        end
        % NoRMCorre
        if strcmpi(mcorr_method,'normcorre')
            params.nonrigid = NoRMCorreSetParms(...
                'd1',512,...        % width of image [default: 512]  *Regardless of user-inputted value, neuroSEE_motioncorrect reads this
                'd2',512,...        % length of image [default: 512] *value from actual image
                'grid_size',[32,32],...     % default: [32,32]
                'overlap_pre',[32,32],...   % default: [32,32]
                'overlap_post',[32,32],...  % default: [32,32]
                'iter',1,...                % default: 1
                'use_parallel',false,...    % default: false
                'max_shift',50,...          % default: 50
                'mot_uf',4,...              % default: 4
                'bin_width',200,...         % default: 200
                'max_dev',3,...             % default: 3
                'us_fac',50,...             % default: 50
                'init_batch',200);          % default: 200
        end
    % ROI segmentation
        params.ROIsegment.cellrad = 10;            % expected radius of a cell (pixels)    [default: 10]
        params.ROIsegment.maxcells = 200;          % estimated number of cells in FOV      [default: 200]
        params.ROIsegment.df_prctile = 5;          % percentile to be used for estimating baseline [default: 5]
        params.ROIsegment.df_medfilt1 = 13;        % degree of smoothing for df_f          [default: 23]
    % neuropil correction
    if dofissa
        params.fissa.ddf_prctile = 5;         % percentile to be used for estimating baseline [default:5]
        params.fissa.ddf_medfilt1 = 17;       % degree of smoothing for ddf_f                 [default: 23]
    end
    % spike extraction
        params.spkExtract.bl_prctile = 85;         % percentile to be used for estimating baseline [default:85]
        params.spkExtract.spk_SNR = 1;             % spike SNR for min spike value            [default: 1]
        params.spkExtract.decay_time = 0.4;        % length of a typical transient in seconds [default: 0.4]
        params.spkExtract.lam_pr = 0.99;           % false positive probability for determing lambda penalty   [default: 0.99]
    % PF mapping
        params.PFmap.Nepochs = 1;             % number of epochs for each 4 min video [default: 1]
        params.PFmap.histsmoothFac = 10;      % Gaussian smoothing window for histogram extraction        [default: 10]
        params.PFmap.Vthr = 20;               % speed threshold (mm/s) Note: David Dupret uses 20 mm/s    [default: 20]
                                                  %                              Neurotar uses 8 mm/s
end

if test
    params.ROIsegment.maxcells = 60;
end


%% Default parameters

if default
    params = load( 'default_params.mat' );

    % Remove irrelevant parameters
    if strcmpi(mcorr_method,'normcorre')
        params = rmfield(params,'fftRigid');
%         fields = {'df_prctile','df_medfilt1'};
%         params.ROIsegment = rmfield(params.ROIsegment,fields);
    else
        params = rmfield(params,'nonrigid');
    end

    if ~dofissa
        params = rmfield(params,'fissa');
    end
end

% Fixed parameters
params.methods.mcorr_method = mcorr_method;
params.methods.segment_method = segment_method;
params.methods.dofissa = dofissa;
if str2double(file(1:4)) > 2018
    params.FOV = 490; % FOV area = FOV x FOV, FOV in um
else
    params.FOV = 330;
end


%% Check if file has been processed. If not, continue processing unless forced to overwrite
% existing processed data
% check(1:6) check for existing data in processing steps 1-6
% check(7) checks for existing mat file pooling all processed data for file


%data_locn = '/home/mc6017/Data/'; % data_locn on linux_box
check = checkforExistingProcData(data_locn, file, mcorr_method, segment_method, dofissa);
if check(7) % All processing steps have been done
    if ~any(force)
        beep
        str = sprintf( '%s: File has been processed with specified options. Skipping processing.\n', file );
        cprintf(str)
        %return
    end
end
if force(1)
    force([2:4,6]) = true; % because redoing motion correction step affects all succeeding steps except step 5
elseif force(2)
    force([3:4,6]) = true;
elseif force(3)
    force([4,6]) = true;
elseif force(4)>0
    force(6) = true;
end


%% Image files
% Load original image files if forced to do motion correction or if motion corrected files don't exist

if force(1) || ~check(1)
    [imG,imR] = load_imagefile( data_locn, file );
else
    imG = []; imR = [];
end


%% (1) Motion correction
% Saved in file folder:
%   motion corrected tif files
%   summary fig & pdf,
%   mat with fields
%       green.[ meanframe, meanregframe ]
%       red.[ meanframe, meanregframe ]
%       template
%       shifts
%       col_shift
%       params

if any([ any(force(1:2)), ~check(2), ~check(1) ])
    [imG, imR, ~, params] = neuroSEE_motionCorrect( imG, imR, data_locn, file, params, force(1) );
else
    fprintf('%s: Motion corrected files found. Skipping motion correction\n', file);
    imG = []; imR = [];
end


%% (2) ROI segmentation
% Saved in file folder:
%   correlation image with ROIs (fig, pdf)
%   mat with fields {tsG, df_f, masks, corr_image, params}


[tsG, df_f, masks, corr_image, params] = neuroSEE_segment( imG, mean(imR,3), data_locn, file, params, force(2) );


%% (3) Run FISSA to extract neuropil-corrected time-series
% Saved in file folder: mat file with fields {tsG, df_f, masks}

if dofissa
    [tsG, dtsG, ddf_f, params] = neuroSEE_neuropilDecon( masks, data_locn, file, params, force(3) );
else
    dtsG = [];
    ddf_f = [];
end


%% (4) Estimate spikes
% Saved in file folder: mat with fields {tsG, dtsG, df_f, ddf_f, spikes, params}

[spikes, params, fname_mat] = neuroSEE_extractSpikes( df_f, ddf_f, data_locn, file, params, force(4) );

% Save spike output if required
if any([ ~check(4), any(force([1,2,4])), and(force(3), dofissa) ])
    spike_output.spikes = spikes;
    spike_output.tsG = tsG;
    spike_output.df_f = df_f;
    if ~isempty(dtsG), spike_output.dtsG = dtsG; end
    if ~isempty(ddf_f), spike_output.ddf_f = ddf_f; end
    spike_output.params = params.spkExtract;
    save(fname_mat,'-struct','spike_output');
end

if manually_refine_spikes
    GUI_manually_refine_spikes( spikes, tsG, dtsG, df_f, ddf_f, data_locn, file, params, corr_image, masks );
    uiwait
end


% %% (5) Find tracking file then load it
% % Saved: fig & pdf of trajectory
% %        mat with fields {time, r, phi, x, y , speed, w, alpha, TTLout, filename}
%
% fname_track = findMatchingTrackingFile( data_locn, file, force(5) );
% trackData = load_trackfile( data_locn,file, fname_track, force(5) );
% if any(trackData.r < 100)
%     params.mode_dim = '2D'; % open field
%     params.PFmap.Nbins = [10, 10]; % number of location bins in [x y]
% else
%     params.mode_dim = '1D'; % circular linear track
%     params.PFmap.Nbins = 30;      % number of location bins
% end
%
%
% %% (6) Generate place field maps
% % Saved: fig & pdf of summary consisting of occMap, infoMap, spikeMap and placeMap
% %        mat file with fields {occMap, spikeMap, infoMap, placeMap, downData, activeData,...
% %                               placeMap_smooth, sorted_placeMap, sortIdx, params}
%
% if manually_refine_spikes
%     force(6) = true;
% end
%
% % if strcmpi(params.mode_dim,'2D')
% %     [occMap, spikeMap, infoMap, placeMap, placeMap_smooth, placeMap_asd, ...
% %      downData, activeData, params] ...
% %       = neuroSEE_mapPF( spikes, trackData, data_locn, file, params, force(6));
% % else
% %     [occMap, spikeMap, infoMap, placeMap, placeMap_smooth, placeMap_asd, ...
% %      downData, activeData, params, ...
% %      sorted_placeMap, sorted_placeMap_smooth, sorted_placeMap_asd, ...
% %      normsorted_placeMap, normsorted_placeMap_smooth, normsorted_placeMap_asd, ...
% %      placeMap_pertrial, normplaceMap_pertrial, PCidx, sortIdx]...
% %       = neuroSEE_mapPF( spikes, trackData, data_locn, file, params, force(6));
% % end
% %
%
% %% Save output. These are all the variables needed for viewing data with GUI
%
% % save(fname_allData,'file','corr_image','masks','tsG','dtsG','df_f','ddf_f','spikes',...
% %                     'fname_track','downData','activeData','occMap','spikeMap','infoMap',...
% %                     'placeMap','placeMap_smooth','params');


t = toc;
str = sprintf('%s: Processing done in %g hrs\n', file, round(t/3600,2));
cprintf(str)
