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
% INPUT
% file : filename in the format 'YYYYMMDD_hh_mm_ss'
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The section labeled "USER-DEFINED INPUT" requires user input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Matlab version requirement for running on hpc: 
%   normcorre and CaImAn only work with Matlab R2017a. 
%   FISSA requires at least Matlab R2018

function frun_pipeline( file )

tic

%% SETTINGS

% Auto-defined FOV
if str2double(file(1,1:4)) > 2018
    FOV = 490;                  % FOV area = FOV x FOV, FOV in um
else
    FOV = 330; 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USER-DEFINED INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Basic settings
test = false;                   % flag to use one of smaller files in test folder)
force = [false;...              % (1) motion correction even if motion corrected images exist
         false;...              % (2) roi segmentation
         false;...              % (3) neuropil decontamination
         false;...              % (4) spike extraction
         false;...              % (5) tracking data extraction
         false];                % (6) place field mapping

mcorr_method = 'normcorre';     % values: [normcorre, normcorre-r, normcorre-nr, fftRigid] 
                                    % CaImAn NoRMCorre method: 
                                    %   normcorre (rigid + nonrigid) 
                                    %   normcorre-r (rigid),
                                    %   normcorre-nr (nonrigid), 
                                    % fft-rigid method (Katie's)
segment_method = 'CaImAn';      % [ABLE,CaImAn]    
dofissa = true;                 % flag to implement FISSA (when false, overrides force(3) setting)
manually_refine_spikes = false; % flag to manually refine spike estimates
doasd = false;                  % flag to do asd pf calculation 
slacknotify = false;            % flag to send Ann slack notifications about processing

% Processing parameters (any parameter that is not set gets a default value)
params = neuroSEE_setparams(...
            'mcorr_method', mcorr_method, ...
            'dofissa', dofissa,...
            'FOV', FOV); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

params.methods.mcorr_method = mcorr_method;
params.methods.segment_method = segment_method;
params.methods.dofissa = dofissa;
params.methods.doasd = doasd;

% Load module folders and define data directory
[data_locn,comp,err] = load_neuroSEEmodules(test);
if ~isempty(err)
    beep
    cprintf('Errors',err);    
    return
end

% Some security measures
if strcmpi(comp,'hpc')
    maxNumCompThreads(32);        % max # of computational threads, must be the same as # of ncpus specified in jobscript (.pbs file)
end

% Matlab version
release = version('-release'); % Find out what Matlab release version is running
MatlabVer = str2double(release(1:4));


%% Check if file has been processed. If not, continue processing unless forced to overwrite 
% existing processed data
% check(1:6) check for existing data in processing steps 1-6
% check(7) checks for existing mat file pooling all processed data for file

check = checkforExistingProcData(data_locn, file, params.methods);

% Some security measures
force = logicalForce(force);    % Only allow combinations of force values that make sense

if ~any(force) && check(7)
    fprintf('%s: File already processed\n', file)
    return
end


%% Image files
% Load original image files if forced to do motion correction or if motion corrected files don't exist 

if force(1) || ~check(1)
    % Continue only if Matlab version is R2017
    if strcmpi(comp,'hpc') && MatlabVer > 2017
        beep
        err = sprintf('%s: Lower Matlab version required; skipping processing.\n', [mouseid '_' expname]);
        cprintf('Errors',err);
        return
    end
    
    % Send Ann slack message if processing has started
    if slacknotify
        slacktext = [file ': Processing started'];
        neuroSEE_slackNotify( slacktext );
    end

    if strcmpi(mcorr_method,'normcorre') 
        filedir = [ data_locn 'Data/' file(1:8) '/Processed/' file '/mcorr_normcorre-r/' ];
        fname_tif_gr_mcorr = [filedir file '_2P_XYT_green_mcorr.tif'];
        fname_tif_red_mcorr = [filedir file '_2P_XYT_red_mcorr.tif'];
        % If mcorr_method is 'normcorre' and the rigid correction has been done,
        % no need to load original image file
        if exist(fname_tif_gr_mcorr,'file') && exist(fname_tif_red_mcorr,'file')
            imG = []; imR = [];
        else
            [ imG, imR ] = load_imagefile( data_locn, file, force(1) );
        end
    else
        [ imG, imR ] = load_imagefile( data_locn, file, force(1) );
    end
else
    imG = []; imR = [];
end

            
%% (1) Motion correction
% Saved in file folder: motion corrected tif files
%                       summary fig & png, 
%                       mat with fields 
%                           green.[ meanframe, meanregframe ]
%                           red.[ meanframe, meanregframe ]
%                           template
%                           shifts
%                           col_shift
%                           params

if any([ force(2), ~check(2), ~check(1) ]) 
    if strcmpi(segment_method,'CaImAn') % CaImAn does not use imR
        [ imG, ~, params.mcorr ] = neuroSEE_motionCorrect( imG, imR, data_locn, file, mcorr_method, params.mcorr, [], force(1) );
        imR = [];
    else
        [ imG, ~, params.mcorr, imR ] = neuroSEE_motionCorrect( imG, imR, data_locn, file, mcorr_method, params.mcorr, [], force(1) );
    end
else 
    fprintf('%s: Motion corrected files found. Skipping motion correction\n', file);
    imG = []; imR = [];
end


%% (2) ROI segmentation
% Saved in file folder: correlation image with ROIs (fig, png) 
%                       summary plots of tsG, df_f (fig, png)
%                       mat with fields {tsG, df_f, masks, corr_image, params}

[tsG, df_f, masks, corr_image, params] = neuroSEE_segment( imG, data_locn, file, params, force(2), mean(imR,3) );
clear imG imR


%% (3) Run FISSA to extract neuropil-corrected time-series
% Saved in file folder: mat file with fields {dtsG, ddf_f, masks}
%                       summary plots of tsG, ddf_f (fig & png)

if dofissa
    release = version('-release'); % Find out what Matlab release version is running
    MatlabVer = str2double(release(1:4));
    if MatlabVer > 2017
        [tsG, dtsG, ddf_f, params] = neuroSEE_neuropilDecon( masks, data_locn, file, params, force(3) );
    else
        fprintf('%s: Higher Matlab version required, skipping FISSA.', file);
        dofissa = false;
        dtsG = [];
        ddf_f = [];
    end
else
    dtsG = [];
    ddf_f = [];
end


%% (4) Estimate spikes
% Saved in file folder: mat with fields {tsG, dtsG, df_f, ddf_f, spikes, params}
%                       summary plot of spikes (fig & png)

[spikes, params] = neuroSEE_extractSpikes( df_f, ddf_f, data_locn, file, params, force(4) );

if manually_refine_spikes
    GUI_manually_refine_spikes( spikes, tsG, dtsG, df_f, ddf_f, data_locn, file, params, corr_image, masks );
    uiwait 
end


%% (5) Find tracking file then load it
% Saved in file folder: trajectory plot (fig & png)
%                       mat with fields {time, r, phi, x, y , speed, w, alpha, TTLout, filename}

trackfile = findMatchingTrackingFile( data_locn, file, force(5) );
if force(5) || ~check(5)
    Trackdata = load_trackfile(data_locn, file, trackfile, force(5));
    downTrackdata = downsample_trackData( Trackdata, size(spikes,2), params.fr );
    save([data_locn 'Data/' file(1:8) '/Processed/' file '/behaviour/' file '_downTrackdata.mat'],'downTrackdata');
else
    M = load([data_locn 'Data/' file(1:8) '/Processed/' file '/behaviour/' file '_downTrackdata.mat']);
    downTrackdata = M.downTrackdata;
end


%% (6) Generate place field maps
% Saved: fig & png of several figures showing occMap, spkMap, pfMap, remapping and spkMap_pertrial
%        mat file with fields {occMap, hist, asd, downData, activeData}

if manually_refine_spikes
    force(6) = true;
end

if any(downTrackdata.r < 100)
    params.mode_dim = '2D';                     % open field
    params.PFmap.Nbins = params.PFmap.Nbins_2D; % number of location bins in [x y]               
else 
    params.mode_dim = '1D';                     % circular linear track
    params.PFmap.Nbins = params.PFmap.Nbins_1D; % number of location bins               
end
fields = {'Nbins_1D','Nbins_2D'};
params.PFmap = rmfield(params.PFmap,fields);

if strcmpi(params.mode_dim,'1D')
    [ hist, asd, pfData, activeData, params ] = neuroSEE_mapPF( spikes, downTrackdata, data_locn, file, params, force(6));
else
    [ occMap, hist, asd, downData, activeData, params, spkMap, spkIdx ] = neuroSEE_mapPF( spikes, trackData, data_locn, file, params, force(6));
end

%% Save output. These are all the variables needed for viewing data with GUI

if dofissa
    str_fissa = 'FISSA';
else
    str_fissa = 'noFISSA';
end
fname_allData = [ data_locn 'Data/' file(1:8) '/Processed/' file '/' file '_' mcorr_method '_' segment_method '_' str_fissa '_allData.mat'];

save(fname_allData,'file','corr_image','masks','tsG','df_f','spikes','trackfile',...
                    'downTrackdata','activeData','pfData','hist','asd','params');
if ~isempty(dtsG), save(fname_allData,'-append','dtsG'); end
if ~isempty(ddf_f), save(fname_allData,'-append','ddf_f'); end
 
t = toc;
str = sprintf('%s: Processing done in %g hrs\n', file, round(t/3600,2));
cprintf(str)

% Send Ann slack message if processing has finished
if slacknotify
    slacktext = [file ': FINISHED in' num2str(round(t/3600,2)) ' hrs. No errors!'];
    neuroSEE_slackNotify( slacktext );
end

end

