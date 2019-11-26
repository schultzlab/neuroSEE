function frun_mcorrTOtracking_batch( array_id, list )

% Written by Ann Go
% 
% This script runs the following processing steps for each file in the
% specified list:
% (1) motion correction (and dezippering)
% (2) roi segmentation 
%
% The following steps are done if the Matlab version is R2018 or higher
% (3) neuropil decontamination and timeseries extraction
% (4) spike extraction
% (5) tracking data extraction
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The section labeled "USER-DEFINED SETTINGS" requires user input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Matlab version requirement: neuroSEE_neuropilDecon requires at least Matlab R2018


%% SETTINGS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USER-DEFINED SETTINGS                         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Basic settings
default = true;             % flag to use default parameters
                            % flag to force
force = [false;...              % (1) motion correction even if motion corrected images exist
         false;...              % (2) roi segmentation
         false;...              % (3) neuropil decontamination
         false;...              % (4) spike extraction
         false];                % (5) tracking data extraction

mcorr_method = 'normcorre';    % values: [normcorre, normcorre-r, normcorre-nr, fftRigid] 
                                    % CaImAn NoRMCorre method: 
                                    %   normcorre (rigid + nonrigid) 
                                    %   normcorre-r (rigid),
                                    %   normcorre-nr (nonrigid), 
                                    % fft-rigid method (Katie's)
segment_method = 'CaImAn';      % values: [ABLE,CaImAn]    
dofissa = true;                 % flag to implement FISSA (when false, overrides force(3) setting)
manually_refine_spikes = false; % flag to manually refine spike estimates
slacknotify = true;            % flag to send Ann slack notification re start and end of processing

            % Not user-defined
            % Load module folders and define data directory
            [data_locn,comp,err] = load_neuroSEEmodules;
            if ~isempty(err)
                beep
                cprintf('Errors',err);    
                return
            end
            % Some security measures
            if ~strcmpi(comp,'mac')
                manually_refine_spikes = false;     % Do not allow gui when running in linuxbox or hpc
            end
            if strcmpi(comp,'hpc')
                maxNumCompThreads(32);      % max # of computational threads, must be the same as # of ncpus specified in jobscript (.pbs file)
            end
            if manually_refine_spikes
                global spikes params
            end

% Processing parameters (if not using default)
if ~default
    params.fr = 30.9;                                % imaging frame rate [default: 30.9]
    % motion correction
        % Katie's method
        if strcmpi(mcorr_method,'fftRigid')
            params.fftRigid.imscale = 1;             % image downsampling factor                                             [default: 1]
            params.fftRigid.Nimg_ave = 10;           % no. of images to be averaged for calculating pixel shift (zippering)  [default: 10]
            params.fftRigid.refChannel = 'green';    % channel to be used for calculating image shift (green,red)            [default: 'green']
            params.fftRigid.redoT = 300;             % no. of frames at start of file to redo motion correction for after 1st iteration [default: 300]
        end
        % NoRMCorre-rigid
        if or(strcmpi(mcorr_method,'normcorre'),strcmpi(mcorr_method,'normcorre-r'))
            params.rigid = NoRMCorreSetParms(...
                'd1',512,...        % width of image [default: 512]  *Regardless of user-inputted value, neuroSEE_motioncorrect reads this 
                'd2',512,...        % length of image [default: 512] *value from actual image    
                'max_shift',20,...          % default: 50
                'bin_width',200,...         % default: 200
                'us_fac',50,...             % default: 50
                'init_batch',200);          % default: 200
        end
        % NoRMCorre-nonrigid
        if or(strcmpi(mcorr_method,'normcorre'),strcmpi(mcorr_method,'normcorre-nr') )    
            params.nonrigid = NoRMCorreSetParms(...
                'd1',512,...        % width of image [default: 512]  *Regardless of user-inputted value, neuroSEE_motioncorrect reads this 
                'd2',512,...        % length of image [default: 512] *value from actual image    
                'grid_size',[32,32],...     % default: [32,32]
                'overlap_pre',[32,32],...   % default: [32,32]
                'overlap_post',[32,32],...  % default: [32,32]
                'iter',1,...                % default: 1
                'use_parallel',false,...    % default: false
                'max_shift',15,...          % default: 50
                'mot_uf',4,...              % default: 4
                'bin_width',200,...         % default: 200
                'max_dev',3,...             % default: 3
                'us_fac',50,...             % default: 50
                'init_batch',200);          % default: 200
        end
    % ROI segmentation 
        params.ROIsegment.df_prctile = 5;     % percentile to be used for estimating baseline   [default: 5]
        params.ROIsegment.df_medfilt1 = 13;   % degree of smoothing for df_f                    [default: 23]
    % neuropil correction
    if dofissa
        params.fissa.ddf_prctile = 5;         % percentile to be used for estimating baseline   [default:5]
        params.fissa.ddf_medfilt1 = 17;       % degree of smoothing for ddf_f                   [default: 23]
    end
    % spike extraction
        params.spkExtract.bl_prctile = 85;    % percentile to be used for estimating baseline   [default:85]
        params.spkExtract.spk_SNR = 1;        % spike SNR for min spike value                   [default: 1]
        params.spkExtract.decay_time = 0.4;   % length of a typical transient in seconds        [default: 0.4]
        params.spkExtract.lam_pr = 0.99;      % false positive probability for determing lambda penalty   [default: 0.99]
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

params.methods.mcorr_method = mcorr_method;
params.methods.segment_method = segment_method;
params.methods.dofissa = dofissa;

% Default parameters
if default
    params = load_defaultparams(params);
end


%% Find details on files to be processed
tic 

listfile = [data_locn 'Digital_Logbook/lists/' list];
files = extractFilenamesFromTxtfile(listfile);
Nfiles = size(files,1); 

if array_id > Nfiles
    beep
    cprintf('Errors','Array_id exceeds number of files\n');    
    return
end
file = files(array_id,:); 

% Some auto-defined parameters
if str2double(files(1,1:4)) > 2018
    params.ROIsegment.cellrad = 6;          % expected radius of a cell (pixels)    
    params.ROIsegment.maxcells = 300;       % estimated number of cells in FOV      
else
    params.ROIsegment.cellrad = 9;            
    params.ROIsegment.maxcells = 200;       
end


%% Processing
% Send Ann slack message if processing has started
if slacknotify
    slacktext = [list ': processing started'];
    neuroSEE_slackNotify( slacktext );
end

% Some security measures
force = logicalForce(force);    % Only allow combinations of force values that make sense

% Check for existing processed data for specific file
check_file = checkforExistingProcData(data_locn, file, params);
if check(1:5)
    fprintf('%s: Processed data found. Skipping processing\n', file);
    return
end


%% (1) Motion correction and dezippering
% Load original image files if forced to do motion correction or if motion corrected files don't exist 
if force(1) || ~check_file(1)
    [imG,imR] = load_imagefile( data_locn, file );
else
    imG = []; imR = [];
end

% motion correction
if any([ any(force(1:2)), ~check_file(2), ~check_file(1) ]) 
    [imG, imR, ~, params] = neuroSEE_motionCorrect( imG, imR, data_locn, file, params, force(1) );
else 
    fprintf('%s: Motion corrected files found. Skipping motion correction\n', file);
    imG = []; imR = [];
end


%% (2) ROI segmentation
if strcmpi(segment_method,'CaImAn')
    clear imR; imR = [];
end
[tsG, df_f, masks, corr_image, params] = neuroSEE_segment( imG, mean(imR,3), data_locn, file, params, force(2) );
clear imG imR


%% Continue with next steps if Matlab version is at least R2018
release = version('-release'); % Find out what Matlab release version is running
MatlabVer = str2double(release(1:4));
if MatlabVer < 2018
    beep
    err = sprintf('%s: Higher Matlab version required; skipping FISSA and the rest of processing steps.\n', file);
    cprintf('Errors',err);    
    return
end


%% (3) Neuropil decontamination and timeseries extraction
if dofissa
    [tsG, dtsG, ddf_f, params] = neuroSEE_neuropilDecon( masks, data_locn, file, params, force(3) );
else
    dtsG = [];
    ddf_f = [];
end


%% (4) Spike estimation
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


%% (5) Tracking data extraction
fname_track = findMatchingTrackingFile( data_locn, file, force(5) );
load_trackfile( data_locn,file, fname_track, force(5) );


t = toc;
str = sprintf('%s: Processing done in %g hrs\n', file, round(t/3600,2));
cprintf(str)

% Send Ann slack message if processing has finished
if slacknotify
    slacktext = [file ': FINISHED in' num2str(round(t/3600,2)) ' hrs. No errors!'];
    neuroSEE_slackNotify( slacktext );
end

end



