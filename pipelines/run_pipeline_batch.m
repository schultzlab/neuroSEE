clear; close all;
tic

%% Basic setup
test = 1;       % set to 1 if testing, this will use one of smaller files in ../test
default = 1;    % set to 1 to use default parameters
force = 0;      % set to 1 to overwrite saved processed files. This will 
                %    force pipeline to redo all steps incl. raw to tif
                %    conversion. If force = 0, processing step is skipped
                %    if output of said step already exists in individual
                %    file directory.

%% Load module folders and define data directory
data_locn = load_neuroSEEmodules(test);


%% Read files

% [files, filesDim, filesFOV] = extractFilenamesFromTxtfile('test.txt');
files = extractFilenamesFromTxtfile_default('list_pipeline.txt');

%% Start processing

parfor i = 1:size(files,1)
    file = files(i,:);
    params = load('default_params.mat');
    params.mode_dim = filesDim(i);
    params.FOV = filesFOV(i);
    
    % Check if file has been processed. If not, continue processing unless forced to overwrite 
    % existing processed data
    fname_allData = fullfile(data_locn,'Data/',file(1:8),'/Processed/',file,'/',file,'_allData.mat');
    if exist(fname_allData,'file')
        if ~force
            str = sprintf( '%s: File has been processed. Skipping processing\n', file );
            cprintf(str)
        end
    end

    if ~exist(fname_allData,'file') || force
        
        % Load raw and save tif files or load tif files if they exist
        [imG,imR] = load_imagefile( data_locn, file, force );

        % Dezipper and do motion correction
        % Saved: tif files, summary fig & pdf, 
        %        mat with fields 
        %           green.template, meanframe, meanregframe, shift
        %           red.template, meanframe, meanregframe, shift
        %           params.imscale, Nimg_ave
        [imG, imR, mcorr_output, params] = neuroSEE_motionCorrect(imG, imR, data_locn, file, params, force );

        % Use ABLE to extract ROIs and raw time series
        % Saved: image with ROIs (fig, pdf), mat with fields {tsG, tsR, masks, mean_imratio, params}

        if test
            params.maxcells = 60;  
        end
        [tsG, tsR, masks, mean_imratio, params] = neuroSEE_segment( imG, mcorr_output.red.meanregframe, ...
                                                                    data_locn, file, params, force );

        % Run FISSA to extract neuropil-corrected time-series


        % Calculate ratiometric Ca time series (R) and extract spikes
        % Saved: mat with fields {R, spikes, params}
        [R, spikes, params] = neuroSEE_extractSpikes( tsG, tsR, data_locn, file, params, force );

        % Find tracking file then load it
        % Saved: fig & pdf of trajectory
        %        mat with fields {time, r, phi, x, y , speed, w, alpha, TTLout, filename}
        fname_track = findMatchingTrackingFile(data_locn, file, force);
        trackdata = load_trackfile(data_locn,file,fname_track,force);

        % Generate place field maps
        % Saved: mat file with fields {occMap, spikeMap, infoMap, placeMap, downData, activeData,...
        %                               placeMap_smooth, sorted_placeMap, sortIdx, params}
        [occMap, spikeMap, infoMap, placeMap, downData, activeData, placeMap_smooth, sorted_placeMap, sortIdx,...
            params] = neuroSEE_mapPF( spikes, trackdata, data_locn, file, params, force);

        % Save all data for 
        parsave(fname_allData,'file','mcorr_output','tsG','tsR','masks','mean_imratio','R','spikes',...
                        'fname_track','occMap','spikeMap','infoMap','placeMap','downData','activeData',...
                        'placeMap_smooth','sorted_placeMap','sortIdx','params');

        t = toc;
        str = sprintf('%s: Processing done in %g hrs', round(t/3600,2));
        cprintf(str)
    end
end