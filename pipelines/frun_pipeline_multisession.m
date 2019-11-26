function frun_pipeline_multisession( list )

% Written by Ann Go
% 
% This script runs the complete data processing pipeline for a list of
% image files typically corresponding to one experiment. Processing steps
% for each file include:
% (1) motion correction (and dezippering)
% (2) roi segmentation 
% (3) neuropil decontamination and timeseries extraction
% (4) spike extraction
% (5) tracking data extraction
% 
% Then across image files, the processing steps are 
% (6) roi registration 
% (7) spike/tracking data consolidation
% (8) place field mapping
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
         false;...              % (5) tracking data extraction
         false;...              % (6) roi registration across images
         false;...              % (7) spike/tracking data consolidation across images
         false];                % (8) place field mapping

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
                'max_shift',20,...          % default: 50
                'bin_width',200,...         % default: 200
                'us_fac',50,...             % default: 50
                'init_batch',200);          % default: 200
        end
        % NoRMCorre-nonrigid
        if or(strcmpi(mcorr_method,'normcorre'),strcmpi(mcorr_method,'normcorre-nr') )    
            params.nonrigid = NoRMCorreSetParms(...
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
    % roi registration
        params.maxthr = [];                     
        params.dist_maxthr = 0.1;             % threshold for turning spatial components into binary masks [default: 0.1]
        params.dist_exp = 0.5;                % power n for distance between masked components: dist = 1 - (and(m1,m2)/or(m1,m2))^n [default: 1]
        params.dist_thr = 0.7;                % threshold for setting a distance to infinity    [default: 0.5]
        params.dist_overlap_thr = 0.5;        % overlap threshold for detecting if one ROI is a subset of another [default: 0.8]
        params.plot_reg = true;
        params.print_msg = true;
        params.nonrigid.print_msg = false;
    % PF mapping
        params.PFmap.Nepochs = 1;             % number of epochs for each 4 min video           [default: 1]
        params.PFmap.histsmoothFac = 7;       % Gaussian smoothing window for histogram extraction        [default: 7]
        params.PFmap.Vthr = 20;               % speed threshold (mm/s) Note: David Dupret uses 20 mm/s    [default: 20]
                                              %                              Neurotar uses 8 mm/s
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

% MouseID and experiment name
[ mouseid, expname ] = find_mouseIDexpname(list);

% Some auto-defined parameters
if str2double(files(1,1:4)) > 2018
    params.FOV = 490;                       % FOV area = FOV x FOV, FOV in um
    params.ROIsegment.cellrad = 6;          % expected radius of a cell (pixels)    
    params.ROIsegment.maxcells = 300;       % estimated number of cells in FOV      
else
    params.FOV = 330; 
    params.ROIsegment.cellrad = 9;            
    params.ROIsegment.maxcells = 200;       
end


%% Processing
% Send Ann slack message if processing has started
if slacknotify
    slacktext = [list ': processing started'];
    neuroSEE_slackNotify( slacktext );
end

% Check for existing processed data
check_list = checkforExistingProcData(data_locn, list, params);
if check_list(4)
    fprintf('%s: Processed data found. Skipping processing\n', list);
    return
end

% Some security measures
force = logicalForce(force);    % Only allow combinations of force values that make sense

for i = 1:Nfiles
    file = files(i,:);
    
    %% (1) Motion correction and dezippering
    % Load original image files if forced to do motion correction or if motion corrected files don't exist 
    if force(1) || ~check_file(1)
        [imG,imR] = load_imagefile( data_locn, file );
    else
        imG = []; imR = [];
        t = load([data_locn 'Data/' file(1:8) '/Processed/' file '/' mcorr_method file '_mcorr_output.mat']);
        templates{i} = t.template;
    end
    
    % motion correction
    if any([ any(force(1:2)), ~check_file(2), ~check_file(1) ]) 
        if ~isempty(imG)
            if or(strcmpi(mcorr_method,'normcorre'),strcmpi(mcorr_method,'normcorre-r'))
                params.rigid.d1 = size(imG,1);
                params.rigid.d2 = size(imG,2);
                params.rigid.grid_size = [params.rigid.d1,params.rigid.d2];
            end
            if or(strcmpi(mcorr_method,'normcorre'),strcmpi(mcorr_method,'normcorre-nr'))
                params.nonrigid.d1 = size(imG,1);
                params.nonrigid.d2 = size(imG,2);
            end
        end
        [imG, imR, ~, params] = neuroSEE_motionCorrect( imG, imR, data_locn, file, params, force(1) );
    else 
        fprintf('%s: Motion corrected files found. Skipping motion correction\n', file);
        imG = []; imR = [];
    end
    
    %% (2) ROI segmentation
    if strcmpi(segment_method,'CaImAn')
        clear imR; imR = [];
    end
    [tsG, df_f, masks{i}, corr_image, params] = neuroSEE_segment( imG, mean(imR,3), data_locn, file, params, force(2) );
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
        [tsG, dtsG, ddf_f, params] = neuroSEE_neuropilDecon( masks{i}, data_locn, file, params, force(3) );
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
    trackData = load_trackfile( data_locn,file, fname_track, force(5) );
    if any(trackData.r < 100)
        params.mode_dim = '2D';         % open field
        params.PFmap.Nbins = [16, 16];  % number of location bins in [x y]               
    else 
        params.mode_dim = '1D';         % circular linear track
        params.PFmap.Nbins = 30;        % number of location bins               
    end
end


%% (6) Load ROIs and templates from each session 
fdir = [ data_locn 'Analysis/' mouseid '/summaries based on registered ROIs/' mouseid '_' expname '_ref_' files(1,:) '/'];
    if ~exist(fdir,'dir'), mkdir(fdir); end
fname = [fdir mouseid '_' expname '_ref_' files(1,:) '_registered_rois.mat' ];

if ~exist(fname,'file') || force
    fprintf('%s: loading ROIs and templates\n',[mouseid '_' expname]);
    for jj = 1:Nfiles
        file = files(jj,:); 
        
        A_temp = zeros(size(masks,1)*size(masks,2),size(masks,3));
        for ii = 1:size(masks,3)
            masks_temp = masks(:,:,ii);
            A_temp(:,ii) = masks_temp(:);
            outlines{ii,jj} = bwboundaries(masks(:,:,ii));    % boundary of each ROI
        end
        A{jj} = sparse(A_temp); % masks

        
    end

    %% Register ROIs across sessions
    params.d1 = size(masks,1);
    params.d2 = size(masks,2);
    params_mc = NoRMCorreSetParms(...
                'd1',params.d1,...        % width of image [default: 512]  *Regardless of user-inputted value, neuroSEE_motioncorrect reads this 
                'd2',params.d2,...        % length of image [default: 512] *value from actual image    
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
    params_mc.print_msg = false;

    fprintf('%s: registering ROIs\n',[mouseid '_' expname]);
    [A_union, assignments, matchings] = register_multisession(A, params, templates, params_mc);
    masks = reshape(full(A_union), params.d1, params.d2, size(A_union,2));
    Nrois = size(masks,3);

    % save masks
    registered_rois.masks = masks;
    registered_rois.outlines = outlines;
    registered_rois.assignments = assignments;
    registered_rois.matchings = matchings;
    registered_rois.params = params;
    registered_rois.params_mc = params_mc;

    % fprintf('%s: saving registered ROIs\n',[mouseid '_' expname]);
    % save(fname, '-struct', 'registered_rois')

    % superposition plots of matched ROIs
    [nRow, nCol] = getnRownCol(Nrois);
    nPlot = nRow*nCol;

    if fsave
        if ~exist([fdir 'registered_rois/'],'dir')
            mkdir([fdir 'registered_rois/']);
        end
    end
    for n = 0:Nrois/nPlot
        fh = figure;
        ha = tight_subplot(nRow,nCol,[.01 .01],[.01 .05],[.01 .01]);
        for ii = 0:nPlot-1
            if (n*nPlot+ii+1) <= Nrois
                axes(ha(ii+1));
                imshow(zeros(512,512)); hold on
                for jj = 1:Nfiles
                    k = assignments(n*nPlot+ii+1,jj);
                    if ~isnan(k)
                        plot(outlines{k,jj}{1}(:,2),outlines{k,jj}{1}(:,1),'Linewidth',1.5); hold on
                    end
                end
                axis off; title(['Cell ' num2str(n*nPlot+jj+1)],'fontsize',15);
            end
        end
%             if Nrois/nPlot <= 1
%                 fname_fig = [fdir 'registered_rois/' mouseid '_' expname '_ref_' files(1,:) '_registered_rois'];
%             else
%                 fname_fig = [fdir 'registered_rois/' mouseid '_' expname '_ref_' files(1,:) '_registered_rois_' num2str(n+1)];
%             end
%             savefig( fh, fname_fig );
%             saveas( fh, fname_fig, 'png' );
%             close( fh );
    end 

else
    fprintf('%s: loading registered ROIs\n',[mouseid '_' expname]);
    load(fname)
    Nrois = size(masks,3);
end

%% Collate the data for tsG, dtsG, df_f, ddf_f, spikes for each roi in A_union
fname = [fdir mouseid '_' expname '_ref_' files(1,:) '_spikes_tracking_data.mat' ];

if ~exist(fname,'file') || force
    fprintf('%s: concatinating timeseries and tracking data\n',[mouseid '_' expname]);
    downr_all = [];
    for jj = 1:Nfiles
        file = files(jj,:);
        ts{jj} = load([data_locn 'Data/' file(1:8) '/Processed/' file '/mcorr_normcorre/CaIman/FISSA/' file '_spikes_output.mat']);
        Nt(jj) = size(ts{jj}.spikes,2);

        trackfile = dir([data_locn 'Data/' file(1:8) '/Processed/' file '/behaviour/*.mat']);
        td = load([data_locn 'Data/' file(1:8) '/Processed/' file '/behaviour/' trackfile.name]);
        x = td.x;
        y = td.y;
        r = td.r;
        phi = td.phi;
        speed = td.speed;
        tracktime = td.time;

        % Pre-process tracking data
        t0 = tracktime(1);                  % initial time in tracking data

        % Convert -180:180 to 0:360
        if min(phi)<0
           phi(phi<0) = phi(phi<0)+360;
        end

        % generate imaging timestamps using known image frame rate
        fr = 30.9;
        dt = 1/fr;
        t = (t0:dt:Nt(jj)*dt)';
        if length(t) ~= Nt(jj)
            t = (t0:dt:(Nt(jj)+1)*dt)';
        end

        % Downsample tracking to Ca trace
        downData{jj}.phi = interp1(tracktime,phi,t,'linear');
        downData{jj}.x = interp1(tracktime,x,t,'linear');
        downData{jj}.y = interp1(tracktime,y,t,'linear');
        downData{jj}.speed = interp1(tracktime,speed,t,'linear'); % mm/s
        downData{jj}.r = interp1(tracktime,r,t,'linear'); % mm/s
        downData{jj}.time = t;
        r_all = [downr_all; downData{jj}.r];
    end
    clear phi x y speed r t

    % Initialise cell arrays
    tsG = cell(Nrois,1);        dtsG = cell(Nrois,1); 
    df_f = cell(Nrois,1);       ddf_f = cell(Nrois,1); 
    spikes = cell(Nrois,1);
    trackData.x = cell(Nrois,1);      trackData.y = cell(Nrois,1);
    trackData.phi = cell(Nrois,1);    trackData.r = cell(Nrois,1); 
    trackData.time = cell(Nrois,1);   trackData.speed = cell(Nrois,1);

    for ii = 1:Nrois
        for jj = 1:Nfiles
            k = assignments(ii,jj);
            if ~isnan(k)
                tsG{ii} = [tsG{ii}; ts{jj}.tsG(k,:)'];
                dtsG{ii} = [dtsG{ii}; ts{jj}.dtsG(k,:)'];
                df_f{ii} = [df_f{ii}; ts{jj}.df_f(k,:)'];
                ddf_f{ii} = [ddf_f{ii}; ts{jj}.ddf_f(k,:)'];
                spikes{ii} = [spikes{ii}; ts{jj}.spikes(k,:)'];

                trackData.phi{ii} = [trackData.phi{ii}; downData{jj}.phi];
                trackData.x{ii} = [trackData.x{ii}; downData{jj}.x];
                trackData.y{ii} = [trackData.y{ii}; downData{jj}.y];
                trackData.speed{ii} = [trackData.speed{ii}; downData{jj}.speed];
                trackData.r{ii} = [trackData.r{ii}; downData{jj}.r];
                trackData.r_all = r_all;
%                     if ~isempty(trackData.time{ii})
%                         tt = trackData.time{ii};
%                     else
%                         tt = 0;
%                     end
                trackData.time{ii} = [trackData.time{ii}; downData{jj}.time];
            end
        end
    end

    % save timeseries and track data
    spike_track_data.tsG = tsG;
    spike_track_data.dtsG = dtsG;
    spike_track_data.df_f = df_f;
    spike_track_data.ddf_f = ddf_f;
    spike_track_data.spikes = spikes;

    spike_track_data.trackData = trackData;

    % save(fname, '-struct', 'spike_track_data')
else
    fprintf('%s: loading timeseries and tracking data\n',[mouseid '_' expname]);
    load(fname)
end

%% Calculate place field maps

% settings
if any(trackData.r_all < 100)
    params.mode_dim = '2D';         % open field
    params.PFmap.Nbins = [16, 16];  % number of location bins in [x y]               
else 
    params.mode_dim = '1D';         % circular linear track
    params.PFmap.Nbins = 30;        % number of location bins               
end
params.PFmap.Nepochs = 1;
params.PFmap.Vthr = 20;
params.PFmap.histsmoothFac = 7;
params.fr = 30.9;

Nepochs = params.PFmap.Nepochs;
fname = [fdir mouseid '_' expname '_ref_' files(1,:) '_PFmap_output.mat'];
if ~exist(fname,'file') || force
    fprintf('%s: generating PFmaps\n', [mouseid '_' expname]);
    if strcmpi(params.mode_dim,'1D')
        % Generate place field maps
        [ occMap, hist, asd, activeData ] = generatePFmap_1d_multisession( spikes, trackData, params );

        % If 1D, sort place field maps 
        [ hist.sort_pfMap, hist.sortIdx ] = sortPFmap_1d( hist.pfMap, hist.infoMap, Nepochs );
        [ asd.sort_pfMap, asd.sortIdx ] = sortPFmap_1d( asd.pfMap, asd.infoMap, Nepochs );
        for en = 1:Nepochs
            hist.sort_pfMap_sm(:,:,en) = hist.pfMap_sm(hist.sortIdx(:,en),:,en);
            hist.sort_normpfMap(:,:,en) = hist.normpfMap(hist.sortIdx(:,en),:,en);
            hist.sort_normpfMap_sm(:,:,en) = hist.normpfMap_sm(hist.sortIdx(:,en),:,en);
            if numel(asd.pcIdx) ~=0
                asd.sort_pfMap(:,:,en) = asd.pfMap(asd.sortIdx(:,en),:,en);
                asd.sort_normpfMap(:,:,en) = asd.normpfMap(asd.sortIdx(:,en),:,en);
            else
                asd.sort_pfMap(:,:,en) = zeros(size(asd.pfMap(:,:,en)));
                asd.sort_normpfMap(:,:,en) = zeros(size(asd.pfMap(:,:,en)));
            end
        end

        % Make plots
        plotPF_1d(occMap, hist, asd, activeData.normspkMap_pertrial, activeData.ytick_files);

        % Save output
        output.occMap = occMap;
        output.hist = hist;
        output.asd = asd;
        output.activeData = activeData;
        output.params = params.PFmap;
        % save(fname,'-struct','output');
    else % '2D'
        [occMap, spkMap, spkIdx, hist, asd, ~, activeData] = generatePFmap_2d(spikes, [], trackData, params, false);

         % Make plots
        plotPF_2d(spkMap, activeData, hist, asd);

        % Save output
        output.occMap = occMap;
        output.spkMap = spkMap;
        output.spkIdx = spkIdx;
        output.hist = hist;
        output.asd = asd;
        output.activeData = activeData;
        output.params = params.PFmap;
        save(fname,'-struct','output');
    end
end

t = toc;
str = sprintf('%s: Processing done in %g hrs\n', [mouseid '_' expname], round(t/3600,2));
cprintf(str)
end



