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
% The section labeled "USER-DEFINED INPUT" requires user input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Matlab version requirement: neuroSEE_neuropilDecon requires at least Matlab R2018


%% SETTINGS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USER-DEFINED INPUT                         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Basic settings
default = true;             % flag to use default parameters
                            % flag to force
force = [false;...              % (1) motion correction even if motion corrected images exist
         false;...              % (2) roi segmentation
         false;...              % (3) neuropil decontamination
         false;...              % (4) spike extraction
         false;...              % (5) tracking data extraction
         true;...              % (6) roi registration across images
         false;...              % (7) spike/tracking data consolidation across images
         false];                % (8) place field mapping

mcorr_method = 'normcorre-nr';    % values: [normcorre, normcorre-r, normcorre-nr, fftRigid] 
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
            [data_locn,comp,err] = load_neuroSEEmodules(true);
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
                'd1', 512,...
                'd2', 512,...
                'max_shift',20,...          % default: 50
                'bin_width',200,...         % default: 200
                'us_fac',50,...             % default: 50
                'init_batch',200);          % default: 200
        end
        % NoRMCorre-nonrigid
        if or(strcmpi(mcorr_method,'normcorre'),strcmpi(mcorr_method,'normcorre-nr') )    
            params.nonrigid = NoRMCorreSetParms(...
                'd1', 512,...
                'd2', 512,...
                'grid_size',[64,64],...     % default: [64,64]
                'overlap_pre',[64,64],...   % default: [64,64]
                'overlap_post',[64,64],...  % default: [64,64]
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
        params.ROIreg.maxthr = [];                     
        params.ROIreg.dist_maxthr = 0.1;             % threshold for turning spatial components into binary masks [default: 0.1]
        params.ROIreg.dist_exp = 0.8;                % power n for distance between masked components: dist = 1 - (and(m1,m2)/or(m1,m2))^n [default: 1]
        params.ROIreg.dist_thr = 0.7;                % threshold for setting a distance to infinity    [default: 0.5]
        params.ROIreg.dist_overlap_thr = 0.8;        % overlap threshold for detecting if one ROI is a subset of another [default: 0.8]
        params.ROIreg.plot_reg = true;
        params.print_msg = true;
        params.nonrigid.print_msg = false;
    % PF mapping
        params.PFmap.Nepochs = 1;             % number of epochs for each 4 min video           [default: 1]
        params.PFmap.histsmoothFac = 7;       % Gaussian smoothing window for histogram extraction        [default: 7]
        params.PFmap.Vthr = 20;               % speed threshold (mm/s) Note: David Dupret uses 20 mm/s    [default: 20]
                                              %                              Neurotar uses 8 mm/s
        params.PFmap.prctile_thr = 95;        % percentile threshold for filtering nonPCs
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

for jj = 1:Nfiles
    file = files(jj,:);
    
    % Check for existing processed data
    check_file = checkforExistingProcData(data_locn, file, params);
    
    %% (1) Motion correction and dezippering
    % Load original image files if forced to do motion correction or if motion corrected files don't exist 
    if force(1) || ~check_file(1)
        [imG,imR] = load_imagefile( data_locn, file );
    else
        imG = []; imR = [];
        t = load([data_locn 'Data/' file(1:8) '/Processed/' file '/mcorr_' mcorr_method '/' file '_mcorr_output.mat']);
        templates{jj} = t.template;
    end
    
    % motion correction
    if any([ any(force(1:2)), ~check_file(2), ~check_file(1) ]) 
        if ~isempty(imG)
            if or(strcmpi(mcorr_method,'normcorre'),strcmpi(mcorr_method,'normcorre-r'))
                params.rigid.d1 = size(imG,1);
                params.rigid.d2 = size(imG,2);
                params.rigid.grid_size = [params.rigid.d1, params.rigid.d2, 1];
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
    release = version('-release'); % Find out what Matlab release version is running
    MatlabVer = str2double(release(1:4));

    if force(2) || ~check_file(2) 
        if strcmpi(segment_method,'CaImAn')
            if MatlabVer > 2017
                beep
                err = sprintf('%s: Lower Matlab version required; skipping ROI segmentation and the rest of processing steps.\n', file);
                cprintf('Errors',err);    
                t = toc;
                str = sprintf('%s: Processing done in %g hrs\n', file, round(t/3600,2));
                cprintf(str)
                return
            end

            clear imR; imR = [];
        end
    end

    [tsG{jj}, df_f{jj}, masks{jj}, ~, params] = neuroSEE_segment( imG, mean(imR,3), data_locn, file, params, force(2) );
    clear imG imR
    
    %% Continue with next steps if Matlab version is at least R2018
    if MatlabVer < 2018
        beep
        err = sprintf('%s: Higher Matlab version required; skipping FISSA and the rest of processing steps.\n', file);
        cprintf('Errors',err);
        t = toc;
        str = sprintf('%s: Processing done in %g hrs\n', file, round(t/3600,2));
        cprintf(str)
        return
    end
    
    %% (3) Neuropil decontamination and timeseries extraction
    if dofissa
        [tsG{jj}, dtsG{jj}, ddf_f{jj}, params] = neuroSEE_neuropilDecon( masks{jj}, data_locn, file, params, force(3) );
        str_fissa = 'FISSA';
    else
        dtsG{jj} = [];
        ddf_f{jj} = [];
        str_fissa = 'noFISSA';
    end

    %% (4) Spike estimation
    [spikes{jj}, params, fname_mat] = neuroSEE_extractSpikes( df_f{jj}, ddf_f{jj}, data_locn, file, params, force(4) );

    % Save spike output if required
    if any([ ~check_file(4), any(force([1,2,4])), and(force(3), dofissa) ])
        spike_output.spikes = spikes{jj};
        spike_output.tsG = tsG{jj};
        spike_output.df_f = df_f{jj};
        if ~isempty(dtsG), spike_output.dtsG = dtsG{jj}; end
        if ~isempty(ddf_f), spike_output.ddf_f = ddf_f{jj}; end
        spike_output.params = params.spkExtract;
        save(fname_mat,'-struct','spike_output');
    end

    if manually_refine_spikes
        GUI_manually_refine_spikes( spikes{jj}, tsG{jj}, dtsG{jj}, df_f{jj}, ddf_f{jj}, data_locn, file, params, corr_image, masks );
        uiwait 
    end

    %% (5) Tracking data extraction
    fname_track = findMatchingTrackingFile( data_locn, file, force(5) );
    trackData{jj} = load_trackfile( data_locn,file, fname_track, force(5) );

end


%% (6) ROI registration across sessions
reffile = files(1,:);
fdir = [data_locn 'Analysis/' mouseid '/environment_PFmaps/' mouseid '_' expname '_ref_' reffile '/'...
        mcorr_method '_' segment_method '_' str_fissa '/'];
    if ~exist(fdir,'dir'), mkdir(fdir); end
fname = [fdir  mouseid '_' expname '_ref_' reffile '_registered_rois.mat'];

if ~check_list(1) || force(6)
    for jj = 1:Nfiles
        masks_n = masks{jj};
        A_temp = zeros( size(masks_n,1)*size(masks_n,2), size(masks_n,3 ));
        for ii = 1:size(masks_n,3)
            mask = masks_n(:,:,ii);
            A_temp(:,ii) = mask(:);
            outlines{jj,ii} = bwboundaries(mask);    % boundary of each ROI
        end
        A{jj} = sparse(A_temp); % masks
    end
    
    params.ROIreg.d1 = size(masks{1},1);
    params.ROIreg.d2 = size(masks{1},2);
    params.ROIreg_mc = params.nonrigid;
    params.ROIreg.figsave = true;
    params.ROIreg_mc.grid_size = [32 32 1];
    params.ROIreg_mc.correct_bidir = false;
    params.ROIreg_mc.print_msg = false;

    fprintf('%s: registering ROIs\n',[mouseid '_' expname]);
    [A_union, assignments, matchings] = register_multisession(A, params.ROIreg, templates, params.ROIreg_mc, fdir, list, reffile);
    Mmasks = reshape(full(A_union), params.ROIreg.d1, params.ROIreg.d2, size(A_union,2));
    Nrois = size(Mmasks,3);

    % save masks
    registered_rois.masks = Mmasks;
    registered_rois.outlines = outlines;
    registered_rois.assignments = assignments;
    registered_rois.matchings = matchings;
    registered_rois.params = params.ROIreg;

    fprintf('%s: saving registered ROIs\n',[mouseid '_' expname]);
    save(fname, '-struct', 'registered_rois')

    do_ROIplots = true;
else
    fprintf('%s: loading registered ROIs\n',[mouseid '_' expname]);
    C = load(fname);
    Mmasks = C.masks;
    outlines = C.outlines;
    assignments = C.assignments;
    params.ROIreg = C.params;
    Nrois = size(Mmasks,3);
    
    if exist([fdir,'registered_rois/'],'dir')
        matfiles = subdir(fullfile(fdir,'registered_rois',['*.','fig']));
        if numel(matfiles) > 0     % superposition plots of matched ROIs
            do_ROIplots = false;
        else
            do_ROIplots = true;
        end
    else
        do_ROIplots = true;
    end
end

if do_ROIplots
    if ~exist([fdir 'registered_rois/'],'dir')
        mkdir([fdir 'registered_rois/']);
    end

    [nRow, nCol] = getnRownCol(Nrois);
    nPlot = nRow*nCol;
    for f = 0:Nrois/nPlot
        fh = figure;
        ha = tight_subplot(nRow,nCol,[.01 .01],[.01 .05],[.01 .01]);
        for ii = 0:nPlot-1
            if (f*nPlot+ii+1) <= Nrois
                axes(ha(ii+1));
                imshow(zeros(512,512)); hold on
                for jj = 1:Nfiles
                    k = assignments(f*nPlot+ii+1,jj);
                    if ~isnan(k)
                        plot(outlines{jj,k}{1}(:,2),outlines{jj,k}{1}(:,1),'Linewidth',1.5); hold on
                    end
                end
                axis off; title(['Cell ' num2str(f*nPlot+ii+1)],'fontsize',15);
            end
        end
        if Nrois/nPlot <= 1
            fname_fig = [fdir 'registered_rois/' mouseid '_' expname '_ref_' files(1,:) '_registered_rois'];
        else
            fname_fig = [fdir 'registered_rois/' mouseid '_' expname '_ref_' files(1,:) '_registered_rois_' num2str(f+1)];
        end
        savefig( fh, fname_fig );
        saveas( fh, fname_fig, 'png' );
        close( fh );
    end 
end
    
%% (7) Spike/tracking data consolidation
fname = [fdir  mouseid '_' expname '_ref_' reffile '_spikes_tracking_data.mat'];

if ~exist(fname,'file') || force(7)
    fprintf('%s: concatinating timeseries and tracking data\n',[mouseid '_' expname]);
    [ downData, r_all ] = downsample_trackData( trackData, spikes, params.fr );    
    clear phi x y speed r t

    % Initialise cell arrays
    MtsG = cell(Nrois,1);               MdtsG = cell(Nrois,1); 
    Mdf_f = cell(Nrois,1);              Mddf_f = cell(Nrois,1); 
    Mspikes = cell(Nrois,1);
    MtrackData.x = cell(Nrois,1);       MtrackData.y = cell(Nrois,1);
    MtrackData.phi = cell(Nrois,1);     MtrackData.r = cell(Nrois,1); 
    MtrackData.time = cell(Nrois,1);    MtrackData.speed = cell(Nrois,1);

    for ii = 1:Nrois
        for jj = 1:Nfiles
            k = assignments(ii,jj);
            if ~isnan(k)
                MtsG{ii} = [MtsG{ii}; tsG{jj}(k,:)'];
                MdtsG{ii} = [MdtsG{ii}; dtsG{jj}(k,:)'];
                Mdf_f{ii} = [Mdf_f{ii}; df_f{jj}(k,:)'];
                Mddf_f{ii} = [Mddf_f{ii}; ddf_f{jj}(k,:)'];
                Mspikes{ii} = [Mspikes{ii}; spikes{jj}(k,:)'];

                MtrackData.phi{ii} = [MtrackData.phi{ii}; downData{jj}.phi];
                MtrackData.x{ii} = [MtrackData.x{ii}; downData{jj}.x];
                MtrackData.y{ii} = [MtrackData.y{ii}; downData{jj}.y];
                MtrackData.speed{ii} = [MtrackData.speed{ii}; downData{jj}.speed];
                MtrackData.r{ii} = [MtrackData.r{ii}; downData{jj}.r];
%                     if ~isempty(trackData.time{ii})
%                         tt = trackData.time{ii};
%                     else
%                         tt = 0;
%                     end
                MtrackData.time{ii} = [MtrackData.time{ii}; downData{jj}.time];
            end
        end
    end

    % save timeseries and track data
    spike_track_data.tsG = MtsG;
    spike_track_data.dtsG = MdtsG;
    spike_track_data.df_f = Mdf_f;
    spike_track_data.ddf_f = Mddf_f;
    spike_track_data.spikes = Mspikes;
    spike_track_data.trackData = MtrackData;

    save(fname, '-struct', 'spike_track_data')
else
    fprintf('%s: loading timeseries and tracking data\n',[mouseid '_' expname]);
    C = load(fname);
    MtrackData = C.trackData;
    Mspikes = C.spikes;
end

%% (8) Place field mapping

% settings
if any(r_all < 100)
    params.mode_dim = '2D';         % open field
    params.PFmap.Nbins = [16, 16];  % number of location bins in [x y]               
else 
    params.mode_dim = '1D';         % circular linear track
    params.PFmap.Nbins = 30;        % number of location bins               
end

Nepochs = params.PFmap.Nepochs;
prctile_thr = 95;                   % percentile threshold for filtering nonPCs
params.PFmap.prctile_thr = prctile_thr;

fname = [fdir mouseid '_' expname '_ref_' reffile '_PFmap_output.mat'];
if ~exist(fname,'file') || force(8)
    fprintf('%s: generating PFmaps\n', [mouseid '_' expname]);
    if strcmpi(params.mode_dim,'1D')
        % Generate place field maps
        [ occMap, hist, asd, activeData ] = generatePFmap_1d_multisession( Mspikes, MtrackData, params );

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

        % Save output
        output.occMap = occMap;
        output.hist = hist;
        output.asd = asd;
        output.activeData = activeData;
        output.params = params.PFmap;
        save(fname,'-struct','output');

    else % '2D'
        [occMap, spkMap, spkIdx, hist, asd, ~, activeData] = generatePFmap_2d(spikes, [], trackData, params, false);

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
            
    do_PFplots = true;

else
    fprintf('%s: loading place field mapping data\n',[mouseid '_' expname]);
    C = load(fname);
    occMap = C.occMap;
    hist = C.hist;
    asd = C.asd;
    activeData = C.activeData;
    
    if exist([fdir,'PFmaps/'],'dir')
        matfiles = subdir(fullfile(fdir,'PFmaps',['*.','fig']));
        if numel(matfiles) > 0     % superposition plots of matched ROIs
            do_PFplots = false;
        else
            do_PFplots = true;
        end
    else
        do_PFplots = true;
    end
end
if do_PFplots
    fprintf('%s: generating place field maps\n',[mouseid '_' expname]);
    if strcmpi(params.mode_dim,'1D')
        plotPF_1d(list, occMap, hist, asd, activeData.normspkMap_pertrial, activeData.ytick_files, true, fdir, files);
    else
        plotPF_2d(spkMap, activeData, hist, asd);
    end
end


t = toc;
str = sprintf('%s: Processing done in %g hrs\n', [mouseid '_' expname], round(t/3600,2));
cprintf(str)
end


