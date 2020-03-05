% Written by Ann Go

% This script runs the complete data processing pipeline for a group of
% image files that have been/have to be REGISTERED (e.g. multiple files for
% one environment).
% Processing steps include:
% (1) image registration for EACH file
% The images are then concatenated.
% For the consolidated GROUP image data, processing steps include
% (2) roi segmentation 
% (3) neuropil decontamination and timeseries extraction
% (4) spike extraction
% (5) tracking data extraction
% (6) place field mapping
%
% INPUTS
% list   : name of text file containing filenames of files to be compared.
%           Typically in the format 'list_m##_expname.txt'.
%   e.g.    'list_m62_fam1nov-fam1.txt'         - all fam1 files in fam1nov experiment
%           'list_m62_fam1nov.txt'              - all files in fam1nov experiments
%           'list_m79_fam1_s1-5.txt'            - all fam1 files across 5 sessions           
%           'list_m86_open_s1-2.txt'            - all open field files across 2 sessions
% reffile   : file to be used as registration template. This file is
%               usually part of 'list' but does not have to be. This file
%               must have already been motion corrected.
% slacknotify : flag to send Slack notification when processing is started
%               or has ended
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The section labeled "USER-DEFINED INPUT" requires user input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Matlab version requirement for running on hpc: 
%   normcorre and CaImAn only work with Matlab R2017a. 
%   FISSA requires at least Matlab R2018

function frun_pipeline_imreg( list, reffile, slacknotify )

if nargin<3, slacknotify = false; end
% if nargin<2, see line 121

tic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USER-DEFINED INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Basic settings
test = false;               % flag to use one of smaller files in test folder
default = false;             % flag to use default parameters
                            % flag to force
force = [false;...              % (1) image registration even if registered images exist
         false;...              % (2) roi segmentation
         false;...              % (3) neuropil decontamination
         false;...              % (4) spike extraction
         false;...              % (5) tracking data consolidation
         true];                % (6) place field mapping
mcorr_method = 'normcorre-nr';  % values: [normcorre, normcorre-r, normcorre-nr, fftRigid] 
                                    % CaImAn NoRMCorre method: 
                                    %   normcorre (rigid + nonrigid) 
                                    %   normcorre-r (rigid),
                                    %   normcorre-nr (nonrigid), 
                                    % fft-rigid method (Katie's)
segment_method = 'CaImAn';      % [ABLE,CaImAn]    
dofissa = true;                 % flag to implement FISSA (when false, overrides force(3) setting)
    if dofissa, str_fissa = 'FISSA'; else, str_fissa = 'noFISSA'; end
    
            % Not user-defined
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

% Processing parameters (if not using default)
if ~default
    params.fr = 30.9;                                % imaging frame rate [default: 30.9]
    % motion correction
        params.mcorr.refChannel = 'green';           % reference channel for motion correction [default: 'green']
        % Katie's method
        if strcmpi(mcorr_method,'fftRigid')
            params.mcorr.fftRigid.imscale = 1;             % image downsampling factor                                             [default: 1]
            params.mcorr.fftRigid.Nimg_ave = 10;           % no. of images to be averaged for calculating pixel shift (zippering)  [default: 10]
            params.mcorr.fftRigid.redoT = 300;             % no. of frames at start of file to redo motion correction for after 1st iteration [default: 300]
        end
        % NoRMCorre-rigid
        if or(strcmpi(mcorr_method,'normcorre'),strcmpi(mcorr_method,'normcorre-r'))
            params.mcorr.normcorre_r = NoRMCorreSetParms(...
                'd1', 512,...
                'd2', 512,...
                'max_shift',30,...          % default: 30
                'bin_width',200,...         % default: 200
                'us_fac',50,...             % default: 50
                'init_batch',200);          % default: 200
            params.mcorr.normcorre_r.print_msg = false;   % default: false
        end
        % NoRMCorre-nonrigid
        if or(strcmpi(mcorr_method,'normcorre'),strcmpi(mcorr_method,'normcorre-nr') )    
            params.mcorr.normcorre_nr = NoRMCorreSetParms(...
                'd1', 512,...
                'd2', 512,...
                'grid_size',[64,64],...     % default: [64,64]
                'overlap_pre',[64,64],...   % default: [64,64]
                'overlap_post',[64,64],...  % default: [64,64]
                'iter',1,...                % default: 1
                'use_parallel',false,...    % default: false
                'max_shift',20,...          % default: 20
                'mot_uf',4,...              % default: 4
                'bin_width',200,...         % default: 200
                'max_dev',3,...             % default: 3
                'us_fac',50,...             % default: 50
                'init_batch',200);          % default: 200
            params.mcorr.normcorre_r.print_msg = false;    % default: false
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
    % PF mapping
        params.PFmap.Nepochs = 1;             % number of epochs for each 4 min video           [default: 1]
        params.PFmap.histsmoothFac = 7;       % Gaussian smoothing window for histogram extraction        [default: 7]
        params.PFmap.Vthr = 20;               % speed threshold (mm/s) Note: David Dupret uses 20 mm/s    [default: 20]
                                              %                              Neurotar uses 8 mm/s
        params.PFmap.prctile_thr = 99;        % percentile threshold for filtering nonPCs       [default: 99]
        params_PFmap.Nbins_1D = 30;           % no. of position bins in 103-cm linear track
        params_PFmap.Nbins_2D = [16 16];      % position bins in 325-mm diameter open field arena
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

params.methods.mcorr_method = mcorr_method;
params.methods.segment_method = segment_method;
params.methods.dofissa = dofissa;

% Default parameters
if default
    params = load_defaultparams(params);
    params_PFmap.Nbins_1D = 30;
    params_PFmap.Nbins_2D = [16 16];
end

% Mouseid, Experiment name, files
[ mouseid, expname ] = find_mouseIDexpname(list);
listfile = [data_locn 'Digital_Logbook/lists/' list];
files = extractFilenamesFromTxtfile( listfile );
if nargin<2, reffile = files(1,:); end
Nfiles = size(files,1);

% Some auto-defined parameters
if str2double(files(1,1:4)) > 2018
    params.FOV = 490;                     % FOV area = FOV x FOV, FOV in um
    params.ROIsegment.cellrad = 6;        % expected radius of a cell (pixels)    
    params.ROIsegment.maxcells = 300;     % estimated number of cells in FOV      
else
    params.FOV = 330; 
    params.ROIsegment.cellrad = 9;            
    params.ROIsegment.maxcells = 200;       
end
release = version('-release'); % Find out what Matlab release version is running
MatlabVer = str2double(release(1:4));


%% Check if list has been processed. If not, continue processing unless forced to overwrite 
% existing processed data
% check(1) check for existing data in processing step  2
% check(2)                                       steps 3-5
% check(3)                                       step  6
% check(4) checks for existing mat file pooling all processed data for list

check_list = checkforExistingProcData(data_locn, list, params, reffile);

% Some security measures
force = logicalForce(force);    % Only allow combinations of force values that make sense

if ~any(force) && check_list(4)
    fprintf('%s: list already processed\n', list)
    return
end


%% Continue with image registration and ROI segmentation if Matlab version is R2017
if strcmpi(comp,'hpc') && MatlabVer > 2017
    beep
    err = sprintf('%s: Lower Matlab version required; skipping processing.\n', [mouseid '_' expname]);
    cprintf('Errors',err);
    return
end

% Send Ann slack message if processing has started
if slacknotify
    slacktext = [mouseid '_' expname ': processing started'];
    neuroSEE_slackNotify( slacktext );
end


%% Location of processed group data for list
grp_sdir = [data_locn 'Analysis/' mouseid '/' mouseid '_' expname '/group_proc/'...
            mcorr_method '_' segment_method '_' str_fissa '/'...
            mouseid '_' expname '_imreg_ref' reffile '/'];
    if ~exist(grp_sdir,'dir'), mkdir(dir); end


%% Image registration
% Load images and do registration if forced to do so or if ROI segmentation data doesn't exist 

if any([ force(1:2), ~check_list(1) ])
    for n = 1:Nfiles
        file = files(n,:);

        if ~strcmpi( file, reffile )
            % Check if file has been registered to reffile.
            check_file = checkfor_mcorrIm( data_locn, file, mcorr_method, reffile );

            if any([ force(1:2), ~check_file, ~check_list(1) ])   
                [fileG,fileR] = load_imagefile( data_locn, file );

                % Send Ann slack message
                if slacknotify
                    if array_id == 1
                        slacktext = [mouseid '_' expname ': registering 1 of ' num2str(size(files,1)) 'files'];
                        neuroSEE_slackNotify( slacktext );
                    end
                end
                
                if strcmpi(segment_method,'CaImAn') % CaImAn does not use imR
                    [ imG{n}, ~, params.mcorr ] = neuroSEE_motionCorrect( fileG, fileR, data_locn, file, ...
                                                            mcorr_method, params.mcorr, reffile, force(1) );
                    imR = [];
                else
                    [ imG{n}, ~, params.mcorr, imR{n} ] = neuroSEE_motionCorrect( fileG, fileR, data_locn, file, ...
                                                            mcorr_method, params.mcorr, reffile, force(1) );
                end
            else 
                fprintf( '%s: Registered images found. Skipping registration\n', [mouseid '_' expname '_' file] );
                imG = []; imR = [];
            end
        else
            if force(2) || ~check_list(1)
                imdir = [data_locn 'Data/' file(1:8) '/Processed/' file '/mcorr_' mcorr_method '/'];
                fprintf('%s: Loading reference image\n', [mouseid '_' expname '_' file]);
                imG{n} = read_file([ imdir file '_2P_XYT_green_mcorr.tif' ]);
                if strcmpi(segment_method,'CaImAn') % CaImAn does not use imR
                    imR = [];
                else
                    imR{n} = read_file([ imdir file '_2P_XYT_red_mcorr.tif' ]);
                end
                
                fname_mat = [grp_sdir mouseid '_' expname '_ref' reffile '_mcorr_template.mat'];
                fname_fig = [grp_sdir mouseid '_' expname '_ref' reffile '_imreg_template.fig'];
                if ~exist(fname_mat, 'file') || ~exist(fname_fig, 'file')
                    c = load([imdir reffile '_mcorr_output.mat']);
                    template_g = c.green.meanregframe;
                    template_r = c.red.meanregframe;
                    fig = figure; 
                    subplot(121); imagesc(template_g); title('GCaMP6');
                    subplot(122); imagesc(template_r); title('mRuby');
                    save(fname_mat,'template_g','template_r')
                    savefig(fig, fname_fig(1:end-4));
                    saveas(fig, fname_fig(1:end-4),'png');
                    close(fig);
                end
            end
        end
    end
end


%% ROI segmentation    
    fprintf('%s: downsampling images\n', [mouseid '_' expname])
    for n = 1:Nfiles
        Yii = imG{n};
        Y(:,:,(n-1)*size(Yii,3)+1:n*size(Yii,3)) = Yii;
    end

    if ~strcmpi(segment_method,'CaImAn') % CaImAn does not use imR
        for n = 1:Nfiles
            Xii = imR{n};
            X(:,:,(n-1)*size(Xii,3)+1:n*size(Xii,3)) = Xii;
        end
        clear Xii
    end
    
    % Downsample images
    if size(files,1) <= 7
        tsub = 5;
    elseif size(files,1) <= 10
        tsub = 7;
    elseif size(files,1) <= 13
        tsub = 9;
    else
        tsub = round( size(files,1)*7420/11000 );
    end
    imG = Y(:,:,1:tsub:end); % downsampled concatenated image
    clear Y Yii
    if ~strcmpi(segment_method,'CaImAn')
        imR = X(:,:,1:tsub:end);
        clear X
    end

    % ROI segmentation
    prevstr = sprintf( '%s: doing ROI segmentation...\n', [mouseid '_' expname] );
    cprintf('Text',prevstr);
    [~, df_f, masks, corr_image, params] = neuroSEE_segment( imG, data_locn, file, params, force(2), mean(imR,3), list, reffile );
    str = sprintf( '%s: ROI segmentation done.\n', [mouseid '_' expname] );
    refreshdisp( str, prevstr );

    
%% Continue with next steps if Matlab version is at least R2018
if strcmpi(comp,'hpc') && MatlabVer < 2018
    beep
    err = sprintf('%s: Higher Matlab version required; skipping FISSA and the rest of processing steps.\n', [mouseid '_' expname]);
    cprintf('Errors',err);
    t = toc;
    str = sprintf('%s: Processing done in %g hrs\n', file, round(t/3600,2));
    cprintf(str)
    return
end
    

%% FISSA, Spike estimation, Tracking data 
for n = 1:Nfiles
    file = files(n,:);
    if dofissa
        [~, dtsG{n}, ddf_f{n}, params] = neuroSEE_neuropilDecon( masks, data_locn, file, params, force(3), list, reffile );
        [ spikes{n}, params ] = neuroSEE_extractSpikes( [], ddf_f{n}, data_locn, file, params, force(4), list, reffile );
    end

    fprintf('%s: loading tracking data\n', [mouseid '_' expname '_' file]);
    if force(5) || ~exist([data_locn 'Data/' file(1:8) '/Processed/' file '/behaviour/' file '_downTrackdata.mat'],'file')
        trackfile = findMatchingTrackingFile(data_locn, file, force(5));
        Trackdata = load_trackfile(data_locn, files(n,:), trackfile, force(5));
        downTrackdata{n} = downsample_trackData( Trackdata, spikes, params.fr );
        save([data_locn 'Data/' file(1:8) '/Processed/' file '/behaviour/' file '_downTrackdata.mat'],'downData');
    else
        M = load([data_locn 'Data/' file(1:8) '/Processed/' file '/behaviour/' file '_downTrackdata.mat']);
        downTrackdata{n} = M.downTrackdata;
    end
end

if ~dofissa
    [ spikes, params ] = neuroSEE_extractSpikes( df_f, [], data_locn, file, params, force(4), list, reffile );
end

%% Concatenate data
% initialise matrices
if dofissa
    SdtsG = []; Sddf_f = []; 
    Sspikes = [];
end
SdownTrackdata.phi = [];
SdownTrackdata.x = [];
SdownTrackdata.y = [];
SdownTrackdata.speed = [];
SdownTrackdata.r = [];
SdownTrackdata.time = [];

for n = 1:Nfiles
    if dofissa
        SdtsG = [SdtsG dtsG{n}];
        Sddf_f = [Sddf_f ddf_f{n}];
        Sspikes = [Sspikes spikes{n}];
    end
    SdownTrackdata.phi = [SdownTrackdata.phi; downTrackdata{n}.phi];
    SdownTrackdata.x = [SdownTrackdata.x; downTrackdata{n}.x];
    SdownTrackdata.y = [SdownTrackdata.y; downTrackdata{n}.y];
    SdownTrackdata.speed = [SdownTrackdata.speed; downTrackdata{n}.speed];
    SdownTrackdata.r = [SdownTrackdata.r; downTrackdata{n}.r];
    SdownTrackdata.time = [SdownTrackdata.time; downTrackdata{n}.time];
end

% Plot and save superset arrays
% raw timeseries & dF/F
if dofissa
    if ~exist([grp_sdir mouseid '_' expname '_ref' reffile '_fissa_result.fig'],'file') ||...
       ~exist([grp_sdir mouseid '_' expname '_ref' reffile '_fissa_df_f.fig'],'file') 
        fname_fig = [grp_sdir mouseid '_' expname '_ref' reffile];
        plotfissa(SdtsG, Sddf_f, fname_fig);
    end
    clear dtsG ddf_f
end

% spikes
% NOTE: plotting spikes results in fatal segmentation fault (core dumped)
% if ~exist([sdir expname str_env '_ref' reffile '_spikes.fig'],'file')
%     if dofissa
%         plotspikes(Sspikes, fname_fig);
%     else
%         plotspikes(spikes, fname_fig);
%     end
% end
if dofissa
    clear spikes
end

% Concatenated data
if dofissa
    dtsG = SdtsG;
    ddf_f = Sddf_f;
    spikes = Sspikes;
end
downTrackdata = SdownTrackdata;

if dofissa
    fprintf('%s: saving fissa, spike, downsampled track data', [mouseid '_' expname]);
    grp_sname = [grp_sdir mouseid '_' expname '_ref' reffile '_fissa_spike_track_data.mat'];
    save(grp_sname,'dtsG','ddf_f','spikes','downTrackdata');
else
    fprintf('%s: saving raw timeseries, spike, downsampled track data', [mouseid '_' expname]);
    grp_sname = [grp_sdir mouseid '_' expname '_ref' reffile '_raw_df_f_spike_track_data.mat'];
    save(grp_sname,'tsG','df_f','spikes','downTrackdata');
end


%% PFmapping
grp_sname = [grp_sdir mouseid '_' expname '_ref' reffile '_PFmap_output.mat'];

if any(downTrackdata.r < 100)
    params.mode_dim = '2D';         % open field
    params.PFmap.Nbins = [16, 16];  % number of location bins in [x y]               
else 
    params.mode_dim = '1D';         % circular linear track
    params.PFmap.Nbins = 30;        % number of location bins       
end

Nepochs = params.PFmap.Nepochs;
if force(6) || ~check_list(3)
    fprintf('%s: generating PFmaps\n', [mouseid '_' expname]);
    if strcmpi(params.mode_dim,'1D')
        % Generate place field maps
        [ hist, asd, activeData, PFdata ] = generatePFmap_1D( spikes, downTrackdata, params, false );
        
        % If 1D, sort place field maps 
        if isfield(hist,'pfMap_MI')
            [ hist.sort_pfMap_MI, hist.sortIdx_MI ] = sortPFmap_1d( hist.pfMap_MI, hist.infoMap_MI );
            for en = 1:Nepochs
                hist.sortIdx_MI_roiNum(:,en) = hist.pcIdx_MI(hist.sortIdx_MI(:,en));
                hist.sort_pfMap_MI_sm(:,:,en) = hist.pfMap_MI_sm(hist.sortIdx_MI(:,en),:,en);
                hist.sort_normpfMap_MI(:,:,en) = hist.normpfMap_MI(hist.sortIdx_MI(:,en),:,en);
                hist.sort_normpfMap_MI_sm(:,:,en) = hist.normpfMap_MI_sm(hist.sortIdx_MI(:,en),:,en);
            end
        end
        if isfield(hist,'pfMap_SIsec')
            [ hist.sort_pfMap_SIsec, hist.sortIdx_SIsec ] = sortPFmap_1d( hist.pfMap_SIsec, hist.infoMap_SIsec );
            for en = 1:Nepochs
                hist.sortIdx_SIsec_roiNum(:,en) = hist.pcIdx_SIsec(hist.sortIdx_SIsec(:,en));
                hist.sort_pfMap_SIsec_sm(:,:,en) = hist.pfMap_SIsec_sm(hist.sortIdx_SIsec(:,en),:,en);
                hist.sort_normpfMap_SIsec(:,:,en) = hist.normpfMap_SIsec(hist.sortIdx_SIsec(:,en),:,en);
                hist.sort_normpfMap_SIsec_sm(:,:,en) = hist.normpfMap_SIsec_sm(hist.sortIdx_SIsec(:,en),:,en);
            end
        end
        if isfield(hist,'pfMap_SIspk')
            [ hist.sort_pfMap_SIspk, hist.sortIdx_SIspk ] = sortPFmap_1d( hist.pfMap_SIspk, hist.infoMap_SIspk );
            for en = 1:Nepochs
                hist.sortIdx_SIspk_roiNum(:,en) = hist.pcIdx_SIspk(hist.sortIdx_SIspk(:,en));
                hist.sort_pfMap_SI_spk_sm(:,:,en) = hist.pfMap_SIspk_sm(hist.sortIdx_SIspk(:,en),:,en);
                hist.sort_normpfMap_SIspk(:,:,en) = hist.normpfMap_SIspk(hist.sortIdx_SIspk(:,en),:,en);
                hist.sort_normpfMap_SIspk_sm(:,:,en) = hist.normpfMap_SIspk_sm(hist.sortIdx_SIspk(:,en),:,en);
            end
        end
        
        if isfield(asd,'pfMap_MI')
            [ asd.sort_pfMap_MI, asd.sortIdx_MI ] = sortPFmap_1d( asd.pfMap_MI, asd.infoMap_MI );
            for en = 1:Nepochs
                asd.sortIdx_MI_roiNum(:,en) = asd.pcIdx_MI(asd.sortIdx_MI(:,en));
                asd.sort_normpfMap_MI(:,:,en) = asd.normpfMap_MI(asd.sortIdx_MI(:,en),:,en);
            end
        end
        if isfield(asd,'pfMap_SIsec')
            [ asd.sort_pfMap_SIsec, asd.sortIdx_SIsec ] = sortPFmap_1d( asd.pfMap_SIsec, asd.infoMap_SIsec );
            for en = 1:Nepochs
                asd.sortIdx_SIsec_roiNum(:,en) = asd.pcIdx_SIsec(asd.sortIdx_SIsec(:,en));
                asd.sort_normpfMap_SIsec(:,:,en) = asd.normpfMap_SIsec(asd.sortIdx_SIsec(:,en),:,en);
            end
        end
        if isfield(asd,'pfMap_SIspk')
            [ asd.sort_pfMap_SIspk, asd.sortIdx_SIspk ] = sortPFmap_1d( asd.pfMap_SIspk, asd.infoMap_SIspk );
            for en = 1:Nepochs
                asd.sortIdx_SIspk_roiNum(:,en) = asd.pcIdx_SIspk(asd.sortIdx_SIspk(:,en));
                asd.sort_normpfMap_SIspk(:,:,en) = asd.normpfMap_SIspk(asd.sortIdx_SIspk(:,en),:,en);
            end
        end

        % Make plots
        plotPF_1d(hist, asd, PFdata, true, [grp_sdir 'PFdata/'], ...
                  [mouseid '_' expname '_ref' reffile], true)
        
        % Save output
        output.PFdata = PFdata;
        output.hist = hist;
        output.asd = asd;
        output.activeData = activeData;
        output.params = params.PFmap;
        save([grp_sdir mouseid '_' expname '_ref' reffile '_PFmap_output.mat'],'-struct','output');
    
    else % '2D'
        [occMap, spkMap, spkIdx, hist, asd, ~, activeData] = generatePFmap_2d(spikes, [], downTrackdata, params, false);

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
        save([grp_sdir mouseid '_' expname '_ref' reffile '_PFmap_output.mat'],'-struct','output');
    end
else
    if ~check_list(4)
        fprintf('%s: loading PF mapping data\n', [mouseid '_' expname]);
        c = load(grp_sname);
        activeData = c.activeData;
        PFdata = c.PFdata;
        hist = c.hist;
        asd = c.asd;
    end
end

sname_allData = [ grp_sdir mouseid '_' expname '_ref' reffile '_allData.mat'];

fprintf('%s: saving all data\n', [mouseid '_' expname]);
save(sname_allData,'list','corr_image','masks','tsG','df_f','dtsG','ddf_f','spikes',...
                    'trackData','activeData','PFdata','hist','asd','params');
if exist('spkMap','var'), save(fname_allData,'-append','spkMap'); end
if exist('spkIdx','var'), save(fname_allData,'-append','spkIdx'); end


% Send Ann slack message if processing has finished
if slacknotify
    slacktext = [expname ': CaImAn FINISHED. No errors!'];
    neuroSEE_slackNotify( slacktext );
end

t = toc;
str = sprintf('%s: Processing done in %g hrs\n', [mouseid '_' expname], round(t/3600,2));
cprintf(str)

end
