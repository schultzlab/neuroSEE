% Written by Ann Go

function frun_pipeline_imreg( list, reffile, slacknotify )

if nargin<3, slacknotify = false; end
% if nargin<2, see line 121

tic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USER-DEFINED INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Basic settings
test = false;               % flag to use one of smaller files in test folder
default = true;             % flag to use default parameters
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
                'max_shift',20,...          % default: 50
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
                'max_shift',15,...          % default: 50
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
        params.PFmap.prctile_thr = 5;        % percentile threshold for filtering nonPCs       [default: 99]
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

params.methods.mcorr_method = mcorr_method;
params.methods.segment_method = segment_method;
params.methods.dofissa = dofissa;

% Default parameters
if default
    params = load_defaultparams(params);
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

% Send Ann slack message if processing has started
if slacknotify
    slacktext = [mouseid '_' expname ': starting CaImAn'];
    neuroSEE_slackNotify( slacktext );
end


%% Check if list has been processed. If not, continue processing unless forced to overwrite 
% existing processed data
% check(1:3) check for existing data in processing steps 2-6 
% check(4) checks for existing mat file pooling all processed data for list

check = checkforExistingProcData(data_locn, list, params, reffile);

% Some security measures
force = logicalForce(force);    % Only allow combinations of force values that make sense

if ~any(force) && check(4)
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
    
%% Image registration
for n = 1:Nfiles
    file = files(n,:);
    if ~strcmpi( file, reffile )
        imdir = [data_locn 'Data/' file(1:8) '/Processed/' file '/mcorr_' mcorr_method '_ref' reffile '/'];
        
        if force(1) || ~exist([imdir file '_imreg_ref' reffile '_output.mat'],'file')
            imG = read_file([ data_locn 'Data/' file(1:8) '/2P/' file '_2P/' file '_2P_XYT_green.tif']);
            imR = read_file([ data_locn 'Data/' file(1:8) '/2P/' file '_2P/' file '_2P_XYT_red.tif']);
            
            [ imG{n}, ~, ~, ~ ] = neuroSEE_imreg( imG, imR, data_locn, file, reffile, params, force );
        else
            if force(2) || ~check(1)
                fprintf('%s: Found registered image. Loading...\n', [mouseid '_' expname '_' file])
                imG{n} = read_file([ imdir file '_2P_XYT_green_imreg_ref' reffile '.tif' ]);
            end
        end
    else
        if force(2) || ~check(1)
            imdir = [data_locn 'Data/' file(1:8) '/Processed/' file '/mcorr_' mcorr_method '/'];
            fprintf('%s: loading reference image\n', [mouseid '_' expname '_' file])
            imG{n} = read_file([ imdir file '_2P_XYT_green_mcorr.tif' ]);
        end
    end
end

%% ROI segmentation
grp_sdir = [data_locn 'Analysis/' mouseid '/' mouseid '_' expname '/group_proc/'...
        mcorr_method '_' segment_method '_' str_fissa '/'...
        mouseid '_' expname '_imreg_ref' reffile '/'];
    if ~exist(grp_sdir,'dir'), mkdir(grp_sdir); end
grp_sname = [grp_sdir mouseid '_' expname '_ref' reffile '_segment_output.mat'];
    
if force(2) || ~check(1)
    fprintf('%s: downsampling images\n', [mouseid '_' expname])
    for n = 1:Nfiles
        Yii = imG{n};
        Y(:,:,(n-1)*size(Yii,3)+1:n*size(Yii,3)) = Yii;
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

    % ROI segmentation
    cellrad = params.ROIsegment.cellrad;
    maxcells = params.ROIsegment.maxcells;
    if strcmpi(segment_method,'CaImAn')
        [df_f, masks, corr_image] = CaImAn( imG, file, maxcells, cellrad );
        df_f = full(df_f);
        
        % Extract raw timeseries
        [d1,d2,T] = size(imG);
        tsG = zeros(size(df_f));
        for ii = 1:size(masks,3)
            maskind = masks(:,:,ii);
            for t = 1:T
                imG_reshaped = reshape( imG(:,:,t), d1*d2, 1);
                tsG(ii,t) = mean( imG_reshaped(maskind) );
            end
        end
        clear imG
        
        fprintf( '%s: ROI segmentation done. Saving output...\n', [mouseid '_' expname] );
    end

    % Output
    % ROIs overlayed on correlation image
    plotopts.plot_ids = 1; % set to 1 to view the ID number of the ROIs on the plot
    fig = plotContoursOnSummaryImage(corr_image, masks, plotopts);
    
    % save masks.mat in individual file folders
    for n = 1:Nfiles
        file = files(n,:);
        if strcmpi(file,reffile)
            sdir = [data_locn 'Data/' file(1:8) '/Processed/' file '/mcorr_' mcorr_method '/' segment_method '_' mouseid '_' expname '/'];
            if ~exist(sdir,'dir'), mkdir(sdir); end
            sname = [sdir mouseid '_' expname '_ref' reffile '_masks.mat'];
        else
            sdir = [data_locn 'Data/' file(1:8) '/Processed/' file '/mcorr_' mcorr_method '_ref' reffile '/' segment_method '_' mouseid '_' expname '/'];
            if ~exist(sdir,'dir'), mkdir(sdir); end
            sname = [sdir mouseid '_' expname '_ref' reffile '_masks.mat'];
        end
        save(sname,'masks','corr_image');
        savefig(fig,[sdir mouseid '_' expname '_ref' reffile '_ROIs']);
        saveas(fig,[sdir mouseid '_' expname '_ref' reffile '_ROIs'],'png');
    end
    
    % save segmentation output in group analysis folder
    output.tsG = tsG;
    output.df_f = df_f;
    output.masks = masks;
    output.corr_image = corr_image;
    output.params = params.ROIsegment;
    save(grp_sname,'-struct','output');
        
    savefig(fig,[grp_sdir mouseid '_' expname '_ref' reffile '_ROIs']);
    saveas(fig,[grp_sdir mouseid '_' expname '_ref' reffile '_ROIs'],'png');
    close(fig);
else
    if any([any(force), ~check(2), ~check(4)])
        fprintf('%s: loading ROI segmentation results\n',[mouseid '_' expname])
        masks = load(grp_sname,'masks');
        masks = masks.masks;
        corr_image = load(grp_sname,'corr_image');
        tsG = load(grp_sname,'tsG');
        df_f = load(grp_sname,'df_f');
    end
end

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
    
%% FISSA, spike extraction, tracking data loading, PF mapping
grp_sname = [grp_sdir mouseid '_' expname '_ref' reffile '_fissa_spike_track_data.mat'];

if any(force(3:5)) || ~check(2)
    % initialise matrices
    SdtsG = []; Sddf_f = []; Sspikes = [];
    SdownData.phi = [];
    SdownData.x = [];
    SdownData.y = [];
    SdownData.speed = [];
    SdownData.r = [];
    SdownData.time = [];

    for n = 1:Nfiles
        file = files(n,:);
        if strcmpi(file,reffile)
            tiffile = [data_locn 'Data/' file(1:8) '/Processed/' file '/mcorr_' mcorr_method '/' file '_2P_XYT_green_mcorr.tif'];
            fissadir = [data_locn 'Data/' file(1:8) '/Processed/' file '/mcorr_' mcorr_method '/' segment_method '_' mouseid '_' expname '/FISSA/'];
        else
            tiffile = [data_locn 'Data/' file(1:8) '/Processed/' file '/mcorr_' mcorr_method '_ref' reffile '/' file '_2P_XYT_green_imreg_ref' reffile '.tif'];
            fissadir = [data_locn 'Data/' file(1:8) '/Processed/' file '/mcorr_' mcorr_method '_ref' reffile '/' segment_method '_' mouseid '_' expname '/FISSA/'];
        end
        if ~exist( fissadir, 'dir' ), mkdir( fissadir ); end

        %% FISSA
        if force(3) || ~exist([fissadir file '_' mouseid '_' expname '_ref' reffile '_fissa_output.mat'],'file')
            fprintf('%s: doing fissa\n',[mouseid '_' expname '_' file]);
            fname_mat_temp = [fissadir 'FISSAout/matlab.mat'];
            if force(3) || and(~exist([fissadir file '_' mouseid '_' expname '_ref' reffile '_fissa_output.mat'],'file'), ~exist(fname_mat_temp,'file'))
                runFISSA( masks, tiffile, fissadir );
            end

            result = load(fname_mat_temp,'result');

            % Convert decontaminated timeseries cell array structure to a matrix
            dtsG = zeros(size(masks,3),size(result.result.cell0.trial0,2));
            for ii = 1:numel(fieldnames(result.result))
                dtsG(ii,:) = result.result.(['cell' num2str(ii-1)]).trial0(1,:);
            end

            % Calculate df_f
            ddf_f = zeros(size(dtsG));
            xddf_prctile = params.fissa.ddf_prctile;
            for ii = 1:size(dtsG,1)
                x = lowpass( medfilt1(dtsG(ii,:),params.fissa.ddf_medfilt1), 1, params.fr );
                fo = ones(size(x)) * prctile(x,xddf_prctile);
                while fo == 0
                    fo = ones(size(x)) * prctile(x,xddf_prctile+5);
                    xddf_prctile = xddf_prctile+5;
                end
                ddf_f(ii,:) = (x - fo) ./ fo;
            end

            fname_fig = [fissadir file '_' mouseid '_' expname '_ref' reffile];
            plotfissa(dtsG, ddf_f, fname_fig);
            save([fissadir file '_' mouseid '_' expname '_ref' reffile '_fissa_output.mat'],'dtsG','ddf_f');
        else
            fprintf('%s: loading fissa output\n',[mouseid '_' expname '_' file]);
            M = load([fissadir file '_' mouseid '_' expname '_ref' reffile '_fissa_output.mat']);
            dtsG = M.dtsG;
            ddf_f = M.ddf_f;
            if ~exist([fissadir file '_' mouseid '_' expname '_ref' reffile '_fissa_result.fig'],'file')
                fname_fig = [fissadir file '_' mouseid '_' expname '_ref' reffile];
                plotfissa(dtsG, ddf_f, fname_fig);
            end
        end

        %% Spike extraction
        if force(4) || ~exist([fissadir file '_' mouseid '_' expname '_ref' reffile '_spikes.mat'],'file')
            fprintf('%s: extracting spikes\n', [mouseid '_' expname '_' file]);
            C = ddf_f;
            N = size(C,1); T = size(C,2);

            for ii = 1:N
                fo = ones(1,T) * prctile(C(ii,:),params.spkExtract.bl_prctile);
                C(ii,:) = (C(ii,:) - fo); % ./ fo;
            end
            spikes = zeros(N,T);
            for ii = 1:N
                spkmin = params.spkExtract.spk_SNR*GetSn(C(ii,:));
                lam = choose_lambda(exp(-1/(params.fr*params.spkExtract.decay_time)),GetSn(C(ii,:)),params.spkExtract.lam_pr);

                [~,spk,~] = deconvolveCa(C(ii,:),'ar2','method','thresholded','optimize_pars',true,'maxIter',20,...
                                        'window',150,'lambda',lam,'smin',spkmin);
                spikes(ii,:) = spk(:);
            end

            fname_fig = [fissadir file '_' mouseid '_' expname '_ref' reffile];
            plotspikes(spikes, fname_fig);
            save([fissadir file '_' mouseid '_' expname '_ref' reffile '_spikes.mat'],'spikes');
        else
            fprintf('%s: loading spike data\n', [mouseid '_' expname '_' file]);
            M = load([fissadir file '_' mouseid '_' expname '_ref' reffile '_spikes.mat']);
            spikes = M.spikes;
            if ~exist([fissadir file '_' mouseid '_' expname '_ref' reffile '_spikes.fig'],'file')
                plotspikes(spikes, [fissadir file '_' mouseid '_' expname '_ref' reffile]);
            end
        end

        %% Tracking data
        fprintf('%s: loading tracking data\n', [mouseid '_' expname '_' file]);
        if force(5) || ~exist([fissadir file '_' mouseid '_' expname '_ref' reffile '_downData.mat'],'file')
            trackfile = findMatchingTrackingFile(data_locn, file, 0);
            c = load_trackfile(data_locn, files(n,:), trackfile, 0);
            downData = downsample_trackData( c, spikes, params.fr );

            save([fissadir file '_' mouseid '_' expname '_ref' reffile '_downData.mat'],'downData');
        else
            M = load([fissadir file '_' mouseid '_' expname '_ref' reffile '_downData.mat']);
            downData = M.downData;
        end

        % Superset arrays
        SdtsG = [SdtsG dtsG];
        Sddf_f = [Sddf_f ddf_f];
        Sspikes = [Sspikes spikes];

        SdownData.phi = [SdownData.phi; downData.phi];
        SdownData.x = [SdownData.x; downData.x];
        SdownData.y = [SdownData.y; downData.y];
        SdownData.speed = [SdownData.speed; downData.speed];
        SdownData.r = [SdownData.r; downData.r];
        SdownData.time = [SdownData.time; downData.time];
    end


    %% Plot and save superset arrays
    % raw timeseries & dF/F
    if ~exist([grp_sdir mouseid '_' expname '_ref' reffile '_fissa_result.fig'],'file') ||...
       ~exist([grp_sdir mouseid '_' expname '_ref' reffile '_fissa_df_f.fig'],'file') 
        fname_fig = [grp_sdir mouseid '_' expname '_ref' reffile];
        plotfissa(SdtsG, Sddf_f, fname_fig);
    end
    clear dtsG ddf_f

    % spikes
    % NOTE: plotting spikes results in fatal segmentation fault (core dumped)
    % if ~exist([sdir expname str_env '_ref' reffile '_spikes.fig'],'file')
    %     fig = figure;
    %     iosr.figures.multiwaveplot(1:size(Sspikes,2),1:size(Sspikes,1),Sspikes,'gain',5); yticks([]); xticks([]);
    %     title('Spikes','Fontweight','normal','Fontsize',12);
    %     savefig(fig,[sdir expname str_env '_ref' reffile '_spikes']);
    %     saveas(fig,[sdir expname str_env '_ref' reffile '_spikes'],'png');
    %     close(fig);
    % end
    clear spikes

    % Concatenated data
    dtsG = SdtsG;
    ddf_f = Sddf_f;
    spikes = Sspikes;
    trackData = SdownData;

    fprintf('%s: saving fissa, spike, track data', [mouseid '_' expname]);
    save(grp_sname,'dtsG','ddf_f','spikes','trackData');

else
    if force(6) || ~check(4)
        fprintf('%s: loading fissa, spike, track data\n', [mouseid '_' expname]);
        c = load(grp_sname);
        dtsG = c.dtsG;
        ddf_f = c.ddf_f;
        spikes = c.spikes;
        trackData = c.trackData;
    end
end

%% PFmapping
grp_sname = [grp_sdir mouseid '_' expname '_ref' reffile '_PFmap_output.mat'];

if any(trackData.r < 100)
    params.mode_dim = '2D';         % open field
    params.PFmap.Nbins = [16, 16];  % number of location bins in [x y]               
else 
    params.mode_dim = '1D';         % circular linear track
    params.PFmap.Nbins = 67;        % number of location bins       
end

Nepochs = params.PFmap.Nepochs;
if force(6) || ~check(3)
    fprintf('%s: generating PFmaps\n', [mouseid '_' expname]);
    if strcmpi(params.mode_dim,'1D')
        % Generate place field maps
        [ hist, asd, activeData, PFdata ] = generatePFmap_1D( spikes, trackData, params, false );
        
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
        save([grp_sdir mouseid '_' expname '_ref' reffile '_PFmap_output.mat'],'-struct','output');
    end
else
    if ~check(4)
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
