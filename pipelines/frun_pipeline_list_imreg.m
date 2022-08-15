% Written by Ann Go

% This script runs the complete data processing pipeline for a group of
% image files that have been/have to be registered (e.g. multiple files for
% one environment).
% Processing steps include:
% (1) image registration for EACH file
% The images are then concatenated.
% For the consolidated GROUP image data, processing steps include
% (2) roi segmentation 
% (3) neuropil decontamination and timeseries estimation
% (4) spike estimation
% (5) tracking data estimation
% (6) place field mapping
%
% INPUTS
% list   : name of text file containing filenames of files to be compared.
%           Typically in the format 'list_m##_expname.txt'.
%   e.g.    'list_m62_fam1nov-fam1.txt'         - all fam1 files in fam1nov experiment
%           'list_m62_fam1nov.txt'              - all files in fam1nov experiments
%           'list_m79_fam1_s1-5.txt'            - all fam1 files across 5 sessions           
%           'list_m86_open_s1-2.txt'            - all open field files across 2 sessions
% reffile   : (optional, default: first file) file to be used as registration template. This file is
%               usually part of 'list' but does not have to be. This file
%               must have already been motion corrected. If no reffile is
%               specified, first file on the list is used.
% conc_runs  : (optional, default: false) flag if rois were segmented from concatenated files from
%               different environments e.g. fam1fam2fam1-fam1 but rois are
%               for fam1fam2fam1. Do not flag for fam1fam2fam1 since in
%               this case it is understood that the rois are from the
%               concatenated environments. However, if you accidentally do
%               this, the script will override the setting and use the
%               correct value.
% dostep    : (optional, default: [1; 1; 1; 1; 1; 1]) array of 6. Set to 0 or 1 to skip (0) or implement
%               (1) a processing step (see list at the top)
% force     : (optional, default: [0; 0; 0; 0; 0; 0]) array of 6. Set to 1 to force the implementation of
%               a processing step (even if output already exists)
% tsub      : (optional, default: 5) temporal downsampling factor for CaImAn. Typically chosen so that
%                no. of files in list x 7420
%               -----------------------------  = 23,000
%                           tsub                    
%               This reduces out-of-memory errors. Keep tsub value below 10.
% min_SNR   : (optional, default: 3 for 330x330um FOV data, 2.5 for 490x490um FOV data) CaImAn parameter. 
%               Minimum SNR for accepting exceptional events. 
% bl_prctile : (optional, default: 85) Parameter for spike estimation.  Percentile to be
%               used for estimating baseline.
% 
% Matlab version requirement for running on hpc: 
%   normcorre only works up to Matlab R2017b.
%   CaImAn needs at least R2017b
%   FISSA requires at least Matlab R2018


function frun_pipeline_list_imreg( list, reffile, conc_runs, numfiles, dostep, force, tsub, min_SNR, bl_prctile )

if nargin<9, bl_prctile = 85; end
% if nargin<8, see line 97
if nargin<7, tsub = 5; end
if nargin<6, force = [0; 0; 0; 0; 0; 0]; end
if nargin<5, dostep = [1; 1; 1; 1; 1; 1]; end
% if nargin<4, see line 156
if nargin<3, conc_runs = false; end
% if nargin<2, see line 121
tic

% Load module folders and define data directory
[data_locn,~,err] = load_neuroSEEmodules;
if ~isempty(err)
    beep
    cprintf('Errors',err);    
    return
end

% Mouseid, Experiment name, files
[ mouseid, expname, fov ] = find_mouseIDexpname(list);
listfile = [data_locn 'Digital_Logbook/lists_imaging/' list];
try
    files = extractFilenamesFromTxtfile( listfile );
catch
    beep
    cprintf('Errors','List does not exist. Check list name.\n');    
    return
end
if nargin<2, reffile = files(1,:); end
Nfiles = size(files,1);

% Auto-defined parameters
% FOV area = FOV x FOV, FOV in um
if str2double(files(1,1:4)) > 2018
    FOV = 490;                     
else
    FOV = 330; 
end
mousenum = str2double(mouseid(2:end));
% min_SNR
if nargin<8 || isempty(min_SNR)
    if mousenum > 76
        min_SNR = 2.5;
    else
        min_SNR = 3;
    end
end
% virus (which determines decay time constant for calcium transient)
if mousenum > 104
    virus = 'jGCaMP7s';
else
    virus = 'GCaMP6s';
end

groupreg_method = 'imreg';      % method for concatenating file data (either register images or rois)
                                % Here, we register the images
                                
%% SETTINGS
% Basic settings
mcorr_method = 'normcorre';     % image registration method
                                % values: [normcorre, normcorre-r, normcorre-nr, fftRigid] 
                                    % CaImAn NoRMCorre method: 
                                    %   normcorre (rigid + nonrigid) 
                                    %   normcorre-r (rigid),
                                    %   normcorre-nr (nonrigid), 
                                    % fft-rigid method (Katie's)
segment_method = 'CaImAn';      % [ABLE,CaImAn]    
runpatches = false;             % for CaImAn processing, flag to run patches (default: false)
dofissa = false;
doasd = false;                  % flag to do asd pf calculation

% Processing parameters (any parameter that is not set gets a default value)
% Add any parameters you want to set after FOV. See neuroSEE_setparams for
% full set of parameters
params = neuroSEE_setparams(...
            'groupreg_method', groupreg_method,...
            'mcorr_method', mcorr_method,...
            'segment_method', segment_method,...
            'runpatches', runpatches,...
            'dofissa', dofissa,...
            'doasd', doasd,...
            'FOV', FOV,...
            'tsub', tsub,...
            'virus', virus,...
            'min_SNR', min_SNR,...
            'bl_prctile', bl_prctile);         

% For slack notifications
slacknotify = false;
% if true, set below
slackURL = '';
slackTarget = '@xxx';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if dofissa, str_fissa = 'FISSA'; else, str_fissa = 'noFISSA'; end

% Matlab version
release = version('-release'); % Find out what Matlab release version is running
MatlabVer = str2double(release(1:4));

%% Check if list has been processed. If not, continue processing unless forced to overwrite 
% existing processed data
% check(1) check for existing data in processing step (2) roi segmentation 
% check(2)                                       step (3) neuropil decontamination
% check(3)                                       step (4) spike estimation
% check(4)                                       step (5) tracking data consolidation
% check(5)                                       step (6) place field mapping
% check(6) checks for existing mat file pooling all processed data for list

% Some security measures
if ~contains(list,'-')                          % conc_runs = true is only an option for sub experiments
    conc_runs = false;                              % e.g. fam1fam2-fam1, fam1novfam1-nov 
    if nargin<4 || isempty(numfiles)
        fprintf('%s: Number of files per run not provided. Timeseries data will not be distributed to run folders.\n', [mouseid '_' expname]);
    end
end
force = logicalForce(force);        % Only allow combinations of force/dostep values that make sense
dostep = logicaldostep(dostep);     % because later steps require earlier ones

check_list = checkforExistingProcData(data_locn, list, params, reffile, conc_runs);

if ~any(force) && check_list(6)
    fprintf('%s: List already processed\n', list)
    return
end


%% Location of processed group data for list
if ~isempty(fov)
    grp_sdir = [data_locn 'Analysis/' mouseid '/' fov '/' mouseid '_' expname '/group_proc/'...
            groupreg_method '_' mcorr_method '_' segment_method '/'...
            mouseid '_' expname '_imreg_ref' reffile '/'];
else
    grp_sdir = [data_locn 'Analysis/' mouseid '/' mouseid '_' expname '/group_proc/'...
            groupreg_method '_' mcorr_method '_' segment_method '/'...
            mouseid '_' expname '_imreg_ref' reffile '/'];
end
if conc_runs
    grp_sdir = [grp_sdir(1:end-1) '_conrunsrois/'];
end
if ~exist(grp_sdir,'dir'), mkdir(grp_sdir); fileattrib(grp_sdir,'+w','g','s'); end


%% 1) Image registration
% Load images and do registration if forced to do so or if ROI segmentation data doesn't exist 

if dostep(1)
    if force(2) || ~check_list(1)
        % Continue only if Matlab version is R2017
%         if strcmpi(comp,'hpc') && MatlabVer > 2017
%             beep
%             err = sprintf('%s: Lower Matlab version required for motion correction. Cannot proceed.\n', [mouseid '_' expname]);
%             cprintf('Errors',err);
%             return
%         end

        % Send slack notification if processing has started
        if slacknotify
            slacktext = [mouseid '_' expname ': Processing started'];
            SendSlackNotification( slackURL, slacktext, slackTarget );
        end
        
        imG = cell(Nfiles,1);
        if ~strcmpi(segment_method,'CaImAn'), imR = cell(Nfiles,1); end
        for n = 1:Nfiles
            file = files(n,:);

            if ~strcmpi( file, reffile )
                % Check if file has been registered to reffile.
                check_file = checkfor_mcorrIm( data_locn, file, mcorr_method, reffile, false );

                if force(1) || ~check_file   
                    [fileG,fileR] = load_imagefile( data_locn, file );
                else
                    fileG = []; fileR = [];
                end

                if strcmpi(segment_method,'CaImAn') % CaImAn does not use imR
                    [ imG{n}, params.mcorr, ~, template_g, template_r ] = neuroSEE_motionCorrect( fileG, fileR, data_locn, file, ...
                                                            mcorr_method, params.mcorr, reffile, force(1), list, false );
                    imR = [];
                else
                    [ imG{n}, params.mcorr, imR{n}, template_g, template_r ] = neuroSEE_motionCorrect( fileG, fileR, data_locn, file, ...
                                                            mcorr_method, params.mcorr, reffile, force(1), list, false );
                end
            else
                imdir = [data_locn 'Data/' file(1:8) '/Processed/' file '/mcorr_' mcorr_method '/'];
                str = sprintf('%s: Loading reference image\n', [mouseid '_' expname '_' file]);
                cprintf('Text',str)
                imG{n} = read_file([ imdir file '_2P_XYT_green_mcorr.tif' ]);
                if strcmpi(segment_method,'CaImAn') % CaImAn does not use imR
                    imR = [];
                else
                    imR{n} = read_file([ imdir file '_2P_XYT_red_mcorr.tif' ]);
                end
                newstr = sprintf('%s: Reference image loaded\n', [mouseid '_' expname '_' file]);
                refreshdisp(newstr, str)

                fname_mat = [grp_sdir mouseid '_' expname '_ref' reffile '_imreg_template.mat'];
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

        % Image concatenation    
        fprintf('%s: Concatenating images\n', [mouseid '_' expname])
        framesperfile = zeros(Nfiles,1);
        for n = 1:Nfiles
            Yii = imG{n};
            framesperfile(n) = size(Yii,3);
            Y(:,:,sum(framesperfile(1:n-1))+1:sum(framesperfile(1:n))) = Yii;
        end
        clear Yii imG
        imG = Y;
        grp_sname = [grp_sdir mouseid '_' expname '_ref' reffile '_framesperfile.mat'];
        save(grp_sname,'framesperfile');
        
        if ~strcmpi(segment_method,'CaImAn') % CaImAn does not use imR
            for n = 1:Nfiles
                Xii = imR{n};
                X(:,:,(n-1)*size(Xii,3)+1:n*size(Xii,3)) = Xii;
            end
            clear Xii imR;
            imR = X;
        end
    else
        fprintf('%s: Registered images found. Loading templates. Skipping image registration.\n', [mouseid '_' expname]);
        imG = []; imR = [];
        c = load([grp_sdir mouseid '_' expname '_ref' reffile '_imreg_template.mat']);
        template_g = c.template_g;
        template_r = c.template_r;
        m = load([grp_sdir mouseid '_' expname '_ref' reffile '_framesperfile.mat']);
        framesperfile = m.framesperfile;
    end
else
    fprintf('%s: No processing steps ticked. Cannot proceed.\n', [mouseid '_' expname]);
    return
end

%% 2) ROI segmentation
if dostep(2)
    [tsG, df_f, masks, corr_image, params] = neuroSEE_segment( imG, data_locn, params, list, reffile, conc_runs, force(2), mean(imR,3) );
                                             
    % Copy files from grp_sdir to folders for individual runs
    if ~contains(list,'-')
        if nargin<4 || isempty(numfiles)
            fprintf('%s: Number of files per run not provided. ROI segmentation data was not distributed to run folders.\n', [mouseid '_' expname]);
        else
            divide = force(2) || ~check_list(1);
            divide_expdata_into_runs( data_locn, list, reffile, numfiles, [], [1,0,0], [divide,0,0] );
        end
    end
else
    fprintf('%s: ROI segmentation step not ticked. Skipping this and later steps.\n', [mouseid '_' expname]);
    t = toc;
    str = sprintf('%s: Processing done in %g hrs\n', [mouseid '_' expname], round(t/3600,2));
    cprintf(str)
    return
end
    
%% 3) FISSA
if dostep(3)
    dtsG = []; ddf_f = []; 
    if dofissa
        if ~force(3) && check_list(2)
            fname_mat = [grp_sdir '/' str_fissa '/' mouseid '_' expname '_ref' reffile '_fissa_output.mat'];
            s = load(fname_mat);
            dtsG = s.dtsG;
            if size(dtsG,1) ~= size(tsG,1)
                force([3,4,6]) = true; % force FISSA step if no. of ROIs in segmentation and fissa outputs don't match
                dtsG = [];
                fprintf('%s: FISSA output found. Redoing FISSA. No. of ROIs in segmentation and FISSA outputs do not match.\n', [mouseid '_' expname]);
            end
        end
        if force(3) || ~check_list(2)
            [ dtsG, ddf_f, params ] = neuroSEE_neuropilDecon( masks, data_locn, params, list, reffile, conc_runs, force(3) );
        end
        
        % Copy files from grp_sdir to folders for individual runs
        if ~contains(list,'-')
            if nargin<4 || isempty(numfiles)
                fprintf('%s: Number of files per run not provided. FISSA output was not distributed to run folders.\n', [mouseid '_' expname]);
            else
                divide = force(3) || ~check_list(2);
                divide_expdata_into_runs( data_locn, list, reffile, numfiles, [], [1,1,0], [0,divide,0] );
            end
        end
    end
else
    if dofissa
        fprintf('%s: FISSA step not ticked. Skipping this and later steps.\n', [mouseid '_' expname]);
        t = toc;
        str = sprintf('%s: Processing done in %g hrs\n', [mouseid '_' expname], round(t/3600,2));
        cprintf(str)
        return
    end
end

%% 4) Spike estimation
if dostep(4)
    if ~force(4) && check_list(3)
        fname_mat = [grp_sdir '/' str_fissa '/bl_prctile' num2str(bl_prctile) '/' mouseid '_' expname '_ref' reffile '_spikes.mat'];
        s = load(fname_mat);
        spikes = s.spikes;
        if size(spikes,1) ~= size(tsG,1)
            force([4,6]) = true; % force spike estimation step if no. of ROIs in segmentation and spike data don't match
%             clear dtsG
            fprintf('%s: Spike data found. Redoing spike estimation. ROIs in segmentation output and spike data do not match.\n', [mouseid '_' expname]);
        end
    end
    
    [ spikes, bl, neuron_sn, g, params ] = neuroSEE_extractSpikes( df_f, ddf_f, data_locn, params, list, reffile, conc_runs, grp_sdir, force(4) );

    % Copy files from grp_sdir to folders for individual runs
    if ~contains(list,'-')
        if nargin<4 || isempty(numfiles)
            fprintf('%s: Number of files per run not provided. Spike estimation output was not distributed to run folders.\n', [mouseid '_' expname]);
        else
            divide = force(4) || ~check_list(3);
            divide_expdata_into_runs( data_locn, list, reffile, numfiles, [], [1,dofissa,1], [0,0,divide] );
        end
    end
else
    fprintf('%s: Spike estimation step not ticked. Skipping this and later steps.\n', [mouseid '_' expname]);
    t = toc;
    str = sprintf('%s: Processing done in %g hrs\n', [mouseid '_' expname], round(t/3600,2));
    cprintf(str)
    return
end

%% 5) Behaviour tracking data
if dostep(5)
    if ~isempty(fov)
        grp_trackdir = [data_locn 'Analysis/' mouseid '/' fov '/' mouseid '_' expname '/group_proc/'];
    else
        grp_trackdir = [data_locn 'Analysis/' mouseid '/' mouseid '_' expname '/group_proc/'];
    end
    if force(5) || ~check_list(4)
        SdownTrackdata.phi = [];            
        SdownTrackdata.x = [];              
        SdownTrackdata.y = [];              
        SdownTrackdata.speed = [];          
        SdownTrackdata.r = [];              
        SdownTrackdata.time = [];           
        cdownTrackdata = cell(Nfiles,1);    
        
        for n = 1:Nfiles
            file = files(n,:);
            file_sdir = [data_locn 'Data/' file(1:8) '/Processed/' file '/behaviour/'];
            if force(5) || ~exist([file_sdir file '_downTrackdata.mat'],'file')
                trackfile = findMatchingTrackingFile(data_locn, file, force(5));
                c = load_trackfile(data_locn, files(n,:), trackfile, force(5));
                downTrackdata = downsample_trackData( c, framesperfile(n), params.PFmap.fr );
                save([file_sdir file '_downTrackdata.mat'],'downTrackdata');
                cdownTrackdata{n} = downTrackdata;
            else
                M = load([file_sdir file '_downTrackdata.mat']);
                cdownTrackdata{n} = M.downTrackdata;
            end
            SdownTrackdata.phi = [SdownTrackdata.phi; cdownTrackdata{n}.phi];
            SdownTrackdata.x = [SdownTrackdata.x; cdownTrackdata{n}.x];
            SdownTrackdata.y = [SdownTrackdata.y; cdownTrackdata{n}.y];
            SdownTrackdata.speed = [SdownTrackdata.speed; cdownTrackdata{n}.speed];
            SdownTrackdata.r = [SdownTrackdata.r; cdownTrackdata{n}.r];
            SdownTrackdata.time = [SdownTrackdata.time; cdownTrackdata{n}.time];
        end
        clear downTrackdata
        downTrackdata = SdownTrackdata;
        clear SdownTrackdata

        % plot consolidated tracking data
        if ~exist([grp_trackdir mouseid '_' expname '_mtrajectory.fig'],'file')
            fig = figure;
            for n = 1:Nfiles
                plot(cdownTrackdata{n}.x, cdownTrackdata{n}.y, 'b'); hold on;
            end
            hold off; axis off;
            title('Mouse trajectory','Fontweight','normal','Fontsize',12);
            savefig(fig,[grp_trackdir mouseid '_' expname '_mtrajectory']);
            saveas(fig,[grp_trackdir mouseid '_' expname '_mtrajectory'],'png');
            close(fig);
        end
        
        % save consolidated tracking data
        fprintf('%s: Saving tracking data\n', [mouseid '_' expname]);
        grp_sname = [grp_trackdir mouseid '_' expname '_downTrackdata.mat'];
        save(grp_sname,'-struct','downTrackdata');
        clear cdownTrackdata 
    else
        grp_sname = [grp_trackdir mouseid '_' expname '_downTrackdata.mat'];
        downTrackdata = load(grp_sname);
        fprintf('%s: Tracking data found and loaded\n', [mouseid '_' expname]);
    end
else
    fprintf('%s: Behaviour tracking step not ticked. Skipping this and later steps.\n', [mouseid '_' expname]);
    t = toc;
    str = sprintf('%s: Processing done in %g hrs\n', [mouseid '_' expname], round(t/3600,2));
    cprintf(str)
    return
end
    
%% 6) PFmapping

if dostep(6) && ~contains(list,'-')                          
    fprintf('%s: Run PFmapping on individual runs, not on whole experiment. Skipping step.\n', [mouseid '_' expname]); 
    t = toc;
    str = sprintf('%s: Processing done in %g hrs\n', [mouseid '_' expname], round(t/3600,2));
    cprintf(str)
    return
end

if dostep(6)
    if any(downTrackdata.r < 100)
        params.mode_dim = '2D';                     % open field
        params.PFmap.Nbins = params.PFmap.Nbins_2D; % number of location bins in [x y]               
    else 
        params.mode_dim = '1D';                     % circular linear track
        params.PFmap.Nbins = params.PFmap.Nbins_1D; % number of location bins         
        
        if ~exist([grp_trackdir mouseid '_' expname '_phiposition.fig'],'file')
            fig = figure('position',[680 678 1000 200]);
            plot(downTrackdata.phi); hold on;
            plot(downTrackdata.time); hold off;
            legend('phi','time'); legend('boxoff');
            title('Mouse phi position','Fontweight','normal','Fontsize',12);
            savefig(fig,[grp_trackdir mouseid '_' expname '_phiposition']);
            saveas(fig,[grp_trackdir mouseid '_' expname '_phiposition'],'png');
            close(fig);
        end
    end
    fields = {'Nbins_1D','Nbins_2D'};
    params.PFmap = rmfield(params.PFmap,fields);

    [ hist, asd, PFdata, hist_epochs, asd_epochs, PFdata_epochs, params ] = ...
        neuroSEE_mapPF( spikes, downTrackdata, data_locn, [], params, force(6), list, reffile, conc_runs);

    %% Saving all data
    sname_allData = [ grp_sdir mouseid '_' expname '_ref' reffile '_' mcorr_method '_' segment_method '_' str_fissa...
                      '_allData_blprctile' num2str(bl_prctile) '.mat' ];

    fprintf('%s: Saving all data\n', [mouseid '_' expname]);
    save(sname_allData, 'list','reffile','template_g','template_r',...
                        'corr_image','masks','tsG','df_f',...
                        'spikes','bl','neuron_sn','g',...
                        'downTrackdata','PFdata','hist','asd','params');
    if ~isempty(dtsG), save(sname_allData,'-append','dtsG'); end
    if ~isempty(ddf_f), save(sname_allData,'-append','ddf_f'); end
else
    fprintf('%s: PF mapping step not ticked. Skipping step.\n', [mouseid '_' expname]);
end

t = toc;
str = sprintf('%s: Processing done in %g hrs. No errors!\n', [mouseid '_' expname], round(t/3600,2));
cprintf(str)

% Send slack notifcation if processing has finished
if slacknotify
    slacktext = [mouseid '_' expname ': FINISHED in' num2str(round(t/3600,2)) ' hrs. No errors!'];
    SendSlackNotification( slacktext );
end

end
