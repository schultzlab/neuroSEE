% Written by Ann Go
% 
% This script runs the complete data processing pipeline for a batch of
% image files. Processing steps for EACH include:
% (1) motion correction (and dezippering)
% (2) roi segmentation 
% (3) neuropil decontamination and timeseries extraction
% (4) spike extraction
% (5) tracking data extraction
% (6) place field mapping
% Steps 3-6 are not necessary at the file level.
%
% This script is designed to be run on the hpc server where array_id loops
% through values specified in the hpc script.
%
% INPUTS
% array_id  : number which serves as array index for files in 'list'
% list      : name of text file containing filenames of files to be processed.
%           May be in the format 'list_m##_expname.txt'. However, when
%           processing tens of lists (hundreds of files), it is easier to
%           create one superset list (e.g. 'list_newfiles.txt') to process.
% dostep    : (optional, default: [1;1;0;0;0;0]) Set to 0 or 1 to skip (0) 
%               or implement (1) a processing step (see list at the top)
% force     : (optional, default: [0;0;0;0;0;0]) Set to 1 to force the implementation 
%               of a processing step (even if output already exists)
% min_SNR   : (optional, default: 2.5) CaImAn parameter. Minimum SNR for accepting exceptional events. 
%               Typically 3 for 330x330um FOV data, 2.5 for 490x490um FOV data.
%
% Matlab version requirement: 
%   On the hpc, normcorre only works with Matlab R2017a. 
%   CaImAn requires at least Matlab R2017b
%   FISSA requires at least Matlab R2018

function frun_pipeline_file_batch( array_id, list, dostep, force, min_SNR )
if nargin<5, min_SNR = 2.5; end
if nargin<4, force  = [0;0;0;0;0;0]; end
if nargin<3, dostep = [1;1;0;0;0;0]; end

tic

% Load module folders and define data directory
[data_locn,comp,err] = load_neuroSEEmodules;
if ~isempty(err)
    beep
    cprintf('Errors',err);    
    return
end

% Image file
listfile = [data_locn 'Digital_Logbook/lists_imaging/' list];
files = extractFilenamesFromTxtfile( listfile );
file = files(array_id,:);

% Auto-defined 
if str2double(file(1,1:4)) > 2018
    FOV = 490;                    % FOV area = FOV x FOV, FOV in um
else
    FOV = 330; 
end
% virus (which determines decay time constant for calcium transient)
imdate = datevec(str2double(file(1:8)));
if datetime(imdate) >= datetime(datevec(str2double('20201208')))
    virus = 'jGCaMP7s';
else
    virus = 'GCaMP6s';
end


%% SETTINGS
% Basic settings
mcorr_method = 'normcorre';     % values: [normcorre, normcorre-r, normcorre-nr, fftRigid] 
                                    % CaImAn NoRMCorre method: 
                                    %   normcorre (rigid + nonrigid) 
                                    %   normcorre-r (rigid),
                                    %   normcorre-nr (nonrigid), 
                                    % fft-rigid method (Katie's)
segment_method = 'CaImAn';      % [ABLE,CaImAn]    
runpatches = false;             % for CaImAn processing, flag to run patches (default: false)
dofissa = true;
manually_refine_spikes = false; % flag to manually refine spikes
doasd = false;                  % flag to do asd pf calculation

% Processing parameters (any parameter that is not set gets a default value)
% Add any parameters you want to set after FOV. See neuroSEE_setparams for
% full set of parameters
params = neuroSEE_setparams(...
            'mcorr_method', mcorr_method,...
            'segment_method', segment_method,...
            'runpatches', runpatches,...
            'dofissa', dofissa,...
            'manually_refine_spikes', manually_refine_spikes,...
            'doasd', doasd,...
            'FOV', FOV,...
            'virus', virus,...
            'min_SNR', min_SNR); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Matlab version
release = version('-release'); % Find out what Matlab release version is running
MatlabVer = str2double(release(1:4));


%% Check if file has been processed. If not, continue processing unless forced to overwrite 
% existing processed data
% check(1:6) check for existing data in processing steps 1-6
% check(7) checks for existing mat file pooling all processed data for file

check = checkforExistingProcData(data_locn, file, params);

% Some security measures
force = logicalForce(force);        % Only allow combinations of force/step values that make sense
dostep = logicaldostep(dostep);     % because later steps require earlier ones

if ~any(force) && check(7)
    fprintf('%s: File already processed\n', file)
    return
end


%% Image files
% Load original image files if forced to do motion correction or if motion corrected files don't exist 

if dostep(1)
    if force(1) || ~check(1)
        % Continue only if Matlab version on hpc is R2017
        if strcmpi(comp,'hpc') && MatlabVer > 2017
            beep
            err = sprintf('%s: Lower Matlab version required for motion correction. Cannot proceed.\n', file);
            cprintf('Errors',err);
            return
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
        if strcmpi(segment_method,'CaImAn') 
            [ imG, ~, params.mcorr ] = neuroSEE_motionCorrect( imG, imR, data_locn, file, mcorr_method, params.mcorr, reffile, force(1) );
            imR = []; % CaImAn does not use imR
        else
            [ imG, ~, params.mcorr, imR ] = neuroSEE_motionCorrect( imG, imR, data_locn, file, mcorr_method, params.mcorr, reffile, force(1) );
        end
    else 
        fprintf('%s: Motion corrected files found. Skipping motion correction\n', file);
        imG = []; imR = [];
    end
else
    fprintf('%s: No processing steps ticked. Cannot proceed.\n', file);
    return
end


%% (2) ROI segmentation
% Saved in file folder: correlation image with ROIs (fig, png) 
%                       summary plots of tsG, df_f (fig, png)
%                       mat with fields {tsG, df_f, masks, corr_image, params}

if dostep(2)
    [tsG, df_f, masks, corr_image, params] = neuroSEE_segment( imG, data_locn, params, file, [], false, force(2), mean(imR,3) );
    clear imG imR
else
    fprintf('%s: ROI segmentation step not ticked. Skipping this and later steps.\n', file);
    t = toc;
    str = sprintf('%s: Processing done in %g hrs\n', file, round(t/3600,2));
    cprintf(str)
    return
end


%% (3) Run FISSA to extract neuropil-corrected time-series
% Saved in file folder: mat file with fields {dtsG, ddf_f, masks}
%                       summary plots of tsG, ddf_f (fig & png)

if dostep(3)
    if dofissa
        release = version('-release'); % Find out what Matlab release version is running
        MatlabVer = str2double(release(1:4));
        if MatlabVer > 2017
            [dtsG, ddf_f, params] = neuroSEE_neuropilDecon( masks, data_locn, file, params, force(3) );
        else
            fprintf('%s: Higher Matlab version required, skipping FISSA.', file);
            dofissa = false; params.methods.dofissa = false;
            dtsG = [];
            ddf_f = [];
        end
    else
        dtsG = [];
        ddf_f = [];
    end
else
    if dofissa
        fprintf('%s: FISSA step not ticked. Skipping this and later steps.\n', file);
        t = toc;
        str = sprintf('%s: Processing done in %g hrs\n', file, round(t/3600,2));
        cprintf(str)
        return
    end
end

%% (4) Estimate spikes
% Saved in file folder: mat with fields {tsG, dtsG, df_f, ddf_f, spikes, params}
%                       summary plot of spikes (fig & png)

if dostep(4)
    [spikes, params] = neuroSEE_extractSpikes( df_f, ddf_f, data_locn, file, params, force(4) );

    if manually_refine_spikes
        GUI_manually_refine_spikes( spikes, tsG, dtsG, df_f, ddf_f, data_locn, file, params, corr_image, masks );
        uiwait 
    end
else
    fprintf('%s: Spike estimation step not ticked. Skipping this and later steps.\n', file);
    t = toc;
    str = sprintf('%s: Processing done in %g hrs\n', file, round(t/3600,2));
    cprintf(str)
    return
end


%% (5) Find tracking file then load it
% Saved in file folder: trajectory plot (fig & png)
%                       mat with fields {time, r, phi, x, y , speed, w, alpha, TTLout, filename}

if dostep(5)
    trackfile = findMatchingTrackingFile( data_locn, file, force(5) );
    if force(5) || ~check(5)
        Trackdata = load_trackfile(data_locn, file, trackfile, force(5));
        downTrackdata = downsample_trackData( Trackdata, size(spikes,2), params.fr );
        save([data_locn 'Data/' file(1:8) '/Processed/' file '/behaviour/' file '_downTrackdata.mat'],'downTrackdata');
    else
        M = load([data_locn 'Data/' file(1:8) '/Processed/' file '/behaviour/' file '_downTrackdata.mat']);
        downTrackdata = M.downTrackdata;
    end
else
    fprintf('%s: Behaviour tracking step not ticked. Skipping this and later steps.\n', file);
    t = toc;
    str = sprintf('%s: Processing done in %g hrs\n', file, round(t/3600,2));
    cprintf(str)
    return
end

%% (6) Generate place field maps
% Saved: fig & png of several figures showing occMap, spkMap, pfMap, remapping and spkMap_pertrial
%        mat file with fields {occMap, hist, asd, downData, activeData}

if dostep(6)
    if manually_refine_spikes
        force(6) = true;
    end

    if any(downTrackdata.r < 100)
        params.mode_dim = '2D';                     % open field
        params.PFmap.Nbins = params.PFmap.Nbins_2D; % number of location bins in [x y]               
    else 
        params.mode_dim = '1D';                     % circular linear track
        params.PFmap.Nbins = params.PFmap.Nbins_1D; % number of location bins   
            
        if ~exist([data_locn 'Data/' file(1:8) '/Processed/' file '/behaviour/' file '_phiposition.fig'],'file')
            fig = figure('position',[680 678 1000 200]);
            plot(downTrackdata.phi); 
            title('Mouse phi position','Fontweight','normal','Fontsize',12);
            savefig(fig, [data_locn 'Data/' file(1:8) '/Processed/' file '/behaviour/' file '_phiposition']);
            saveas(fig, [data_locn 'Data/' file(1:8) '/Processed/' file '/behaviour/' file '_phiposition'],'png');
            close(fig);
        end
    end
    fields = {'Nbins_1D','Nbins_2D'};
    params.PFmap = rmfield(params.PFmap,fields);

    if strcmpi(params.mode_dim,'1D')
        [ hist, asd, PFdata, params ] = neuroSEE_mapPF( spikes, downTrackdata, data_locn, file, params, force(6));
    else
        [ occMap, hist, asd, downData, params, spkMap, spkIdx ] = neuroSEE_mapPF( spikes, trackData, data_locn, file, params, force(6));
    end

    %% Save output. These are all the variables needed for viewing data with GUI

    if dofissa
        str_fissa = 'FISSA';
    else
        str_fissa = 'noFISSA';
    end
    fname_allData = [ data_locn 'Data/' file(1:8) '/Processed/' file '/' file '_' mcorr_method '_' segment_method '_' str_fissa '_allData.mat'];

    save(fname_allData,'file','corr_image','masks','tsG','df_f','spikes','trackfile',...
                        'downTrackdata','pfData','hist','asd','params');
    if ~isempty(dtsG), save(fname_allData,'-append','dtsG'); end
    if ~isempty(ddf_f), save(fname_allData,'-append','ddf_f'); end
else
    fprintf('%s: PF mapping step not ticked. Skipping step.\n', file);
end

t = toc;
str = sprintf('%s: Processing done in %g hrs\n', file, round(t/3600,2));
cprintf(str)

end

