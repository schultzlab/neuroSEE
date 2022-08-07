% Written by Ann Go
%
% This script implements motion correction or image registration.
% When no reference file (reffile) is specified, motion correction 
% of all files within list is done. Otherwise, each file within the 
% list is registered to the reference file (reffile). 
% This script is designed to be run on the hpc server where array_id 
% loops through values specified in the hpc script.
%
% INPUTS
% array_id  : number which serves as array index for files in 'list'
% list      : name of text file containing filenames of files to be compared.
%      Typically in the format 'list_m##_expname.txt'.
%   e.g.    'list_m62_fam1nov-fam1.txt'         - all fam1 files in fam1nov experiment
%           'list_m62_fam1nov.txt'              - all files in fam1nov experiments
%           'list_m79_fam1_s1-5.txt'            - all fam1 files across 5 sessions           
%           'list_m86_open_s1-2.txt'            - all open field files across 2 sessions
% force     : (optional, default: false) flag to force motion correction or image 
%               registration even though processed data exist 
% reffile   : (optional) file to be used as registration template. When empty, motion 
%               correction of file is done, with a template computed from first 200 frames.
%               This file is usually part of 'list' (i.e. for best results, choose a 
%               reference file from the same experiment) but does not have to be. This 
%               file must have already been motion corrected.
% refChannel : (optional, default: 'green') channel (red or green) to be used as registation
%               template
% max_dev    : (optional, default: 5) Normcorre parameter. Maximum deviation of patch shift 
%               from rigid shift. Typical <14. Check for image distortion
%               when using high values.
% smallerpatch : (optional, default: false) Uses smaller patch sizes
%                       for hard-to-motion correct images.
% imregmode  : (optional, default: 2)
%                   1 - frame-by-frame registration with template, ...
%                   2 - registration of templates then application of shifts to other frames

function frun_mcorr_batch( array_id, list, force, reffile, refChannel, max_dev, smallerpatch, imregmode )

if nargin<8, imregmode = 2; end
if nargin<7, smallerpatch = false; end
if nargin<6, max_dev = 5; end
if nargin<5, refChannel = 'green'; end
if nargin<4, reffile = []; end
if nargin<3, force = false; end
tic

%% Load module folders and define data directory

test = false;                      % flag to use one of smaller files in test folder)
[data_locn,comp,err] = load_neuroSEEmodules(test);
if ~isempty(err)
    beep
    cprintf('Errors',err);    
    return
end
if strcmpi(comp,'hpc')
    maxNumCompThreads(32);      % max # of computational threads, must be the same as # of ncpus specified in jobscript (.pbs file)
end


%% Continue with motion correction or image registration if Matlab version is R2017
release = version('-release'); % Find out what Matlab release version is running
MatlabVer = str2double(release(1:4));
if strcmpi(comp,'hpc') && MatlabVer > 2017
    beep
    err = sprintf('%s: Lower Matlab version required; skipping processing.\n', list);
    cprintf('Errors',err);
    return
end

% image registration cannot be done with fftRigid method at the moment
if strcmpi(mcorr_method,'fftRigid')
    if ~isempty(reffile)
        str = sprintf('%s: Cannot do image registration with fftRigid method. Choose another method.\n', list);
        beep
        cprintf('Errors',str);    
        return;
    end
end


%% USER: Set parameters 
% Any parameter that is not set gets a default value

if smallerpatch
    grid_size_nr = [64,64];
    overlap_pre = [16,16]; 
    min_patch_size = [16,16];    
    min_diff = [8,8];
else
    grid_size_nr = [128,128];
    overlap_pre = [64,64]; 
    min_patch_size = [64,64];    
    min_diff = [32,32];

end

params = neuroSEE_setparams(...
            'mcorr_method', 'normcorre',...
            'refChannel', refChannel,...        % reference channel for motion correction
            'max_shift_r', 30,...       % maximum rigid shift
            'max_shift_nr', 30,...     % maximum nonrigid shift
            'max_dev', max_dev,...
            'grid_size_nr', grid_size_nr,...
            'overlap_pre', overlap_pre,... 
            'min_patch_size', min_patch_size,...      
            'min_diff', min_diff,...
            'imregmode', imregmode);       
        
params_mcorr = params.mcorr;

% For slack notifications
slacknotify = false;
% if true, set below
slackURL = '';
slackTarget = '@xxx';


%% Motion correction/Image registration 
listfile = [data_locn 'Digital_Logbook/lists_imaging/' list];
files = extractFilenamesFromTxtfile( listfile );

% Image to be registered
file = files(array_id,:);

% Send user slack message
if slacknotify
    if array_id == 1
        slacktext = [list(6:end-4) ': registering 1 of ' num2str(size(files,1)) 'files'];
        SendSlackNotification( slackURL, slacktext, slackTarget );
    end
end

% Check if file has been processed 
check = checkfor_mcorrIm( data_locn, file, mcorr_method, reffile );

if force || ~check    
    if ~isempty(reffile) % image registration
        if strcmpi(file,reffile)
            fprintf('%s: Same as reference file. Skipping image registration\n', file);    
            return
        end
    
        check_mcorr = checkfor_mcorrIm( data_locn, file, mcorr_method );
        if check_mcorr % motion-corrected image exists
            [imG,imR] = load_imagefile( data_locn, file, false, '_mcorr' );
            neuroSEE_motionCorrect2( imG, imR, data_locn, file, mcorr_method, params_mcorr, reffile, force, list, true, imregmode );
        else 
            [imG,imR] = load_imagefile( data_locn, file, false, [] );
            neuroSEE_motionCorrect2( imG, imR, data_locn, file, mcorr_method, params_mcorr, reffile, force, list, true, 1 );
        end
    else % motion correction
        [imG,imR] = load_imagefile( data_locn, file, false, [] );
        neuroSEE_motionCorrect2( imG, imR, data_locn, file, mcorr_method, params_mcorr, reffile, force, list, true, 1 );
    end
    
    if slacknotify
        if array_id == size(files,1)
            slacktext = [list(6:end-4) ': done registering ' num2str(size(files,1)) ' of ' num2str(size(files,1)) 'files'];
            neuroSEE_slackNotify( slacktext );
        end
    end

else 
    if isempty(reffile)
        fprintf('%s: Motion corrected files found. Skipping motion correction\n', file);
    else
        fprintf('%s: Registered images found. Skipping image registration to %s\n', file, reffile);
    end
end

t = toc;
str = sprintf('%s: Processing done in %g hrs\n', file, round(t/3600,2));
cprintf(str)
end
