% Written by Ann Go
%
% This script implements motion correction of a single file or the
% registration of a file to a reference file (if reference file is specified). 
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
% force     : (optional) flag to force generation of comparison figures even though they
%               already exist. Default: false
% reffile   : (optional) file to be used as registration template. When empty, motion 
%               correction of file is done, with a template computed from first 200 frames.
%               This file is usually part of 'list' (i.e. for best results, choose a 
%               reference file from the same experiment) but does not have to be. This 
%               file must have already been motion corrected.
% refChannel : (optional) channel (red or green) to be used as registation
%               template. Default: 'green'
% slacknotify : (optional) flag to send Slack notification when processing is started
%               or has ended. Default: false


function frun_imreg_batch( array_id, list, force, reffile, refChannel, slacknotify )

if nargin<6, slacknotify = false; end
if nargin<5, refChannel = 'green'; end
if nargin<4, reffile = []; end
if nargin<3, force = false; end
tic

%% Load module folders and define data directory

test = false;                      % flag to use one of smaller files in test folder)
default = true;                    % flag to use default motion correction parameters
mcorr_method = 'normcorre-nr';     % [fftRigid, normcorre, normcorre-nr, normcorre-r]

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


%% USER: Set parameters (if not using default)
if ~default
    params_mcorr.refChannel = 'green';           % reference channel for motion correction [default: 'green']
    % NoRMCorre-rigid
    if or(strcmpi(mcorr_method,'normcorre'),strcmpi(mcorr_method,'normcorre-r'))
        params_mcorr.normcorre_r = NoRMCorreSetParms(...
            'd1', 512,...
            'd2', 512,...
            'max_shift',30,...          % default: 30
            'bin_width',200,...         % default: 200
            'us_fac',50,...             % default: 50
            'init_batch',200);          % default: 200
        params_mcorr.normcorre_r.print_msg = false;   % default: false
    end
    % NoRMCorre-nonrigid
    if or(strcmpi(mcorr_method,'normcorre'),strcmpi(mcorr_method,'normcorre-nr') )    
        params_mcorr.normcorre_nr = NoRMCorreSetParms(...
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
        params_mcorr.normcorre_nr.print_msg = false;    % default: false
    end
end
if default
    c = load('default_params.mat');
    if strcmpi(mcorr_method,'normcorre')
        params_mcorr.normcorre_r = c.mcorr.normcorre_r;
        params_mcorr.normcorre_nr = c.mcorr.normcorre_nr;
    elseif strcmpi(mcorr_method,'normcorre-r')
        params_mcorr.normcorre_r = c.mcorr.normcorre_r;
    elseif strcmpi(mcorr_method,'normcorre-nr')
        params_mcorr.normcorre_nr = c.mcorr.normcorre_nr; 
    elseif strcmpi(mcorr_method,'fftRigid')
        if ~isempty(reffile)
            str = sprintf('%s: Cannot do image registration with fftRigid method. Choose another method.', list);
            beep
            cprintf('Errors',str);    
            return;
        else
            params_mcorr.fftRigid = c.mcorr.fftRigid;
        end
    else
        beep
        cprintf('Errors','Invalid motion correction method.');    
        return
    end
end
params_mcorr.refChannel = refChannel;


%% Files
listfile = [data_locn 'Digital_Logbook/lists/' list];
files = extractFilenamesFromTxtfile( listfile );

% image to be registered
file = files(array_id,:);

% Send Ann slack message
if slacknotify
    if array_id == 1
        slacktext = [list(6:end-4) ': registering 1 of ' num2str(size(files,1)) 'files'];
        neuroSEE_slackNotify( slacktext );
    end
end


%% Check if file has been processed. 

check = false;
if isempty(reffile)
    filedir = [ data_locn 'Data/' file(1:8) '/Processed/' file '/mcorr_' mcorr_method '/' ];
        if ~exist( filedir, 'dir' ), mkdir( filedir ); end
    fname_tif_gr = [filedir file '_2P_XYT_green_mcorr.tif'];
    fname_tif_red = [filedir file '_2P_XYT_red_mcorr.tif'];
    fname_mat = [filedir file '_mcorr_output.mat'];
else
    filedir = [ data_locn 'Data/' file(1:8) '/Processed/' file '/mcorr_' mcorr_method '_ref' reffile '/' ];
        if ~exist( filedir, 'dir' ), mkdir( filedir ); end
    fname_tif_gr = [filedir file '_2P_XYT_green_imreg_ref' reffile '.tif'];
    fname_tif_red = [filedir file '_2P_XYT_red_imreg_ref' reffile '.tif'];
    fname_mat = [filedir file '_imreg_ref' reffile '_output.mat'];
end

if all( [exist(fname_tif_gr,'file'),...
         exist(fname_tif_red,'file'),...
         exist(fname_mat,'file')] )
    check = true;
end

if force || ~check    
    if ~isempty(reffile) && strcmpi(file,reffile)
        fprintf('%s: Same as reference file. Skipping image registration\n', file);    
        return
    end

    [imG,imR] = load_imagefile( data_locn, file );
    
    %% Motion correction/Image registration 
    neuroSEE_motionCorrect( imG, imR, data_locn, file, mcorr_method, params_mcorr, reffile, force );
    
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
