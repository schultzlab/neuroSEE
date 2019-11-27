function frun_mcorr_batch( array_id, list, force )

% Written by Ann Go
% 
% This script runs the motion correction algorithm for a list of image 
% files. 

% clear; % close all;
tic

if nargin<3, force = false; end

%% USER: Set basic settings
                            
default = true;                 % flag to use default parameters
mcorr_method = 'normcorre';     % [normcorre,fftRigid] CaImAn NoRMCorre method, fft-rigid method (Katie's)

%% Load module folders and define data directory

[data_locn,comp,err] = load_neuroSEEmodules;
if ~isempty(err)
    beep
    cprintf('Errors',err);    
    return
end
if strcmpi(comp,'hpc')
    maxNumCompThreads(32);      % max # of computational threads, must be the same as # of ncpus specified in jobscript (.pbs file)
end

% Send Ann slack message 
if array_id == 1 
    slacktext = [list(1:end-5) ': processing 1st file'];
    neuroSEE_slackNotify( slacktext );
end

%% Image file

listfile = [data_locn 'Digital_Logbook/lists/' list];
files = extractFilenamesFromTxtfile(listfile);

file = files(array_id,:); 


%% USER: Set parameters (if not using default)

if ~default
    % motion correction
        % neurosee method
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
end


%% Default and non-user defined parameters

if default
    params = load( 'default_params.mat', 'fftRigid', 'nonrigid', 'rigid' );
    
    % Remove irrelevant parameters 
    if strcmpi(mcorr_method,'normcorre')
        params = rmfield(params,'fftRigid');
    elseif strcmpi(mcorr_method,'normcorre-r')
        fields = {'fftRigid','nonrigid'};
        params = rmfield(params,fields);
    elseif strcmpi(mcorr_method,'normcorre-nr')
        fields = {'fftRigid','rigid'};
        params = rmfield(params,fields);
    elseif strcmpi(mcorr_method,'fftRigid')
        fields = {'rigid','nonrigid'};
        params = rmfield(params,fields);
    else
        beep
        cprintf('Errors','Invalid motion correction method.');    
        return
    end
end


%% Check if file has been processed. 

dir_proc = [data_locn 'Data/' file(1:8) '/Processed/' file '/'];
        
check = 0;
if exist(dir_proc,'dir')
    dir_mcorr = [dir_proc 'mcorr_' mcorr_method '/'];
    if all( [exist([dir_mcorr file '_2P_XYT_green_mcorr.tif'],'file'),...
             exist([dir_mcorr file '_2P_XYT_red_mcorr.tif'],'file'),...
             exist([dir_mcorr file '_mcorr_output.mat'],'file')] )
        check = 1;
    end
end


%% Motion correction
% Saved in file folder: motion corrected tif files
%                       summary fig & png, 
%                       mat with fields 
%                           green.[ meanframe, meanregframe ]
%                           red.[ meanframe, meanregframe ]
%                           template
%                           shifts
%                           col_shift
%                           params

if force || ~check(1)
    [imG,imR] = load_imagefile( data_locn, file );
    params.methods.mcorr_method = mcorr_method;
    if or(strcmpi(mcorr_method,'normcorre'),strcmpi(mcorr_method,'normcorre-r'))
        params.rigid.d1 = size(imG,1);
        params.rigid.d2 = size(imG,2);
        params.rigid.grid_size = [params.rigid.d1, params.rigid.d2, 1];
    end
    if or(strcmpi(mcorr_method,'normcorre'),strcmpi(mcorr_method,'normcorre-nr'))
        params.nonrigid.d1 = size(imG,1);
        params.nonrigid.d2 = size(imG,2);
    end
    [~, ~, ~, ~] = neuroSEE_motionCorrect( imG, imR, data_locn, file, params, force );
else 
    fprintf('%s: Motion corrected files found. Skipping motion correction\n', file);
end

%% Send Ann slack message 
if array_id == size(files,1)
    slacktext = [list(1:end-5) ': processing LAST file'];
    neuroSEE_slackNotify( slacktext );
end

t = toc;
str = sprintf('%s: Processing done in %g hrs\n', file, round(t/3600,2));
cprintf(str)

end
