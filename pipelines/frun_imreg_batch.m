% Written by Ann Go

function frun_imreg_batch( array_id, list, reffile, refChannel, slacknotify, force )

if nargin<6, force = false; end
if nargin<5, slacknotify = false; end
if nargin<4, refChannel = 'green'; end
% if nargin<3, see line 29
tic

%% Load module folders and define data directory

test = false;                       % flag to use one of smaller files in test folder)
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

%% Files

listfile = [data_locn 'Digital_Logbook/lists/' list];
files = extractFilenamesFromTxtfile( listfile );
if nargin<3, reffile = files(1,:); end

% Send Ann slack message
if slacknotify
    if array_id == 1
        slacktext = [list(6:end-4) ': registering 1 of ' num2str(size(files,1)) 'files'];
        neuroSEE_slackNotify( slacktext );
    end
end

% image to be registered
file = files(array_id,:);
params.methods.mcorr_method = mcorr_method;
params.refChannel = refChannel;

if ~strcmpi( file, reffile )

    filedir = [data_locn 'Data/' file(1:8) '/Processed/' file '/mcorr_' mcorr_method '_ref' reffile '/'];
    if force || ~exist([filedir file '_imreg_ref' reffile '_output.mat'],'file')
        imG = read_file([ data_locn 'Data/' file(1:8) '/2P/' file '_2P/' file '_2P_XYT_green.tif']);
        imR = read_file([ data_locn 'Data/' file(1:8) '/2P/' file '_2P/' file '_2P_XYT_red.tif']);
    else
        imG = zeros(512,512);
        imR = zeros(512,512);
    end


    %% Image registration 
    params.nonrigid = NoRMCorreSetParms(...
                'd1',size(imG,1),...
                'd2',size(imG,2),...
                'grid_size',[64,64],...
                'overlap_pre',[32,32],...
                'overlap_post',[32,32],...
                'iter',1,...
                'use_parallel',false,...
                'max_shift',20,...
                'mot_uf',4,...
                'bin_width',200,...
                'max_dev',3,...
                'us_fac',50,...
                'init_batch',200);

    try
        [ ~, ~, ~, ~ ] = neuroSEE_imreg( imG, imR, data_locn, file, reffile, params, force );
    catch
        fprintf('%s: Error in image registration, changing max_shift to 15\n', file);
        params.nonrigid.max_shift = 15;
        [ ~, ~, ~, ~ ] = neuroSEE_imreg( imG, imR, data_locn, file, reffile, params, force );
    end

    % Send Ann slack message
    if slacknotify
        if array_id == size(files,1)
            slacktext = [list(6:end-4) ': registering ' num2str(size(files,1)) ' of ' num2str(size(files,1)) 'files'];
            neuroSEE_slackNotify( slacktext );
        end
    end
    
    t = toc;
    str = sprintf('%s: Processing done in %g hrs\n', file, round(t/3600,2));
    cprintf(str)

end

end