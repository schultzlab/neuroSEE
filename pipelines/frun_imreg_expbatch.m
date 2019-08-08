% Written by Ann Go

function frun_imreg_expbatch( array_id, expname, list, reffile, slacknotify, force )

if nargin<6, force = false; end
if nargin<5, slacknotify = false; end
tic

%% Load module folders and define data directory

test = false;                   % flag to use one of smaller files in test folder)
mcorr_method = 'normcorre';     % [fftRigid, normcorre]
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

files = extractFilenamesFromTxtfile( list );

% Send Ann slack message
if slacknotify
    if array_id == 1
        slacktext = [expname ': registering 1 of ' num2str(size(files,1)) 'files'];
        SendSlackNotification('https://hooks.slack.com/services/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', ...
           slacktext, '@m.go', ...
           [], [], [], []);   
    end
end

% image to be registered
file = files(array_id,:);
params.methods.mcorr_method = mcorr_method;

if ~strcmpi(file,reffile)

    filedir = [data_locn 'Data/' file(1:8) '/Processed/' file '/mcorr_normcorre_ref' reffile '/'];
    if force || ~exist([filedir file '_imreg_ref' reffile '_output.mat'],'file')
        [ imG, ~ ] = load_imagefile( data_locn, file );
    else
        imG = [];
    end

    % template file
    filedir = [data_locn 'Data/' reffile(1:8) '/Processed/' reffile '/mcorr_' mcorr_method '/'];
    c = load([filedir reffile '_mcorr_output.mat']);
    template = c.template;
    clear c


    %% Image registration 
    params.nonrigid = NoRMCorreSetParms(...
                'd1',size(template,1),...
                'd2',size(template,2),...
                'grid_size',[32,32],...
                'overlap_pre',[32,32],...
                'overlap_post',[32,32],...
                'iter',1,...
                'use_parallel',false,...
                'max_shift',50,...
                'mot_uf',4,...
                'bin_width',200,...
                'max_dev',3,...
                'us_fac',50,...
                'init_batch',200);

    try
        [ ~, ~, ~ ] = neuroSEE_imreg( imG, data_locn, file, reffile, template, params, force );
    catch
        try
            fprintf('%s: Error in image registration, changing max_shift to 40\n', file);
            params.nonrigid.max_shift = 40;
            [ ~, ~, ~ ] = neuroSEE_imreg( imG, data_locn, file, reffile, template, params, force );
        catch
            try
                fprintf('%s: Error in image registration, changing max_shift to 30\n', file);
                params.nonrigid.max_shift = 30;
                [ ~, ~, ~ ] = neuroSEE_imreg( imG, data_locn, file, reffile, template, params, force );
            catch
                fprintf('%s: Error in image registration, changing max_shift to 25\n', file);
                params.nonrigid.max_shift = 25;
                [ ~, ~, ~ ] = neuroSEE_imreg( imG, data_locn, file, reffile, template, params, force );
            end
        end
    end

    % Send Ann slack message
    if slacknotify
        if array_id == size(files,1)
            slacktext = [expname ': registering ' num2str(size(files,1)) ' of ' num2str(size(files,1)) 'files'];
            SendSlackNotification('https://hooks.slack.com/services/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', ...
               slacktext, '@m.go', ...
               [], [], [], []);   
        end
    end
    
    t = toc;
    str = sprintf('%s: Processing done in %g hrs\n', file, round(t/3600,2));
    cprintf(str)

end

end
