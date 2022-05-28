function frun_redoFISSA(list, reffile)

[data_locn,~,err] = load_neuroSEEmodules;
if ~isempty(err)
    beep
    cprintf('Errors',err);    
    return
end

%% Compare number of ROIs in segment_output and FISSA_output
[mouseid,expname,fov] = find_mouseIDexpname( list );
dir_concenv = [data_locn 'Analysis/' mouseid '/' fov '/' mouseid '_' expname ...
                '/group_proc/imreg_normcorre_CaImAn/' mouseid '_' expname '_imreg_ref' reffile '/'];
if ~exist(dir_concenv,'dir')            
    cprintf('Errors','Experiment has not been processed.\n');
    return
end

so = load([dir_concenv mouseid '_' expname '_ref' reffile '_segment_output.mat']);
if exist([dir_concenv '/FISSA/' mouseid '_' expname '_ref' reffile '_fissa_output.mat'],'file')
    fissa = load([dir_concenv '/FISSA/' mouseid '_' expname '_ref' reffile '_fissa_output.mat']);
    num_rois_so = size(so.masks,3);
    num_rois_fissa = size(fissa.dtsG,1);

    if num_rois_so ~= num_rois_fissa 
        fprintf('%s: ROI numbers not the same for segment_output and fissa_output. Redoing FISSA.\n', [mouseid '_' expname '_ref' reffile]);
        frun_pipeline_imreg( list, reffile, true, [0; 0; 1; 0; 0; 0], [1; 1; 1; 1; 0; 0] )
    else
        fprintf('%s: ROI numbers are the same for segment_output and fissa_output. No need to redo FISSA.\n', [mouseid '_' expname '_ref' reffile]);
    end
else
    fprintf('%s: FISSA_output does not exist. Now doing FISSA.\n', [mouseid '_' expname '_ref' reffile]);
    frun_pipeline_imreg( list, reffile, true, [0; 0; 1; 0; 0; 0], [1; 1; 1; 1; 0; 0] )
end
end