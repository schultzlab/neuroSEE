function frun_redoFISSA(list, reffile)

[data_locn,~,err] = load_neuroSEEmodules;
if ~isempty(err)
    beep
    cprintf('Errors',err);    
    return
end

%% Find mouse ID and environment names
[mouseid,expname,fov] = find_mouseIDexpname( list );
env{1} = 'fam1';
ind1 = strfind(expname,'fam1');
if numel(ind1) > 2
    env{3} = 'fam1r2';
    if numel(strfind(expname,'fam1rev')) > 0
        env{2} = 'fam1rev';
    end
elseif numel(ind1) == 2
    if numel(strfind(expname,'fam1rev')) > 0
        env{2} = 'fam1rev';
    else
        env{3} = 'fam1r2';
    end
end
if numel(strfind(expname,'fam2')) > 0
    env{2} = 'fam2';
end
if numel(strfind(expname,'nov')) > 0
    env{2} = 'nov';
end

%% Compare number of ROIs in segment_output and FISSA_output
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
        frun_pipeline_imreg( list, reffile, true, [0; 0; 1; 0; 0; 0], [1; 1; 1; 1; 1; 1] )
    end
else
    frun_pipeline_imreg( list, reffile, true, [0; 0; 1; 0; 0; 0], [1; 1; 1; 1; 1; 1] )
end
end