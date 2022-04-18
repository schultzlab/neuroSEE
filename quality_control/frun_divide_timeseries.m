function frun_divide_timeseries( list, reffile, numfiles, force )

if nargin<4, force = false; end

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
else
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

%% Copy files from dir_concenv to folder for sub experiments
dir_concenv = [data_locn 'Analysis/' mouseid '/' fov '/' mouseid '_' expname ...
                '/group_proc/imreg_normcorre_CaImAn/' mouseid '_' expname '_imreg_ref' reffile '/'];
if ~exist(dir_concenv,'dir')            
    cprintf('Errors','Experiment has not been processed.\n');
    return
end

fpf = load([dir_concenv mouseid '_' expname '_ref' reffile '_framesperfile.mat']);
so = load([dir_concenv mouseid '_' expname '_ref' reffile '_segment_output.mat']);
fissa = load([dir_concenv '/FISSA/' mouseid '_' expname '_ref' reffile '_fissa_output.mat']);

ff = 1; a = 1; b = 0;
for j = 1:numel(env)
    dir_env = [data_locn 'Analysis/' mouseid '/' fov '/' mouseid '_' expname '-' env{j}...
                '/group_proc/imreg_normcorre_CaImAn/' mouseid '_' expname '-' env{j} ...
                '_imreg_ref' reffile '_concenvrois/'];
    if ~exist(dir_env,'dir') 
        mkdir( dir_env ); fileattrib(dir_env,'+w','g','s'); 
        mkdir( [dir_env 'FISSA/'] ); fileattrib(dir_env,'+w','g','s'); 
    end
    if force || ~all([exist([dir_env mouseid '_' expname '-' env{j} '_ref' reffile '_framesperfile.mat'],'file'),...
                      exist([dir_env mouseid '_' expname '-' env{j} '_ref' reffile '_segment_output.mat'],'file'),...
                      exist([dir_env 'FISSA/' mouseid '_' expname '-' env{j} '_ref' reffile '_fissa_output.mat'],'file')])
            
        fprintf('Copying files to %s.\n', [mouseid '_' expname '-' env{j} '_ref' reffile]);
        copyfile([dir_concenv mouseid '_' expname '_ref' reffile '_df_f.fig'],...
                 [dir_env mouseid '_' expname '-' env{j} '_ref' reffile '_df_f.fig'])
        copyfile([dir_concenv mouseid '_' expname '_ref' reffile '_df_f.png'],...
                 [dir_env mouseid '_' expname '-' env{j} '_ref' reffile '_df_f.png'])
        copyfile([dir_concenv mouseid '_' expname '_ref' reffile '_elimROIs.fig'],...
                 [dir_env mouseid '_' expname '-' env{j} '_ref' reffile '_elimROIs.fig'])
        copyfile([dir_concenv mouseid '_' expname '_ref' reffile '_elimROIs.png'],...
                 [dir_env mouseid '_' expname '-' env{j} '_ref' reffile '_elimROIs.png'])
        copyfile([dir_concenv mouseid '_' expname '_ref' reffile '_imreg_template.fig'],...
                 [dir_env mouseid '_' expname '-' env{j} '_ref' reffile '_imreg_template.fig'])
        copyfile([dir_concenv mouseid '_' expname '_ref' reffile '_imreg_template.mat'],...
                 [dir_env mouseid '_' expname '-' env{j} '_ref' reffile '_imreg_template.mat'])
        copyfile([dir_concenv mouseid '_' expname '_ref' reffile '_imreg_template.png'],...
                 [dir_env mouseid '_' expname '-' env{j} '_ref' reffile '_imreg_template.png'])
        copyfile([dir_concenv mouseid '_' expname '_ref' reffile '_raw_timeseries.fig'],...
                 [dir_env mouseid '_' expname '-' env{j} '_ref' reffile '_raw_timeseries.fig'])
        copyfile([dir_concenv mouseid '_' expname '_ref' reffile '_raw_timeseries.png'],...
                 [dir_env mouseid '_' expname '-' env{j} '_ref' reffile '_raw_timeseries.png'])
        copyfile([dir_concenv mouseid '_' expname '_ref' reffile '_ROIs.fig'],...
                 [dir_env mouseid '_' expname '-' env{j} '_ref' reffile '_ROIs.fig'])
        copyfile([dir_concenv mouseid '_' expname '_ref' reffile '_ROIs.png'],...
                 [dir_env mouseid '_' expname '-' env{j} '_ref' reffile '_ROIs.png'])
        copyfile([dir_concenv '/FISSA/' mouseid '_' expname '_ref' reffile '_fissa_df_f.fig'],...
                 [dir_env '/FISSA/' mouseid '_' expname '-' env{j} '_ref' reffile '_fissa_df_f.fig'])
        copyfile([dir_concenv '/FISSA/' mouseid '_' expname '_ref' reffile '_fissa_df_f.png'],...
                 [dir_env '/FISSA/' mouseid '_' expname '-' env{j} '_ref' reffile '_fissa_df_f.png'])
        copyfile([dir_concenv '/FISSA/' mouseid '_' expname '_ref' reffile '_fissa_timeseries.fig'],...
                 [dir_env '/FISSA/' mouseid '_' expname '-' env{j} '_ref' reffile '_fissa_timeseries.fig'])
        copyfile([dir_concenv '/FISSA/' mouseid '_' expname '_ref' reffile '_fissa_timeseries.png'],...
                 [dir_env '/FISSA/' mouseid '_' expname '-' env{j} '_ref' reffile '_fissa_timeseries.png'])

        % framesperfile
        framesperfile = fpf.framesperfile;
        framesperfile = framesperfile(ff:ff+numfiles(j)-1);
        ff = ff + numfiles(j);
        save([dir_env mouseid '_' expname '-' env{j} '_ref' reffile '_framesperfile.mat'],'framesperfile');

        % segment output
        so2 = so;
        b = b + sum(framesperfile);
        so2.tsG = so.tsG(a:b);
        so2.elim_tsG = so.elim_tsG(a:b);
        so2.df_f = so.df_f(a:b);
        so2.elim_df_f = so.elim_df_f(a:b);
        fprintf('Dividing tsG and df and saving.\n');
        save([dir_env mouseid '_' expname '-' env{j} '_ref' reffile '_segment_output.mat'],'-struct','so2');

        % fissa output
        fissa2 = fissa;
        fissa2.dtsG = fissa.dtsG(a:b);
        fissa2.ddf_f = fissa.ddf_f(a:b);
        a = b+1;
        fprintf('Dividing fissa timeseries and saving.\n');
        save([dir_env 'FISSA/' mouseid '_' expname '-' env{j} '_ref' reffile '_fissa_output.mat'],'-struct','fissa2');
    else
        fprintf('Timeseries distribution has been processed.\n')
    end
end

