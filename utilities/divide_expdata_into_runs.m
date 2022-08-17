% Written by Ann Go
% Script for copying files from experiment folder to folder for individual runs
% Segmentation and fissa outputs and extracted spike data for 
%  e.g. mXX_fam1fam2fam1/group_proc/imreg_normcorre_CaImAn/mXX_fam1fam2fam1
% are copied over to
%       mXX_fam1fam2fam1-fam1/group_proc/imreg_normcorre_CaImAn/mXX_fam1fam2fam1_concrunsrois
%       mXX_fam1fam2fam1-fam2/group_proc/imreg_normcorre_CaImAn/mXX_fam1fam2fam1_concrunsrois
%       mXX_fam1fam2fam1-fam1r2/group_proc/imreg_normcorre_CaImAn/mXX_fam1fam2fam1_concrunsrois

function divide_expdata_into_runs( data_locn, list, reffile, numfiles, dostep, force, dofissa )

if nargin<7, dofissa = false; end
    if dofissa
        str_fissa = 'FISSA';
    else
        str_fissa = 'noFISSA';
    end
if nargin<6, force = [0,0,0]; end
if nargin<5, dostep = [1,1,1]; end

% Mouseid, Experiment name, runs
[ mouseid, expname, fov ] = find_mouseIDexpname(list);

run{1} = 'fam1';
ind1 = strfind(expname,'fam1');
if numel(ind1) > 2
    run{3} = 'fam1r2';
    if numel(strfind(expname,'fam1rev')) > 0
        run{2} = 'fam1rev';
    end
elseif numel(ind1) == 2
    if numel(strfind(expname,'fam1rev')) > 0
        run{2} = 'fam1rev';
    else
        run{3} = 'fam1r2';
    end
end
if numel(strfind(expname,'fam2')) > 0
    run{2} = 'fam2';
end
if numel(strfind(expname,'nov')) > 0
    run{2} = 'nov';
end

% load data from experiment folder
if ~isempty(fov)
    dir_exp = [data_locn 'Analysis/' mouseid '/' fov '/' mouseid '_' expname ...
            '/group_proc/imreg_normcorre_CaImAn/' mouseid '_' expname '_imreg_ref' reffile '/'];
else
    dir_exp = [data_locn 'Analysis/' mouseid '/' mouseid '_' expname ...
            '/group_proc/imreg_normcorre_CaImAn/' mouseid '_' expname '_imreg_ref' reffile '/'];
end

fpf = load([dir_exp mouseid '_' expname '_ref' reffile '_framesperfile.mat']);
if dostep(1)
    so = load([dir_exp mouseid '_' expname '_ref' reffile '_segment_output.mat']);
end
if dostep(2)
    fissa = load([dir_exp str_fissa '/' mouseid '_' expname '_ref' reffile '_fissa_output.mat']);
end
if dostep(3)
    spikes = load([dir_exp str_fissa '/' mouseid '_' expname '_ref' reffile '_spikes.mat']);
end

% copy experiment files to individual run folders
ff = 1; a1 = 1; a2 = 1; a3 = 1; 
        b1 = 0; b2 = 0; b3 = 0;
for r = 1:numel(run)
    dir_run = [data_locn 'Analysis/' mouseid '/' fov '/' mouseid '_' expname '-' run{r}...
                '/group_proc/imreg_normcorre_CaImAn/' mouseid '_' expname '-' run{r} ...
                '_imreg_ref' reffile '_concrunsrois/'];
    if ~exist(dir_run,'dir') 
        mkdir( dir_run ); fileattrib(dir_run,'+w','g','s'); 
    end
    
    framesperfile = fpf.framesperfile;
    framesperfile = framesperfile(ff:ff+numfiles(r)-1);
    ff = ff + numfiles(r);
    
    % copy framesperfile
    if any(force) || ~exist([dir_run mouseid '_' expname '-' run{r} '_ref' reffile '_framesperfile.mat'],'file')
        save([dir_run mouseid '_' expname '-' run{r} '_ref' reffile '_framesperfile.mat'],'framesperfile');
    end
    
    % copy roi segmentation output
    if dostep(1)
        if force(1) || ~exist([dir_run mouseid '_' expname '-' run{r} '_ref' reffile '_segment_output.mat'],'file')
            fprintf('Copying ROI segmentation files to %s.\n', [mouseid '_' expname '-' run{r} '_ref' reffile]);
            copyfile([dir_exp mouseid '_' expname '_ref' reffile '_elimROIs.fig'],...
                     [dir_run mouseid '_' expname '-' run{r} '_ref' reffile '_elimROIs.fig'])
            copyfile([dir_exp mouseid '_' expname '_ref' reffile '_elimROIs.png'],...
                     [dir_run mouseid '_' expname '-' run{r} '_ref' reffile '_elimROIs.png'])
            copyfile([dir_exp mouseid '_' expname '_ref' reffile '_imreg_template.fig'],...
                     [dir_run mouseid '_' expname '-' run{r} '_ref' reffile '_imreg_template.fig'])
            copyfile([dir_exp mouseid '_' expname '_ref' reffile '_imreg_template.mat'],...
                     [dir_run mouseid '_' expname '-' run{r} '_ref' reffile '_imreg_template.mat'])
            copyfile([dir_exp mouseid '_' expname '_ref' reffile '_imreg_template.png'],...
                     [dir_run mouseid '_' expname '-' run{r} '_ref' reffile '_imreg_template.png'])
            copyfile([dir_exp mouseid '_' expname '_ref' reffile '_ROIs.fig'],...
                     [dir_run mouseid '_' expname '-' run{r} '_ref' reffile '_ROIs.fig'])
            copyfile([dir_exp mouseid '_' expname '_ref' reffile '_ROIs.png'],...
                     [dir_run mouseid '_' expname '-' run{r} '_ref' reffile '_ROIs.png'])
            
            so2 = so;
            b1 = b1 + sum(framesperfile);
            so2.tsG = so.tsG(:,a1:b1);
            so2.elim_tsG = so.elim_tsG(:,a1:b1);
            so2.df_f = so.df_f(:,a1:b1);
            so2.elim_df_f = so.elim_df_f(:,a1:b1);
            fprintf('Dividing tsG and df and saving to %s.\n', [mouseid '_' expname '-' run{r} '_ref' reffile]');
            save([dir_run mouseid '_' expname '-' run{r} '_ref' reffile '_segment_output.mat'],'-struct','so2');
            
            % save new plots for raw timeseries and dF/F
            multiplot_ts(so2.tsG, [dir_run mouseid '_' expname '-' run{r} '_ref' reffile '_raw_timeseries'], 'Raw timeseries')
            multiplot_ts(so2.df_f, [dir_run mouseid '_' expname '-' run{r} '_ref' reffile '_df_f'], 'dF/F')
        else
            fprintf('ROI files in %s already exist, skipping file copying.\n', [mouseid '_' expname '-' run{r} '_ref' reffile]');
        end
    end
    
    % copy fissa output
    if dostep(2)
        if force(2) || ~exist([dir_run 'FISSA/' mouseid '_' expname '-' run{r} '_ref' reffile '_fissa_output.mat'],'file')
            if exist([dir_exp 'FISSA/'],'dir') && ~exist([dir_run 'FISSA/'],'dir') 
                mkdir( [dir_run 'FISSA/'] ); fileattrib(dir_run,'+w','g','s'); 
            end
            
            fissa2 = fissa;
            fissa2.dtsG = fissa.dtsG(:,a2:b2);
            fissa2.ddf_f = fissa.ddf_f(:,a2:b2);
            a2 = b2+1;
            fprintf('Dividing fissa timeseries and saving to %s.\n', [mouseid '_' expname '-' run{r} '_ref' reffile]');
            save([dir_run 'FISSA/' mouseid '_' expname '-' run{r} '_ref' reffile '_fissa_output.mat'],'-struct','fissa2');
            
            % save new plots for fissa output
            multiplot_ts(fissa2.dtsG, [dir_run 'FISSA/' mouseid '_' expname '-' run{r} '_ref' reffile '_fissa_timeseries'],...
                'Fissa-corrected raw timeseries')
            multiplot_ts(fissa2.ddf_f, [dir_run 'FISSA/' mouseid '_' expname '-' run{r} '_ref' reffile '_fissa_df_f'],...
                'Fissa-corrected dF/F')
        else
            fprintf('FISSA output in %s already exists, skipping file copying.\n', [mouseid '_' expname '-' run{r} '_ref' reffile]');
        end
    end
    
    % copy spikes
    if dostep(3)
        if force(3) || ~exist([dir_run str_fissa '/' mouseid '_' expname '-' run{r} '_ref' reffile '_spikes.mat'],'file')            
            spikes2 = spikes;
            spikes2.spikes = spikes2.spikes(:,a3:b3);
            a3 = b3+1;
            fprintf('Dividing spikes and saving to %s.\n', [mouseid '_' expname '-' run{r} '_ref' reffile]');
            save([dir_run str_fissa '/' mouseid '_' expname '-' run{r} '_ref' reffile '_spikes.mat'],'-struct','fissa2');
            
            % save new plot for spikes
            plotSpikes(spikes2.spikes, [dir_run str_fissa '/' mouseid '_' expname '-' run{r} '_ref' reffile '_spikes'])
        else
            fprintf('Spikes in %s already exists, skipping file copying.\n', [mouseid '_' expname '-' run{r} '_ref' reffile]');
        end
    end
end