% Written by Ann Go
% 
% This script registers the ROIs from two environments and plots the place
% field tuning for each sorted according to their own cells and according
% to the cells of the other environment.
%
% INPUTS
% mouseid   : 'm##' e.g. 'm62'
% env1      : environment 1
% env2      : environment 2
%   e.g. 'fam1', 'fam2', 'nov', 'fam1rev'
% ref1      : image registration template file for env1
% ref2      : image registration template file for env2
%   format: 'YYYYMMDD_hh_mm_ss'
% force     : flag to force generation of comparison figures even though they
%               already exist
% figclose    : flag to close figures that will be generated (they are
%               automatically saved)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The section labeled "USER-DEFINED INPUT" requires user input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function frun_showRemapping_2env_multiAnimals( mouseid_array, env1, ref1_array, env2, ref2_array, force, figclose )

if nargin<6, force = false; end
if nargin<7, figclose = true; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USER-DEFINED INPUT                         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Basic settings
mcorr_method = 'normcorre';
segment_method = 'CaImAn';
dofissa = true;
    if dofissa, str_fissa = 'FISSA'; else, str_fissa = 'noFISSA'; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load module folders and define data directory
[data_locn,~,err] = load_neuroSEEmodules(false);
if ~isempty(err)
    beep
    cprintf('Errors',err);    
    return
end

tic
% Output directory
mousegrp = '';
for n = 1:size(mouseid_array,1)
    if ~isempty(mousegrp)
        mousegrp = [mousegrp '_' mouseid_array(n,2:end)];
    else
        mousegrp = mouseid_array(n,:);
    end
end
sdir = [data_locn 'Analysis/remapping_multiAnimals/' env1 env2 '/' mousegrp '/'];

fname_remap = [sdir mousegrp '_' env1 env2 '_remapping_output.mat'];
fname_remapfig = [sdir mousegrp '_' env1 env2 '_remapping_summary.fig'];

if ~exist(fname_remap,'file') || force 
    % Pooling remapping data
    env1PF = [];
    env2PF_env1Sorting = [];
    env1PF_env2Sorting = [];
    env2PF = [];
    for n = 1:size(mouseid_array,1)
        % mouse identity
        mouseid = mouseid_array(n,:);

        fdir = [data_locn 'Analysis/' mouseid '/' mouseid '_' env1 env2 '/remapping/imreg_' mcorr_method '_' ...
               segment_method '_' str_fissa '/' mouseid '_' env1 env2 '_imreg_ref' ref1_array(n,:) '-' ref2_array(n,:) '/'];
        
        % Check if data exist for mouse in env1 and env2. Quit if data does not exist
        fname = [fdir  mouseid '_' env1 env2 '_remapping_output.mat'];
        if ~exist(fname,'file')   
            beep
            err = sprintf('No processed data for %s in %s & %s\n', mouseid, env1, env2);
            cprintf('Errors',err);    
            return
        end

        c = load(fname);
        env1PF = [env1PF; c.env1PF];
        env2PF_env1Sorting = [env2PF_env1Sorting; c.env2PF_env1Sorting];
        env1PF_env2Sorting = [env1PF_env2Sorting; c.env1PF_env2Sorting];
        env2PF = [env2PF; c.env2PF];
    end

    [ ~, maxLoc_env1 ] = prefLoc( env1PF );
    [ ~, maxLoc_env2 ] = prefLoc( env2PF );

    [ ~, sortIdx_env1 ] = sort( maxLoc_env1 );
    [ ~, sortIdx_env2 ] = sort( maxLoc_env2 );

    sort_env1PF = env1PF(sortIdx_env1,:);
    sort_env2PF_env1Sorting = env2PF_env1Sorting(sortIdx_env1,:);
    sort_env1PF_env2Sorting = env1PF_env2Sorting(sortIdx_env2,:);
    sort_env2PF = env2PF(sortIdx_env2,:);
    
    % Save output
    if ~exist(sdir,'dir'), mkdir(sdir); end
    output.sort_env1PF = sort_env1PF;
    output.sort_env2PF_env1Sorting = sort_env2PF_env1Sorting;
    output.sort_env1PF_env2Sorting = sort_env1PF_env2Sorting;
    output.sort_env2PF = sort_env2PF;
    fprintf('%s: saving remapping output\n',[mousegrp '_' env1 env2]);
    save(fname_remap,'-struct','output');
else
    if ~exist(fname_remapfig,'file')
        c = load(fname_remap);
        sort_env1PF = c.sort_env1PF;
        sort_env2PF_env1Sorting = c.sort_env2PF_env1Sorting;
        sort_env1PF_env2Sorting = c.sort_env1PF_env2Sorting;
        sort_env2PF = c.sort_env2PF;
    end
end

if ~exist(fname_remapfig,'file') || force
    fh = figure;
    fontsize = 16;
    Nbins = size(sort_env1PF,2);
    subplot(1,4,1);
        cmap = viridisMap;
        imagesc(sort_env1PF); 
        colormap(cmap); %colorbar
        title(env1,'Fontweight','normal','Fontsize',fontsize);
        yticks([1 size(sort_env1PF,1)]); yticklabels([1 size(sort_env1PF,1)]);
        xticks([1 Nbins]); xticklabels([1 100]);
        xlabel('Position (cm)'); %ylabel('Cell no.');
    subplot(1,4,2);
        imagesc(sort_env2PF_env1Sorting); 
        title(env2,'Fontweight','normal','Fontsize',fontsize);
        yticks([]);
        xticks([1 Nbins]); xticklabels([1 100]);
        xlabel('Position (cm)');
    subplot(1,4,3);
        imagesc(sort_env1PF_env2Sorting); 
        title(env1,'Fontweight','normal','Fontsize',fontsize); 
        yticks([1 size(sort_env2PF,1)]); yticklabels([1 size(sort_env2PF,1)]);
        xticks([1 Nbins]); xticklabels([1 100]);
        xlabel('Position (cm)');
    subplot(1,4,4);
        imagesc(sort_env2PF); 
        title(env2,'Fontweight','normal','Fontsize',fontsize); 
        yticks([]);
        xticks([1 Nbins]); xticklabels([1 100]);
        xlabel('Position (cm)');

        
    fprintf('%s: saving remapping summary figure\n',[mousegrp '_' env1 env2]);
    savefig( fh, fname_remapfig(1:end-4) );
    saveas( fh, fname_remapfig(1:end-4), 'png' );
    if figclose, close( fh ); end   
end

    
t = toc;
str = sprintf('%s: Processing done in %g hrs\n', [mousegrp '_' env1 env2], round(t/3600,2));
cprintf(str)
end



