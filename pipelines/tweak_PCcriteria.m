function tweak_PCcriteria( list, reffile, bl_prctile )

%% USER INPUT
% list = 'list_m62_fam1nov-fam1.txt';
% reffile = '20181015_09_37_54';
mcorr_method = 'normcorre';
segment_method = 'CaImAn';
dofissa = true;
groupreg_method = 'imreg';
% bl_prctile = 85;

% parameters to tweak for PC selection
independent_testing = false;
doplot = true;
pfactivet_thr = 0.05;      % fraction of dwell time in place field cell is required to be active (default: 0.05)
activetrials_thr = 0.5;    % fraction of trials cell is required to be active (default: 0.2)
prctile_thr = 99.99;
Nrand = 1000;


%% Load module folders and define data directory
[data_locn,~,err] = load_neuroSEEmodules;
if ~isempty(err)
    beep
    cprintf('Errors',err);    
    return
end

%% MouseID and experiment name
[ mouseid, expname ] = find_mouseIDexpname(list);

%% Location of downsampled tracking data for list
if dofissa, str_fissa = 'FISSA'; else, str_fissa = 'noFISSA'; end

fdir = [data_locn 'Analysis/' mouseid '/' mouseid '_' expname '/group_proc/'];

load([fdir mouseid '_' expname '_downTrackdata.mat'])    

%% Location of processed group data for list
if isempty(bl_prctile)
    sdir = [fdir groupreg_method '_' mcorr_method '_' segment_method '/'...
            mouseid '_' expname '_imreg_ref' reffile '/' str_fissa '/'];
else
    sdir = [fdir groupreg_method '_' mcorr_method '_' segment_method '/'...
            mouseid '_' expname '_imreg_ref' reffile '/' str_fissa '/bl_prctile' num2str(bl_prctile) '/'];
end

load([sdir mouseid '_' expname '_ref' reffile '_spikes.mat'])
load([sdir mouseid '_' expname '_ref' reffile '_PFmap_output.mat'])

%% Manually eliminate rois from ROI segmentation data
fprintf('%s: %g total cells\n', [mouseid '_' expname], size(pfData.spkMap,1))

% PC selection via bootstrap shuffle test
Vthr = params.Vthr;
Nbins = params.Nbins; 
activephi   = phi(speed > Vthr);
activespk   = spikes(:,speed > Vthr);
activet     = time(speed > Vthr);
[bin_phi,~] = discretize(activephi,Nbins);

if independent_testing
    % bootstrap shuffle test
    [ ~, inc_idx1, ~, exc_idx1 ] = identifyPCs_1d( ...
        bin_phi, activespk, hist.infoMap, hist.pf_activet, pfData.activetrials, prctile_thr, ...
        0, 0, Nrand, 'hist', 2 );
    fprintf('%s: %g pcs after bootstrap shuffle test\n', [mouseid '_' expname], numel(inc_idx1))
    if ~isempty(exc_idx1) && doplot
        disp('independent testing')
        plotSpikeRasterTrials( pfData.normspkRaster(exc_idx1), pfData.ytick_files, 'NonPC', false, [], false )
    end

    % PC selection through pf_activet criterion
    % pfactivet_thr = 0.04;
    inc_idx2 = []; exc_idx2 = []; 
    Ncells = size(pfData.spkMap,1);
    for c = 1:Ncells
        if hist.pf_activet(c) >= pfactivet_thr
            inc_idx2 = [inc_idx2; c];
        else
            exc_idx2 = [exc_idx2; c];
        end
    end
    fprintf('%s: %g pcs after pf_activet criterion\n', [mouseid '_' expname], numel(inc_idx2))
    if ~isempty(exc_idx2) && doplot
        plotSpikeRasterTrials( pfData.normspkRaster(exc_idx2), pfData.ytick_files, 'NonPC', false, [], false )
    end
    if ~isempty(inc_idx2) && doplot
        plotSpikeRasterTrials( pfData.normspkRaster(inc_idx2), pfData.ytick_files, 'PC', false, [], false )
    end

    % PC selection through activetrials_thr criterion
    % activetrials_thr = 0.5;
    inc_idx3 = []; exc_idx3 = []; 
    Ncells = size(pfData.spkMap,1);
    for c = 1:Ncells
        if pfData.activetrials(c) >= activetrials_thr
            inc_idx3 = [inc_idx3; c];
        else
            exc_idx3 = [exc_idx3; c];
        end
    end
    fprintf('%s: %g pcs after activetrials_thr criterion\n', [mouseid '_' expname], numel(inc_idx3))
    if ~isempty(exc_idx3) && doplot
        plotSpikeRasterTrials( pfData.normspkRaster(exc_idx3), pfData.ytick_files, 'NonPC', false, [], false )
    end
    if ~isempty(inc_idx3) && doplot
        plotSpikeRasterTrials( pfData.normspkRaster(inc_idx2), pfData.ytick_files, 'PC', false, [], false )
    end
end

%% Sequential enforcement of PC selection criteria
% bootstrap shuffle test
[ pcIdx_SIsec1, pcIdx_SIspk1, nonpcIdx_SIsec1, nonpcIdx_SIspk1 ] = identifyPCs_1d( ...
    bin_phi, activespk, hist.infoMap, hist.pf_activet, pfData.activetrials, prctile_thr, ...
    0, 0, Nrand, 'hist', 2 );
fprintf('%s: %g pcs after bootstrap shuffle test\n', [mouseid '_' expname], numel(pcIdx_SIspk1))
if ~isempty(nonpcIdx_SIspk1) && doplot
    plotSpikeRasterTrials( pfData.normspkRaster(nonpcIdx_SIspk1), pfData.ytick_files, 'NonPC', false, [], false )
end

% PC selection through pf_activet criterion
% pfactivet_thr = 0.05;
pcIdx_SIsec2 = []; 
nonpcIdx_SIsec2 = nonpcIdx_SIsec1; 
for n = 1:numel(pcIdx_SIsec1)
    c = pcIdx_SIsec1(n);
    if hist.pf_activet(c) >= pfactivet_thr
        pcIdx_SIsec2 = [pcIdx_SIsec2; c];
    else
        nonpcIdx_SIsec2 = [nonpcIdx_SIsec2; c];
    end
end

pcIdx_SIspk2 = []; 
nonpcIdx_SIspk2 = nonpcIdx_SIspk1; 
for n = 1:numel(pcIdx_SIspk1)
    c = pcIdx_SIspk1(n);
    if hist.pf_activet(c) >= pfactivet_thr
        pcIdx_SIspk2 = [pcIdx_SIspk2; c];
    else
        nonpcIdx_SIspk2 = [nonpcIdx_SIspk2; c];
    end
end
fprintf('%s: %g pcs after pf_activet criterion\n', [mouseid '_' expname], numel(pcIdx_SIspk2))
if ~isempty(nonpcIdx_SIspk2) && doplot
    plotSpikeRasterTrials( pfData.normspkRaster(nonpcIdx_SIspk2), pfData.ytick_files, 'NonPC', false, [], false )
end

% PC selection through activetrials_thr criterion
% activetrials_thr = 0.6;
pcIdx_SIsec3 = []; nonpcIdx_SIsec3 = nonpcIdx_SIsec2; 
for n = 1:numel(pcIdx_SIsec2)
    c = pcIdx_SIsec2(n);
    if pfData.activetrials(c) >= activetrials_thr
        pcIdx_SIsec3 = [pcIdx_SIsec3; c];
    else
        nonpcIdx_SIsec3 = [nonpcIdx_SIsec3; c];
    end
end

pcIdx_SIspk3 = []; nonpcIdx_SIspk3 = nonpcIdx_SIspk2; 
for n = 1:numel(pcIdx_SIspk2)
    c = pcIdx_SIspk2(n);
    if pfData.activetrials(c) >= activetrials_thr
        pcIdx_SIspk3 = [pcIdx_SIspk3; c];
    else
        nonpcIdx_SIspk3 = [nonpcIdx_SIspk3; c];
    end
end
fprintf('%s: %g pcs after activetrials_thr criterion\n', [mouseid '_' expname], numel(pcIdx_SIspk3))
if ~isempty(nonpcIdx_SIspk3) && doplot
    plotSpikeRasterTrials( pfData.normspkRaster(nonpcIdx_SIspk3), pfData.ytick_files, 'NonPC', false, [], false )
end
if ~isempty(pcIdx_SIspk3) && doplot
    plotSpikeRasterTrials( pfData.normspkRaster(pcIdx_SIspk3), pfData.ytick_files, 'PC', false, [], false )
end

hist.SIsec.pcIdx = pcIdx_SIsec3; hist.SIsec.nonpcIdx = nonpcIdx_SIsec3;
hist.SIsec = sortPFmaps(hist.rateMap, hist.rateMap_sm, hist.normrateMap_sm, hist.pfLoc, hist.SIsec);

hist.SIspk.pcIdx = pcIdx_SIspk3; hist.SIspk.nonpcIdx = nonpcIdx_SIspk3;
hist.SIspk = sortPFmaps(hist.rateMap, hist.rateMap_sm, hist.normrateMap_sm, hist.pfLoc, hist.SIspk);
if doplot    
    cmap0 = [0.9 0.9 0.9];
    cmap1 = [0 0 1];
    cmap = zeros(50,3);
    for j=1:3
        cmap(:,j) = linspace(cmap0(j),cmap1(j),50);
    end
    figure; imagesc(hist.SIspk.sort_normpfMap_sm); colormap(cmap);
end

% save data
output.hist = hist;
output.pfData = pfData;
output.params = params;
save([sdir mouseid '_' expname '_ref' reffile '_PFmap_output.mat'],'-struct','output');

fig_sdir = [sdir '/PFdata/'];
rmdir(fig_sdir,'s'); % mkdir(fig_sdir);
plotPF_1d(hist, [], pfData, true, true, sdir, [mouseid '_' expname '_imreg_ref' reffile ])

