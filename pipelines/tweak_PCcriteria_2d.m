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
fprintf('%s: %g total cells\n', [mouseid '_' expname], size(PFdata.spkMap,1))

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
        bin_phi, activespk, hist.infoMap, hist.pf_activet, PFdata.activetrials, prctile_thr, ...
        0, 0, Nrand, 'hist', 2 );
    fprintf('%s: %g pcs after bootstrap shuffle test\n', [mouseid '_' expname], numel(inc_idx1))
    if ~isempty(exc_idx1) && doplot
        disp('independent testing')
        plotSpikeRasterTrials( PFdata.normspkRaster(exc_idx1), PFdata.ytick_files, 'NonPC', false, [], false )
    end

    % PC selection through pf_activet criterion
    % pfactivet_thr = 0.04;
    inc_idx2 = []; exc_idx2 = []; 
    Ncells = size(PFdata.spkMap,1);
    for c = 1:Ncells
        if hist.pf_activet(c) >= pfactivet_thr
            inc_idx2 = [inc_idx2; c];
        else
            exc_idx2 = [exc_idx2; c];
        end
    end
    fprintf('%s: %g pcs after pf_activet criterion\n', [mouseid '_' expname], numel(inc_idx2))
    if ~isempty(exc_idx2) && doplot
        plotSpikeRasterTrials( PFdata.normspkRaster(exc_idx2), PFdata.ytick_files, 'NonPC', false, [], false )
    end
    if ~isempty(inc_idx2) && doplot
        plotSpikeRasterTrials( PFdata.normspkRaster(inc_idx2), PFdata.ytick_files, 'PC', false, [], false )
    end

    % PC selection through activetrials_thr criterion
    % activetrials_thr = 0.5;
    inc_idx3 = []; exc_idx3 = []; 
    Ncells = size(PFdata.spkMap,1);
    for c = 1:Ncells
        if PFdata.activetrials(c) >= activetrials_thr
            inc_idx3 = [inc_idx3; c];
        else
            exc_idx3 = [exc_idx3; c];
        end
    end
    fprintf('%s: %g pcs after activetrials_thr criterion\n', [mouseid '_' expname], numel(inc_idx3))
    if ~isempty(exc_idx3) && doplot
        plotSpikeRasterTrials( PFdata.normspkRaster(exc_idx3), PFdata.ytick_files, 'NonPC', false, [], false )
    end
    if ~isempty(inc_idx3) && doplot
        plotSpikeRasterTrials( PFdata.normspkRaster(inc_idx2), PFdata.ytick_files, 'PC', false, [], false )
    end
end

%% Sequential enforcement of PC selection criteria
% bootstrap shuffle test
[pcIdx_SIsec1, pcIdx_SIspk1, nonpcIdx_SIsec1, nonpcIdx_SIspk1] ...
    = identifyPCs_2d( activespk, xh, yh, hist.infoMap, hist.pf_activet, hist.Nbins, prctile_thr, 0, Nrand, 'hist', 1, gaussfiltSigma );
fprintf('%s: %g pcs after bootstrap shuffle test\n', [mouseid '_' expname], numel(pcIdx_SIsec1))
if ~isempty(pcIdx_SIsec1) && doplot
    plotpfMaps_2d(activeData, hist.rMap(:,:,pcIdx_SIsec1), hist.rMap_sm(:,:,pcIdx_SIsec1), [], pcIdx_SIsec1, 'PC');
end
if ~isempty(nonpcIdx_SIsec1) && doplot
    plotpfMaps_2d(activeData, hist.rMap(:,:,nonpcIdx_SIsec1), hist.rMap_sm(:,:,nonpcIdx_SIsec1), [], nonpcIdx_SIsec1, 'NonPC');
end

% PC selection through pf_activet criterion
pfactivet_thr = 0.01;
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

fprintf('%s: %g pcs after activetrials_thr criterion\n', [mouseid '_' expname], numel(pcIdx_SIsec2))
if ~isempty(pcIdx_SIsec2) && doplot
    plotpfMaps_2d(activeData, hist.rMap(:,:,pcIdx_SIsec2), hist.rMap_sm(:,:,pcIdx_SIsec2), [], pcIdx_SIsec2, 'PC');
end
if ~isempty(nonpcIdx_SIsec2) && doplot
    plotpfMaps_2d(activeData, hist.rMap(:,:,nonpcIdx_SIsec2), hist.rMap_sm(:,:,nonpcIdx_SIsec2), [], nonpcIdx_SIsec2, 'NonPC');
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
output.PFdata = PFdata;
output.params = params;
save([sdir mouseid '_' expname '_ref' reffile '_PFmap_output.mat'],'-struct','output');

fig_sdir = [sdir '/PFdata/'];
rmdir(fig_sdir,'s'); mkdir(fig_sdir);
plotPF_1d(hist, [], PFdata, true, true, fig_sdir, [mouseid '_' expname '_imreg_ref' reffile ])

