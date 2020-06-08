% Written by Ann Go (some parts adapted from Giuseppe's PF_ASD_1d.m)
%
% This function maps place fields
%
% INPUTS:
%   spikes      : spike estimates obtained with oasisAR2
%   imtime      : imaging timestamps
%   downTrackdata   : cell of tracking data with fields x, y, r, phi, w,
%                  speed, time, alpha, TTLout
%   params.
%     fr                    : imaging frame rate [default: 30.9 Hz]
%     PFmap.Nbins           : number of location bins
%     PFmap.Nepochs         : number of epochs for each 4 min video [default: 1]
%     PFmap.Vthr            : speed threshold (mm/s) [default: 20]
%     PFmap.histsmoothFac   : Gaussian smoothing window for histogram
%                               estimation [default: 10]

% OUTPUTS:
%   occMap                  : occupancy map
%   hist., asd.
%     spkMap                : spike map (Ncells rows x Nbins columns)
%     normspkMap
%     infoMap               : information map 
%     pfMap                 : place field map obtained with histogram estimation 
%     pfMap_sm              : (hist only) smoothed version of placeMap 
%     normpfMap             : place field map obtained with histogram estimation 
%     normpfMap_sm          : (hist only) smoothed version of placeMap 
%     spkRaster
%     normspkRaster
%     pcIdx                 : row indices of original spikes corresponding
%                               to place cells
%   downData    : tracking data downsampled to imaging frequency, fields are
%                 x, y, r, phi, speed, t
%   activeData  : downsampled tracking data for when animal was moving, fields are
%                 x, y, r, phi, speed, t, spikes, spikes_pc 

function [ hist, asd, activeData, PFdata, varargout ] = generatePFmap_1d( spikes, downTrackdata, params, doasd )
    
if nargin < 4, doasd = false; end
fr = params.fr;
Nrand = params.PFmap.Nrand;
Nbins = params.PFmap.Nbins;
Nepochs = params.PFmap.Nepochs;
Vthr = params.PFmap.Vthr;
histsmoothWin = params.PFmap.histsmoothWin;
prctile_thr = params.PFmap.prctile_thr;
Ncells = size(spikes,1);

% Input data
x = downTrackdata.x;
y = downTrackdata.y;
phi = downTrackdata.phi;
r = downTrackdata.r;
speed = downTrackdata.speed;
t = downTrackdata.time;
    % find out how many files data is from
    ind = find(abs(diff(t))>200);
    if numel(ind)>1 % more than 1 image file
        dt = mean(diff(t(1:ind(1))));
    else
        dt = mean(diff(t));
    end

% Consider only samples when the mouse is active
activex     = x(speed > Vthr);
activey     = y(speed > Vthr);
activephi   = phi(speed > Vthr);
activespk   = spikes(:,speed > Vthr);
activet     = t(speed > Vthr);
activespeed = speed(speed > Vthr);
activer     = r(speed > Vthr);
clear x y phi r speed t

% Bin phi data
[bin_phi,~] = discretize(activephi,Nbins);


%% ALL CELLS: calculate PF data for entire duration (one epoch)
[PFdata, hist, asd] = calcPFdata_1d(bin_phi, activephi, activespk, activet, Nbins, fr);
% PFdata is a structure with the following fields
%   phi_trials, spk_trials, bintime_trials, bintime
%   spkRaster, normspkRaster, activetrials
%   spk_rate, spk_amplrate, bin_activity
%   occMap, spkMap = PFdata.spkMap, normspkMap = PFdata.normspkMap;

% hist and asd are structures with the following fields
%   rMap, rMap_sm (hist only), normrMap_sm (hist only)
%   infoMap, pfLoc, fieldSize, pfBins


%% Identify PLACE CELLS
% Cells are sorted in descending order of info content
[hist.SIsec.pcIdx, hist.SIspk.pcIdx, hist.SIsec.nonpcIdx, hist.SIspk.nonpcIdx] ...
    = identifyPCs_1d( bin_phi, activespk, hist.infoMap, prctile_thr, Nrand );
if doasd
    [asd.SIsec.pcIdx, asd.SIspk.pcIdx, asd.SIsec.nonpcIdx, asd.SIspk.nonpcIdx] ...
    = identifyPCs_1d( bin_phi, activespk, asd.infoMap,  prctile_thr, Nrand, 'asd');
end


%% Organise pc and nonpc data structures
% histogram estimation
hist.SIsec = organisefields(spkRaster, normspkRaster, spkMap, normspkMap, ...
                            spk_rate, spk_amplrate, hist, hist.SIsec, 'SIbitspersec');
hist.SIspk = organisefields(spkRaster, normspkRaster, spkMap, normspkMap, ...
                            spk_rate, spk_amplrate, hist, hist.SIspk, 'SIbitsperspk');

% asd
if doasd
    asd.SIsec = organisefields(spkRaster, normspkRaster, spkMap, normspkMap, ...
                                spk_rate, spk_amplrate, asd, asd.SIsec, 'SIbitspersec');
    asd.SIspk = organisefields(spkRaster, normspkRaster, spkMap, normspkMap, ...
                                spk_rate, spk_amplrate, asd, asd.SIspk, 'SIbitsperspk');
end

% sort pf maps
hist.SIsec = sortPFmaps(hist.SIsec);
hist.SIspk = sortPFmaps(hist.SIspk);


%% ADDITIONAL PROCESSING IF Nepochs > 1
% Calculate place field maps for each epoch 
if Nepochs > 1
    % Initialise matrices
    bintime_e = zeros(Nepochs, Nbins);                         
    occMap_e = zeros(Nepochs, Nbins);                         
    spkMap_e = zeros(Ncells, Nbins, Nepochs);
    normspkMap_e = zeros(Ncells, Nbins, Nepochs);   
    bin_phi_e = zeros(Nepochs, Nbins);
    activephi_e = zeros(Ncells, Nbins, Nepochs);
    activespk_e = zeros(Ncells, Nbins, Nepochs);
    activet_e = zeros(Ncells, Nbins, Nepochs);
    spk_rate_e = zeros(Ncells, Nepochs);
    spk_amplrate_e = zeros(Ncells, Nepochs);
    bin_activity_e = zeros(Ncells, Nbins, Nepochs);
    
    hist_epochs.rMap = zeros(Ncells, Nbins, Nepochs);               
    hist_epochs.rMap_sm = zeros(Ncells, Nbins, Nepochs);            
    hist_epochs.normrMap_sm = zeros(Ncells, Nbins, Nepochs);     
    hist_epochs.infoMap = zeros(Ncells, 2, Nepochs); 
    hist_epochs.pfLoc = zeros(Ncells, Nepochs);
    hist_epochs.fieldSize = zeros(Ncells, Nepochs);
    hist_epochs.pfBins = zeros(Ncells, Nepochs);
    
    if doasd
        asd_epochs.rMap = zeros(Ncells, Nbins, Nepochs);               
        asd_epochs.rMap_sm = zeros(Ncells, Nbins, Nepochs);            
        asd_epochs.normrMap_sm = zeros(Ncells, Nbins, Nepochs);     
        asd_epochs.infoMap = zeros(Ncells, 2, Nepochs); 
        asd_epochs.pfLoc = zeros(Ncells, Nepochs);
        asd_epochs.fieldSize = zeros(Ncells, Nepochs);
        asd_epochs.pfBins = zeros(Ncells, Nepochs);
    end
    
    % Calculate PF maps
    e_bound = round( linspace(1,size(activespk,2),Nepochs+1) );
    
    % separate exploration into smaller intervals
    for e = 1:Nepochs
        bin_phi_e(e,:) = bin_phi(e_bound(e):e_bound(e+1));
        activephi_e(:,:,e) = activephi(:,e_bound(e):e_bound(e+1));
        activespk_e(:,:,e) = activespk(:,e_bound(e):e_bound(e+1));
        activet_e(:,:,e) = activet(:,e_bound(e):e_bound(e+1));
            
        [PFdata_e, hist_e, asd_e] = calcPFdata_1d(bin_phi_e(e,:), activephi_e(:,:,e), activespk_e(:,:,e), activet_e(:,:,e), Nbins, fr);
        phi_trials_e{e} = PFdata_e.phi_trials;
        spk_trials_e{e} = PFdata_e.spk_trials;
        bintime_trials{e} = PFdata_e.bintime_trials;
        bintime_e(e,:) = PFdata_e.bintime;
        spkRaster = PFdata_e.spkRaster;
        normspkRaster = PFdata_e.normspkRaster;
        activetrials = PFdata_e.activetrials;
        spk_rate = PFdata_e.spk_rate;
        spk_amplrate = PFdata_e.spk_amplrate;
        bin_activity = PFdata_e.bin_activity;
        occMap = PFdata_e.occMap;
        spkMap = PFdata_e.spkMap;
        normspkMap = PFdata_e.normspkMap;
    end
end


%% Organise pc and nonpc data structures for Nepochs > 1
if Nepochs > 1
    % histogram estimation
    hist_e.SIsec = organisefields([], [], spkMap_e, normspkMap_e, ...
                            spk_rate_e, spk_amplrate_e, hist_e, hist_e.SIsec, 'SIbitspersec');
    hist_e.SIspk = organisefields([], [], spkMap_e, normspkMap_e, ...
                            spk_rate_e, spk_amplrate_e, hist_e, hist_e.SIspk, 'SIbitsperspk');

    %asd
    if doasd
        asd_e.SIsec = organisefields([], [], spkMap_e, normspkMap_e, ...
                                spk_rate_e, spk_amplrate_e, asd_e, asd_e.SIsec, 'SIbitspersec');
        asd_e.SIspk = organisefields([], [], spkMap_e, normspkMap_e, ...
                                spk_rate_e, spk_amplrate_e, asd_e, asd_e.SIspk, 'SIbitsperspk');
    end
end

% sort pf maps
hist_e.SIsec = sortPFmaps(hist_e.SIsec);
hist_e.SIspk = sortPFmaps(hist_e.SIspk);


%% Outputs
activeData.x = activex;
activeData.y = activey;
activeData.r = activer;
activeData.phi = activephi;
activeData.speed = activespeed;
activeData.time = activet;
activeData.spikes = activespk;

if Nepochs>1
    PFdata_epochs.occMap = occMap_e;                         
    PFdata_epochs.spkMap = spkMap_e;
    PFdata_epochs.normspkMap = normspkMap_e;   
    PFdata_epochs.bin_phi = bin_phi_e;
    PFdata_epochs.activespk = activespk_e;
    PFdata_epochs.spk_rate = spk_rate_e;
    PFdata_epochs.spk_amplrate = spk_amplrate_e;
    PFdata_epochs.bin_activity = bin_activity_e;
    
    hist_epochs = hist_e;
    if doasd, asd_epochs = asd_e; end
else
    hist_epochs = [];
    if doasd, asd_epochs = []; end
end

if ~doasd, asd = []; end

varargout(1) = hist_epochs;
varargout(2) = asd_epochs;
end

function hSIstruct = organisefields(spkRaster, normspkRaster, spkMap, normspkMap, ...
                                    bintime_trials, bin_activity, activetrials, ...
                                    spk_rate, spk_amplrate, hstruct, hSIstruct, infotype)
    if nargin<12, infotype = 2; end
    if strcmpi(infotype, 'SIbitspersec'), info = 1; end
    if strcmpi(infotype, 'SIbitsperspk'), info = 2; end
    
    Nepochs = size(spkmap,3);
    if Nepochs == 1
        pcIdx = hSIstruct.pcIdx;
        if ~isempty(pcIdx)
            hSIstruct.spkRaster_pc = spkRaster(pcIdx);
            hSIstruct.normspkRaster_pc = normspkRaster(pcIdx);
            hSIstruct.activetrials_pc = activetrials(pcIdx);
            hSIstruct.bintime_trials_pc = bintime_trials(pcIdx,:);
            hSIstruct.bin_activity_pc = bin_activity(pcIdx,:);
            hSIstruct.spk_rate_pc = spk_rate(pcIdx);
            hSIstruct.spk_amplrate_pc = spk_amplrate(pcIdx);
            hSIstruct.infoMap_pc = hstruct.infoMap(pcIdx, info);

            hSIstruct.spkMap = spkMap(pcIdx,:);
            hSIstruct.normspkMap = normspkMap(pcIdx,:);
            hSIstruct.pfMap = hstruct.rMap(pcIdx,:);
            if isfield(hstruct,'rMap_sm')
                hSIstruct.pfMap_sm = hstruct.rMap_sm(pcIdx,:);
                hSIstruct.normpfMap_sm = hstruct.normrMap_sm(pcIdx,:);
            end
            hSIstruct.pfLoc = hstruct.pfLoc(pcIdx);
            hSIstruct.pfSize = hstruct.fieldSize(pcIdx);
        end
        nonpcIdx = hSIstruct.nonpcIdx;
        if ~isempty(nonpcIdx)
            hSIstruct.spkRaster_nonpc = spkRaster(nonpcIdx);
            hSIstruct.normspkRaster_nonpc = normspkRaster(nonpcIdx);
            hSIstruct.activetrials_nonpc = activetrials(nonpcIdx);
            hSIstruct.bintime_trials_nonpc = bintime_trials(nonpcIdx,:);
            hSIstruct.bin_activity_nonpc = bin_activity(nonpcIdx,:);
            hSIstruct.spk_rate_nonpc = spk_rate(nonpcIdx);
            hSIstruct.spk_amplrate_nonpc = spk_amplrate(nonpcIdx);
            hSIstruct.infoMap_nonpc = hstruct.infoMap(nonpcIdx, info);
        end
    else
        for e = 1:Nepochs
            pcIdx = hSIstruct.pcIdx{e};
            if ~isempty(pcIdx)
                hSIstruct.spk_rate_pc{e} = spk_rate(pcIdx,e);
                hSIstruct.spk_amplrate_pc{e} = spk_amplrate(pcIdx,e);
                hSIstruct.infoMap_pc{e} = hstruct.infoMap(pcIdx,info,e);

                hSIstruct.spkMap{e} = spkMap(pcIdx,:,e);
                hSIstruct.normspkMap{e} = normspkMap(pcIdx,:,e);
                hSIstruct.pfMap{e} = hstruct.rMap(pcIdx,:,e);
                if isfield(hstruct,'rMap_sm')
                    hSIstruct.pfMap_sm{e} = hstruct.rMap_sm(pcIdx,:,e);
                    hSIstruct.normpfMap_sm{e} = hstruct.normrMap_sm(pcIdx,:,e);
                end
                hSIstruct.pfLoc{e} = hstruct.pfLoc(pcIdx,e);
                hSIstruct.pfSize{e} = hstruct.fieldSize(pcIdx,e);
            end
            nonpcIdx = hSIstruct.nonpcIdx{e};
            if ~isempty(nonpcIdx)
                hSIstruct.spk_rate_nonpc{e} = spk_rate(nonpcIdx,e);
                hSIstruct.spk_amplrate_nonpc{e} = spk_amplrate(nonpcIdx,e);
                hSIstruct.infoMap_nonpc{e} = hstruct.infoMap(nonpcIdx,info,e);
            end
        end
    end
end

function hSIstruct = sortPFmaps(hSIstruct)
    % Sort place field maps
    Nepochs = size(hSIstruct.pfMap,3);
    if Nepochs == 1
        if ~isempty(hSIstruct.pcIdx)
            [ ~, sortIdx ] = sort( hSIstruct.pfLoc );
            hSIstruct.sortIdx = sortIdx;
            hSIstruct.sort_pfMap = hSIstruct.pfMap(sortIdx,:);
            if isfield(hSIstruct,'pfMap_sm')
                hSIstruct.sort_pfMap_sm = hSIstruct.pfMap_sm(sortIdx,:);
                hSIstruct.sort_normpfMap_sm = hSIstruct.normpfMap_sm(sortIdx,:);
            end
        end
    else
        for e = 1:Nepochs
            if ~isempty(hSIstruct.pcIdx{e})
                [ ~, sortIdx ] = sort( hSIstruct.pfLoc{e} );
                hSIstruct.sortIdx{e} = sortIdx;
                hSIstruct.sort_pfMap{e} = hSIstruct.pfMap{e}(sortIdx,:);
                if isfield(hSIstruct,'pfMap_sm')
                    hSIstruct.sort_pfMap_sm{e} = hSIstruct.pfMap_sm{e}(sortIdx,:);
                    hSIstruct.sort_normpfMap_sm{e} = hSIstruct.normpfMap_sm{e}(sortIdx,:);
                end
            end
        end
    end
end

