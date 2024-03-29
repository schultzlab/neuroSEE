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

function [ hist, asd, PFdata, hist_epochs, asd_epochs, PFdata_epochs ] = generatePFmap_1d( spikes, downTrackdata, params )
if nargin<3, params = neuroSEE_setparams; end

doasd = params.methods.doasd; 
fr = params.PFmap.fr;
Nrand = params.PFmap.Nrand;
Nbins = params.PFmap.Nbins;
Nepochs = params.PFmap.Nepochs;
Vthr = params.PFmap.Vthr;
histsmoothWin = params.PFmap.histsmoothWin;
prctile_thr = params.PFmap.prctile_thr;
pfactivet_thr = params.PFmap.pfactivet_thr;
activetrials_thr = params.PFmap.activetrials_thr;
fieldrate_thr = params.PFmap.fieldrate_thr;

% Input data
% x = downTrackdata.x;
% y = downTrackdata.y;
phi = downTrackdata.phi;
% r = downTrackdata.r;
speed = downTrackdata.speed;
t = downTrackdata.time;

% Consider only samples when the mouse is active
% activex     = x(speed > Vthr);
% activey     = y(speed > Vthr);
activephi   = phi(speed > Vthr);
activespk   = spikes(:,speed > Vthr);
activet     = t(speed > Vthr);
% activespeed = speed(speed > Vthr);
% activer     = r(speed > Vthr);
clear x y phi r speed t

% Bin phi data
[bin_phi,~] = discretize(activephi,Nbins);


%% ALL CELLS: calculate PF data for entire duration (one epoch)
[PFdata, hist, asd] = calcPFdata_1d(bin_phi, activephi, activespk, activet, Nbins, histsmoothWin, fr, doasd);
% PFdata is a structure with the following fields
%   phi_trials, spk_trials, bintime_trials, bintime
%   spkRaster, normspkRaster, activetrials
%   spk_rate, spk_amplrate, bin_activity
%   occMap, spkMap, normspkMap

% hist and asd are structures with the following fields
%   rMap, rMap_sm (hist only), normrMap_sm (hist only)
%   infoMap, pfLoc, fieldSize, pfBins


%% Identify PLACE CELLS
% Cells are sorted in descending order of info content
[hist.SIsec.pcIdx, hist.SIspk.pcIdx, hist.SIsec.nonpcIdx, hist.SIspk.nonpcIdx] = identifyPCs_1d( ...
    bin_phi, activespk, hist.infoMap, hist.pf_activet, PFdata.activetrials, hist.pfBins, prctile_thr, pfactivet_thr, activetrials_thr, fieldrate_thr, Nrand );
if doasd
    [asd.SIsec.pcIdx, asd.SIspk.pcIdx, asd.SIsec.nonpcIdx, asd.SIspk.nonpcIdx] = identifyPCs_1d( ...
    bin_phi, activespk, asd.infoMap, asd.pf_activet, PFdata.activetrials, asd.pfBins, prctile_thr, pfactivet_thr, activetrials_thr, fieldrate_thr, Nrand, 'asd');
end

% sort pf maps
hist.SIsec = sortPFmaps(hist.rateMap, hist.rateMap_sm, hist.normrateMap_sm, hist.pfLoc, hist.SIsec);
hist.SIspk = sortPFmaps(hist.rateMap, hist.rateMap_sm, hist.normrateMap_sm, hist.pfLoc, hist.SIspk);

if doasd
    asd.SIsec = sortPFmaps(asd.rateMap, [], asd.normrateMap_sm, asd.pfLoc, asd.SIsec);
    asd.SIspk = sortPFmaps(asd.rateMap, [], asd.normrateMap_sm, asd.pfLoc, asd.SIspk);
end


%% ADDITIONAL PROCESSING IF Nepochs > 1
% Calculate place field maps for each epoch 
if Nepochs > 1
    
    % divide timeseries into equally spaced epochs
    % e_bound = round( linspace(1,size(activespk,2),Nepochs+1) );
    
    % divide timeseries into equal-number-laps epochs
    Ntrials = PFdata.ytick_files(end);
    tr_div = round( linspace(1,Ntrials,Nepochs+1) );  
    e_bound = zeros(length(tr_div),1);
    for i = 1:length(tr_div)
        e_bound(i) = PFdata.idx_trials{tr_div(i)}(1);
    end
    
    % separate exploration into smaller intervals
    activephi_e = cell(Nepochs,1);
    bin_phi_e = cell(Nepochs,1);
    activespk_e = cell(Nepochs,1);
    activet_e = cell(Nepochs,1);
    PFdata_e = cell(Nepochs,1);
    hist_e = cell(Nepochs,1);
    asd_e = cell(Nepochs,1);
    for e = 1:Nepochs
        activephi_e{e} = activephi(e_bound(e):e_bound(e+1));
        bin_phi_e{e} = bin_phi(e_bound(e):e_bound(e+1));
        activespk_e{e} = activespk(:,e_bound(e):e_bound(e+1));
        activet_e{e} = activet(e_bound(e):e_bound(e+1));
            
        [PFdata_e{e}, hist_e{e}, asd_e{e}] = calcPFdata_1d(bin_phi_e{e}, activephi_e{e}, activespk_e{e}, activet_e{e}, Nbins, histsmoothWin, fr, doasd);
    end
    
    % combine epoch data into one structure
    f = fields(PFdata_e{1});
    for i = 1:length(f)
        for e = 1:Nepochs
            PFdata_epochs.(f{i}){e} = PFdata_e{e}.(f{i});
        end
    end
    f = fields(hist_e{1});
    for i = 1:length(f)
        for e = 1:Nepochs
            hist_epochs.(f{i}){e} = hist_e{e}.(f{i});
        end
    end
    if doasd
        f = fields(asd_e{1});
        for i = 1:length(f)
            for e = 1:Nepochs
                asd_epochs.(f{i}){e} = asd_e{e}.(f{i});
            end
        end
    end

    for e = 1:Nepochs
        % Identify PLACE CELLS per epoch
        [hist_epochs.SIsec.pcIdx{e}, hist_epochs.SIspk.pcIdx{e}, hist_epochs.SIsec.nonpcIdx{e}, hist_epochs.SIspk.nonpcIdx{e}] = identifyPCs_1d( ...
            bin_phi_e{e}, activespk_e{e}, hist_epochs.infoMap{e}, hist_epochs.pf_activet{e}, PFdata_epochs.activetrials{e}, hist_epochs.pfBins{e}, prctile_thr, pfactivet_thr, activetrials_thr, fieldrate_thr, Nrand );
        if doasd
            [asd_epochs.SIsec.pcIdx{e}, asd_epochs.SIspk.pcIdx{e}, asd_epochs.SIsec.nonpcIdx{e}, asd_epochs.SIspk.nonpcIdx{e}] = identifyPCs_1d( ...
            bin_phi_e{e}, activespk_e{e}, asd_epochs.infoMap{e}, asd_epochs.pf_activet{e}, PFdata_epochs.activetrials{e}, asd_epochs.pfBins{e}, prctile_thr, pfactivet_thr, activetrials_thr, fieldrate_thr, Nrand, 'asd');
        end
    end

    % sort pf maps per epoch
    hist_epochs.SIsec = sortPFmaps(hist_epochs.rateMap, hist_epochs.rateMap_sm, hist_epochs.normrateMap_sm, hist_epochs.pfLoc, hist_epochs.SIsec);
    hist_epochs.SIspk = sortPFmaps(hist_epochs.rateMap, hist_epochs.rateMap_sm, hist_epochs.normrateMap_sm, hist_epochs.pfLoc, hist_epochs.SIspk);

    if doasd
        asd_epochs.SIsec = sortPFmaps(asd_epochs.rateMap, [], asd_epochs.normrateMap_sm, asd_epochs.pfLoc, asd_epochs.SIsec);
        asd_epochs.SIspk = sortPFmaps(asd_epochs.rateMap, [], asd_epochs.normrateMap_sm, asd_epochs.pfLoc, asd_epochs.SIspk);
    end
end

%% Outputs
if Nepochs == 1
    hist_epochs = [];
    PFdata_epochs = [];
    asd_epochs = []; 
end
if ~doasd, asd = []; end

end
