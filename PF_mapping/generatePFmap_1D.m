% Written by Ann Go (some parts adapted from Giuseppe's PF_ASD_1d.m)
%
% This function maps place fields
%
% INPUTS:
%   spikes      : spike estimates obtained with oasisAR2
%   imtime      : imaging timestamps
%   trackData   : cell of tracking data with fields x, y, r, phi, w,
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
%     spkMap_pertrial
%     normspkfMap_pertrial
%     pcIdx                 : row indices of original spikes corresponding
%                               to place cells
%   downData    : tracking data downsampled to imaging frequency, fields are
%                 x, y, r, phi, speed, t
%   activeData  : downsampled tracking data for when animal was moving, fields are
%                 x, y, r, phi, speed, t, spikes, spikes_pc 

function [ occMap, hist, asd, activeData ] = generatePFmap_1D( spikes, trackData, params, dsample )
    
if nargin<5, dsample = false; end

fr = params.fr;
Nbins = params.PFmap.Nbins;
Nepochs = params.PFmap.Nepochs;
Vthr = params.PFmap.Vthr;
histsmoothFac = params.PFmap.histsmoothFac;
prctile_thr = params.PFmap.prctile_thr;
Ncells = size(spikes,1);

%% Pre-process tracking data
if dsample
    downData = downsample_trackData(trackData, spikes, fr);
    downx = downData.x;
    downy = downData.y;
    downphi = downData.phi;
    downr = downData.r;
    downspeed = downData.speed;
    t = downData.time;
else
    downx = trackData.x;
    downy = trackData.y;
    downphi = trackData.phi;
    downr = trackData.r;
    downspeed = trackData.speed;
    t = trackData.time;
end

% Consider only samples when the mouse is active
activex     = downx(downspeed > Vthr);
activey     = downy(downspeed > Vthr);
activephi   = downphi(downspeed > Vthr);
activespk   = spikes(:,downspeed > Vthr);
activet     = t(downspeed > Vthr);
activespeed = downspeed(downspeed > Vthr);
activer     = downr(downspeed > Vthr);

% Bin phi data
[bin_phi,~] = discretize(activephi,Nbins);


%% Identify place cells by first calculating PF maps for entire session
% (i.e. Nepochs = 1)

% Initialise matrices
spkMap = zeros(Ncells, Nbins);         % spike map
pfMap = zeros(Ncells, Nbins);            % place field map
pfMap_asd = zeros(Ncells, Nbins);        % place field map for asd
infoMap = zeros(Ncells, 2);              % mutual info
infoMap_asd = zeros(Ncells, 2);          % mutual info for asd

% Calculate PF maps
occMap = histcounts(bin_phi,Nbins);
for id = 1:Ncells
    z = activespk(id,:);

    % Spike rate maps
    for n = 1:Nbins
        spkMap(id,n) = sum(z(bin_phi == n));
    end

    % histogram estimation
    pfMap(id,:) = spkMap(id,:)./occMap;
    [infoMap(id,1), infoMap(id,2)] = infoMeasures(pfMap(id,:),occMap,0);

    % ASD estimation
    [pfMap_asd(id,:),~] = runASD_1d(bin_phi,z',Nbins);
    [infoMap_asd(id,1), infoMap_asd(id,2)] =...
        infoMeasures(pfMap_asd(id,:)',ones(Nbins,1),0);
end

% Identify place cells
[hist.pcIdx,asd.pcIdx] = filter_nonPC(bin_phi, activespk, infoMap, infoMap_asd, Nbins, prctile_thr);
activespk_hist = activespk(hist.pcIdx,:);
activespk_asd = activespk(asd.pcIdx,:);
Npcs = length(hist.pcIdx);
Npcs_asd = length(asd.pcIdx);


%% Finalise place field maps, recalculate if Nepochs > 1
if Nepochs == 1
    hist.spkMap = spkMap(hist.pcIdx,:);
    hist.pfMap = pfMap(hist.pcIdx,:);
    for id = 1:Npcs
        hist.normspkMap(id,:) = hist.spkMap(id,:)./max(hist.spkMap(id,:));
        hist.pfMap_sm(id,:) = smoothdata(hist.pfMap(id,:),'gaussian',Nbins/histsmoothFac);
        hist.normpfMap(id,:) = hist.pfMap(id,:)./max(hist.pfMap(id,:));
        hist.normpfMap_sm(id,:) = hist.pfMap_sm(id,:)./max(hist.pfMap_sm(id,:));
    end
    hist.infoMap = infoMap(hist.pcIdx,:);
    
    asd.spkMap = spkMap(asd.pcIdx,:);
    asd.pfMap = pfMap_asd(asd.pcIdx,:);
    for id = 1:Npcs_asd
        asd.normspkMap(id,:) = asd.spkMap(id,:)./max(asd.spkMap(id,:));
        asd.normpfMap(id,:) = asd.pfMap(id,:)./max(asd.pfMap(id,:));
    end
    asd.infoMap = infoMap_asd(asd.pcIdx,:);

else
    % Calculate PF maps for each epoch
    % Initialise matrices
    occMap = zeros(Nepochs, Nbins);                         
    hist.spkMap = zeros(Npcs, Nbins, Nepochs);            
    hist.normspkMap = zeros(Npcs, Nbins, Nepochs);            
    hist.pfMap = zeros(Npcs, Nbins, Nepochs);               
    hist.pfMap_sm = zeros(Npcs, Nbins, Nepochs);            
    hist.normpfMap = zeros(Npcs, Nbins, Nepochs);        
    hist.normpfMap_sm = zeros(Npcs, Nbins, Nepochs);     
    hist.infoMap = zeros(Npcs, 2, Nepochs);              
    asd.spkMap = zeros(Npcs_asd, Nbins, Nepochs);            
    asd.normspkMap = zeros(Npcs_asd, Nbins, Nepochs);            
    asd.pfMap = zeros(Npcs_asd, Nbins, Nepochs);              
    asd.normpfMap = zeros(Npcs_asd, Nbins, Nepochs);          
    asd.infoMap = zeros(Npcs_asd, 2, Nepochs);               

    % Calculate PF maps
    e_bound = round( linspace(1,size(activespk,2),Nepochs+1) );
    for id = 1:Npcs
        z = activespk_hist(id,:);

        % separate exploration in smaller intervals
        for e = 1:Nepochs
            bin_phi_e = bin_phi(e_bound(e):e_bound(e+1));
            spike_e = z(e_bound(e):e_bound(e+1));

            % Occupancy and spike rate maps
            occMap(e,:) = histcounts(bin_phi_e,Nbins);
            for n = 1:Nbins
                hist.spkMap(id,n,e) = sum(spike_e(bin_phi_e == n));
            end
            hist.normspkMap(id,:,e) = hist.spkMap(id,:,e)./max(hist.spkMap(id,:,e));
            
            % histogram estimation
            hist.pfMap(id,:,e) = hist.spkMap(id,:,e)./occMap(e,:);
            hist.pfMap(isnan(hist.pfMap)) = 0;
            hist.pfMap_sm(id,:,e) = smoothdata(hist.pfMap(id,:,e),'gaussian',Nbins/histsmoothFac);

            hist.normpfMap(id,:,e) = hist.pfMap(id,:,e)./max(hist.pfMap(id,:,e));
            hist.normpfMap_sm(id,:,e) = hist.pfMap_sm(id,:,e)./max(hist.pfMap_sm(id,:,e));
            [hist.infoMap(id,1,e), hist.infoMap(id,2,e)] = infoMeasures(hist.pfMap(id,:,e),occMap(e,:),0);
        end
    end
    for id = 1:Npcs_asd
        z = activespk_asd(id,:);

        % separate exploration in smaller intervals
        for e = 1:Nepochs
            bin_phi_e = bin_phi(e_bound(e):e_bound(e+1));
            spike_e = z(e_bound(e):e_bound(e+1));

            % Occupancy and spike rate maps
            for n = 1:Nbins
                asd.spkMap(id,n,e) = sum(spike_e(bin_phi_e == n));
            end
            asd.normspkMap(id,:,e) = asd.spkMap(id,:,e)./max(asd.spkMap(id,:,e));

            % asd estimation
            [asd.pfMap(id,:,e),~] = runASD_1d(bin_phi_e,(spike_e)',Nbins);
            asd.normpfMap(id,:,e) = asd.pfMap(id,:,e)./max(asd.pfMap(id,:,e));
            [asd.infomap(id,1,e), asd.infomap(id,2,e)] = ...
                infoMeasures(squeeze(asd.pfMap(id,:,e))',ones(Nbins,1),0);
        end
    end
end

% Calculate spike maps per trial
dthr = 10;
for ii = 1:Ncells
    % find the delineations for the video: find t = 0
    idx_file = find(diff(activet) < 0);
    idx_file = [0; idx_file; numel(activet)] +1; 
    p = bin_phi;
    s = activespk(ii,:);
    Ntrial = 1;
    ytick_files = 1;
    
    for jj = 1:numel(idx_file)-1
        % find the delineations per trial (i.e. loop)
        p_tr = p(idx_file(jj):idx_file(jj+1)-1);
        idx_tr = find( abs(diff(p_tr)) > dthr );
        for k = 1:numel(idx_tr)-1
            if idx_tr(k+1) == idx_tr(k)+1
                idx_tr(k) = 0;
            end
        end
        idx_tr = idx_tr( idx_tr > 0 );
        
        for k = 1:numel(idx_tr)-1
            phi{Ntrial} = p_tr(idx_tr(k)+1:idx_tr(k+1));
            spike{ii}{Ntrial} = s(idx_tr(k)+1:idx_tr(k+1));
            Ntrial = Ntrial + 1;
        end
        
        Ntrials(jj) = numel(idx_tr)-1;
        if jj == numel(idx_file)-1
            ytick_files = [ytick_files; sum(Ntrials(1:jj))];
        else
            ytick_files = [ytick_files; sum(Ntrials(1:jj))+1];
        end
    end
end

for ii = 1:Ncells
    for tr = 1:numel(phi)
        phi_tr = phi{tr};
        spike_tr = spike{ii}{tr};

        for n = 1:Nbins
            spkMap_pertrial{ii}(tr,n) = sum(spike_tr(phi_tr == n));
        end

        normspkMap_pertrial{ii}(tr,:) = spkMap_pertrial{ii}(tr,:)./max(spkMap_pertrial{ii}(tr,:));
        normspkMap_pertrial{ii}(isnan(normspkMap_pertrial{ii})) = 0;
    end
end

if ~isempty(hist.pcIdx)
    hist.spkMap_pertrial = spkMap_pertrial(hist.pcIdx);
    hist.normspkMap_pertrial = normspkMap_pertrial(hist.pcIdx);
end

if ~isempty(asd.pcIdx)
    asd.spkMap_pertrial = spkMap_pertrial(asd.pcIdx);
    asd.normspkMap_pertrial = normspkMap_pertrial(asd.pcIdx);
end


% Outputs
activeData.x = activex;
activeData.y = activey;
activeData.r = activer;
activeData.phi = activephi;
activeData.speed = activespeed;
activeData.t = activet;
activeData.spikes = activespk;
activeData.spikes_hist = activespk_hist;
activeData.spikes_asd = activespk_asd; 
activeData.spkMap_pertrial = spkMap_pertrial;
activeData.normspkMap_pertrial = normspkMap_pertrial;
activeData.ytick_files = ytick_files;

end


