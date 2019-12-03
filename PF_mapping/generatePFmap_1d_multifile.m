% Written by Ann Go 

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

function [ occMap, hist, asd, activeData ] = generatePFmap_1d_multifile( spikes, trackData, params )
    
Nbins = params.PFmap.Nbins;
Nepochs = params.PFmap.Nepochs;
Vthr = params.PFmap.Vthr;
histsmoothFac = params.PFmap.histsmoothFac;
prctile_thr = params.PFmap.prctile_thr;
Ncells = size(spikes,1);

%% Tracking data
downphi   = trackData.phi;
downx     = trackData.x;
downy     = trackData.y;
downspeed = trackData.speed;
downr     = trackData.r;
t         = trackData.time;

activex = cell(Ncells,1);       activey = cell(Ncells,1);
activephi = cell(Ncells,1);     activespk = cell(Ncells,1);
activet = cell(Ncells,1);       activespeed = cell(Ncells,1);
activer = cell(Ncells,1);       
bin_phi = cell(Ncells,1);       occMap = zeros(Ncells,Nbins);

for ii = 1:Ncells
    % Consider only samples when the mouse is active
    ind = find(downspeed{ii} > Vthr);
    activex{ii}     = downx{ii}(ind);
    activey{ii}     = downy{ii}(ind);
    activephi{ii}   = downphi{ii}(ind);
    activespk{ii}   = spikes{ii}(ind);
    activet{ii}     = t{ii}(ind);
    activespeed{ii} = downspeed{ii}(ind);
    activer{ii}     = downr{ii}(ind);
    
    % Bin phi data
    [bin_phi{ii},~] = discretize(activephi{ii},Nbins);
    occMap(ii,:) = histcounts(bin_phi{ii},Nbins);
end

%% Identify place cells by first calculating PF maps for entire session
% (i.e. Nepochs = 1)

% Initialise matrices
spkMap = zeros(Ncells, Nbins);         % spike map
pfMap = zeros(Ncells, Nbins);            % place field map
pfMap_asd = zeros(Ncells, Nbins);        % place field map for asd
infoMap = zeros(Ncells, 2);              % mutual info
infoMap_asd = zeros(Ncells, 2);          % mutual info for asd

% Calculate PF maps
for ii = 1:Ncells
    z = activespk{ii};

    % Spike rate maps
    for n = 1:Nbins
        spkMap(ii,n) = sum(z(bin_phi{ii} == n));
    end

    % histogram estimation
    pfMap(ii,:) = spkMap(ii,:)./occMap(ii,:);
    [infoMap(ii,1), infoMap(ii,2)] = infoMeasures(pfMap(ii,:),occMap(ii,:),0);

    % ASD estimation
    [pfMap_asd(ii,:),~] = runASD_1d(bin_phi{ii},z,Nbins);
    [infoMap_asd(ii,1), infoMap_asd(ii,2)] =...
        infoMeasures(pfMap_asd(ii,:)',ones(Nbins,1),0);
end

% Identify place cells
[hist.pcIdx,asd.pcIdx] = filter_nonPC(bin_phi, activespk, infoMap, infoMap_asd, Nbins, prctile_thr);
activespk_hist = activespk(hist.pcIdx,:);
activespk_asd = activespk(asd.pcIdx,:);
Npcs = length(hist.pcIdx);
Npcs_asd = length(asd.pcIdx);

% Calculate spike maps per trial
dthr = 10;
for ii = 1:Ncells
    % find the delineations for the video: find t = 0
    idx_file = find(diff(activet{ii}) < 0);
    idx_file = [0; idx_file; numel(activet{ii})] +1; 
    p = bin_phi{ii};
    s = activespk{ii};
    Ntrial = 1;
    ytick_files{ii} = 1;
    
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
            phi{ii}{Ntrial} = p_tr(idx_tr(k)+1:idx_tr(k+1));
            spike{ii}{Ntrial} = s(idx_tr(k)+1:idx_tr(k+1));
            Ntrial = Ntrial + 1;
        end
        
        Ntrials(ii,jj) = numel(idx_tr)-1;
        if jj == numel(idx_file)-1
            ytick_files{ii} = [ytick_files{ii}; sum(Ntrials(ii,1:jj))];
        else
            ytick_files{ii} = [ytick_files{ii}; sum(Ntrials(ii,1:jj))+1];
        end
    end
end

for ii = 1:Ncells
    for tr = 1:numel(phi{ii})
        phi_tr = phi{ii}{tr};
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
    hist.ytick_files = ytick_files(hist.pcIdx);
end

if ~isempty(asd.pcIdx)
    asd.spkMap_pertrial = spkMap_pertrial(asd.pcIdx);
    asd.normspkMap_pertrial = normspkMap_pertrial(asd.pcIdx);
    asd.ytick_files = ytick_files(asd.pcIdx);
end


%% Finalise place field maps, recalculate if Nepochs > 1
if Nepochs == 1
    hist.spkMap = spkMap(hist.pcIdx,:);
    hist.pfMap = pfMap(hist.pcIdx,:);
    for ii = 1:Npcs
        hist.normspkMap(ii,:) = hist.spkMap(ii,:)./max(hist.spkMap(ii,:));
        hist.pfMap_sm(ii,:) = smoothdata(hist.pfMap(ii,:),'gaussian',Nbins/histsmoothFac);
        hist.normpfMap(ii,:) = hist.pfMap(ii,:)./max(hist.pfMap(ii,:));
        hist.normpfMap_sm(ii,:) = hist.pfMap_sm(ii,:)./max(hist.pfMap_sm(ii,:));
    end
    hist.infoMap = infoMap(hist.pcIdx,:);
    
    asd.spkMap = spkMap(asd.pcIdx,:);
    asd.pfMap = pfMap_asd(asd.pcIdx,:);
    for ii = 1:Npcs_asd
        asd.normspkMap(ii,:) = asd.spkMap(ii,:)./max(asd.spkMap(ii,:));
        asd.normpfMap(ii,:) = asd.pfMap(ii,:)./max(asd.pfMap(ii,:));
    end
    asd.infoMap = infoMap_asd(asd.pcIdx,:);

else
    % Calculate PF maps for each epoch
    % Initialise matrices
    occMap = zeros(Npcs, Nbins, Nepochs);                         
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
    e_bound = round( linspace(1,size(activespk{1},1),Nepochs+1) );
    for ii = 1:Npcs
        z = activespk_hist{ii};

        % separate exploration in smaller intervals
        for e = 1:Nepochs
            bin_phi_e = bin_phi{ii}(e_bound(e):e_bound(e+1));
            spike_e = z(e_bound(e):e_bound(e+1));

            % Occupancy and spike rate maps
            occMap(ii,:,e) = histcounts(bin_phi_e,Nbins);
            for n = 1:Nbins
                hist.spkMap(ii,n,e) = sum(spike_e(bin_phi_e == n));
            end
            hist.normspkMap(ii,:,e) = hist.spkMap(ii,:,e)./max(hist.spkMap(ii,:,e));
            
            % histogram estimation
            hist.pfMap(ii,:,e) = hist.spkMap(ii,:,e)./occMap(e,:);
            hist.pfMap(isnan(hist.pfMap)) = 0;
            hist.pfMap_sm(ii,:,e) = smoothdata(hist.pfMap(ii,:,e),'gaussian',Nbins/histsmoothFac);

            hist.normpfMap(ii,:,e) = hist.pfMap(ii,:,e)./max(hist.pfMap(ii,:,e));
            hist.normpfMap_sm(ii,:,e) = hist.pfMap_sm(ii,:,e)./max(hist.pfMap_sm(ii,:,e));
            [hist.infoMap(ii,1,e), hist.infoMap(ii,2,e)] = infoMeasures(hist.pfMap(ii,:,e),occMap(e,:),0);
        end
    end
    for ii = 1:Npcs_asd
        z = activespk_asd{ii};

        % separate exploration in smaller intervals
        for e = 1:Nepochs
            bin_phi_e = bin_phi{ii}(e_bound(e):e_bound(e+1));
            spike_e = z(e_bound(e):e_bound(e+1));

            % Occupancy and spike rate maps
            for n = 1:Nbins
                asd.spkMap(ii,n,e) = sum(spike_e(bin_phi_e == n));
            end
            asd.normspkMap(ii,:,e) = asd.spkMap(ii,:,e)./max(asd.spkMap(ii,:,e));

            % asd estimation
            [asd.pfMap(ii,:,e),~] = runASD_1d(bin_phi_e,(spike_e)',Nbins);
            asd.normpfMap(ii,:,e) = asd.pfMap(ii,:,e)./max(asd.pfMap(ii,:,e));
            [asd.infomap(ii,1,e), asd.infomap(ii,2,e)] = ...
                infoMeasures(squeeze(asd.pfMap(ii,:,e))',ones(Nbins,1),0);
        end
    end
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


