function [PFdata, hist, asd] = calcPFdata_1d(bin_phi, activephi, activespk, activet, Nbins, histsmoothWin, fr, doasd)
if nargin<8, doasd = false; end     % flag to do asd estimation of pf
if nargin<7, fr = 30.9; end         % imaging frame rate
dt = 1/fr;

%% separate phi and spike data into trials (laps)
bp = bin_phi;
p = activephi;    
ntrial = 1;
ytick_files = 1;
Ncells = size(activespk,1);

% find the delineations for diff videos: find t = 0
idx_file = find(diff(activet) < 0);
idx_file = [0; idx_file; numel(activet)] +1; 
    
% initialise matrices
files_Ntrials = zeros(numel(idx_file)-1, 1);
tcount = (1:length(activet));

for f = 1:numel(idx_file)-1
    % bin_phi for video jj
    bp_file = bp(idx_file(f):idx_file(f+1)-1);
    
    % find the delineations for diff trials (laps)
    idx_tr = [];
    for li = 1:numel(bp_file)
        if bp_file(li) <= bp_file(1)+1 && bp_file(li) >= bp_file(1)-1
            idx_tr = [idx_tr; li];
        end
    end
    % Eliminate false delineations - if indices are less than Nbins
    % apart, they are probably due to animal staying in the same place
    % or moving back and forth (if animal moves Nbins (102 cm) in Nbins
    % time points (30.9 Hz recording), this corresponds to >600 mm/s!
    % Animal doesn't run this fast, so Nbins is a safe threshold.
    temp = idx_tr;
    for l = 1:numel(idx_tr)-1
        if (temp(l+1) - temp(l)) <= Nbins 
            idx_tr(l+1) = 0;
        end
    end
    idx_tr = idx_tr( idx_tr > 0 );
    idx1 = []; idx2 = []; idx3 = []; idx4 = [];
    temp = idx_tr;
    for l = 1:numel(idx_tr)-1
        bp_tr = bp_file(temp(l):temp(l+1)); 
        for li = 1:numel(bp_tr)
            % This trial must contain a bin close to the 360th degree position
            if bp_tr(li) <= Nbins && bp_tr(li) >= Nbins-4
                idx1 = [idx1; li];
            end
            % This trial must contain a bin close to the 270th degree position
            if bp_tr(li) <= (3*Nbins/4)+2 && bp_tr(li) >= (3*Nbins/4)-2
                idx2 = [idx2; li];
            end
            % This trial must contain a bin close to the 180th degree position
            if bp_tr(li) <= (Nbins/2)+2 && bp_tr(li) >= (Nbins/2)-2
                idx3 = [idx3; li];
            end
            % This trial must contain a bin close to the 90th degree position
            if bp_tr(li) <= (Nbins/4)+2 && bp_tr(li) >= (Nbins/4)-2
                idx4 = [idx4; li];
            end
        end
        if any([isempty(idx1) isempty(idx2) isempty(idx3) isempty(idx4)])
            idx_tr(l+1) = 0;
        end
    end
    idx_tr = idx_tr( idx_tr > 0 );
    files_Ntrials(f) = numel(idx_tr)-1;
    
    if files_Ntrials(f) > 0
        p_file = p(idx_file(f):idx_file(f+1)-1);
        tcount_file = tcount(idx_file(f):idx_file(f+1)-1);
        idx_tr(1) = 0;
        for l = 1:numel(idx_tr)-1
            bp_trials{ntrial} = bp_file(idx_tr(l)+1:idx_tr(l+1));
            phi_trials{ntrial} = p_file(idx_tr(l)+1:idx_tr(l+1));
            idx_trials{ntrial} = tcount_file(idx_tr(l)+1:idx_tr(l+1));
            for c = 1:Ncells
                s = activespk(c,:);
                s_file = s(idx_file(f):idx_file(f+1)-1);
                spk_trials{c}{ntrial} = s_file(idx_tr(l)+1:idx_tr(l+1));
            end
            ntrial = ntrial + 1;
        end

        if f == numel(idx_file)-1
            ytick_files = [ytick_files; sum(files_Ntrials(1:f))];
        else
            ytick_files = [ytick_files; sum(files_Ntrials(1:f))+1];
        end
    end
end

%% find time spent in bin for different trials
% bintime_trials:   [trials x bins] dwell time in bin per trial
% bintime :         [1 x bins]      total dwell time in bin

Ntrials = sum(files_Ntrials);
bintime_trials = zeros(Ntrials,Nbins);
for tr = 1:Ntrials
    for bin = 1:Nbins
        bintime_trials(tr,bin) = length(find(bp_trials{tr} == bin))/fr;
    end
end
bintime = sum(bintime_trials,1);

%% generate spike rastergrams and find active laps
spkRaster = cell(Ncells, 1);
normspkRaster = cell(Ncells, 1);
activetrials = zeros(Ncells, 1);
for c = 1:Ncells
    for l = 1:numel(bp_trials)
        bp_tr = bp_trials{l};
        spk_tr = spk_trials{c}{l};

        for n = 1:Nbins
            spkRaster{c}(l,n) = sum(spk_tr(bp_tr == n));
        end

        normspkRaster{c}(l,:) = spkRaster{c}(l,:)./max(spkRaster{c}(l,:));
        normspkRaster{c}(isnan(normspkRaster{c})) = 0;
    end
    meanspkRaster = mean(spkRaster{c},2);
    % number of trials (laps) for which cell was active
    activetrials(c) = numel(find(meanspkRaster))/numel(bp_trials);
end

%% calculate total event amplitude per bin, event and activity rates
L = size(activespk,2);
T = L / fr;
spk_eventrate = zeros(Ncells,1);                   % average event ("spike" rate)
spk_rate = zeros(Ncells,1);               % average event amplitude rate
bin_activet = zeros(Ncells, Nbins);         % fraction of dwell time in bin for which cell was active   
spkMap = zeros(Ncells, Nbins);              % cumulative spike map
normspkMap = zeros(Ncells, Nbins);          % normalised spike map

for c = 1:Ncells
   % get list of spikes for this cell
   z = activespk(c,:);
   spk_idx = z>0;
   spk_ampl = z(spk_idx);
   spk_eventrate(c) = length(spk_ampl)/T;
   spk_rate(c) = sum(spk_ampl)/T;
   for bin = 1:Nbins
       bin_idx = find(bin_phi == bin);
       bin_activet(c,bin) = length(find(z(bin_idx)>0))/length(bin_idx);
       spkMap(c,bin) = sum(z(bin_phi == bin));
   end
    normspkMap(c,:) = spkMap(c,:)./max(spkMap(c,:));
end


%% calculate rate maps for entire session 
hist.rateMap = zeros(Ncells, Nbins);          % place field map
hist.rateMap_sm = zeros(Ncells, Nbins);       % smoothened place field map
hist.normrateMap_sm = zeros(Ncells, Nbins);   % smoothened normalised place field map
hist.infoMap = zeros(Ncells,2);             % info values
hist.pfLoc = zeros(Ncells,1);
hist.fieldSize = zeros(Ncells,1);
hist.pfBins = zeros(Ncells,1);

if doasd
    asd.rateMap = zeros(Ncells, Nbins);           % place field map for asd
    asd.rateMap = zeros(Ncells, Nbins);       % normalised place field map for asd
    asd.infoMap = zeros(Ncells,2);              % info values
    asd.pfLoc = zeros(Ncells,1);
    asd.fieldSize = zeros(Ncells,1);
    asd.pfBins = zeros(Ncells,1);
end

occMap = histcounts(bin_phi,Nbins);
for c = 1:Ncells
    % histogram estimation
    hist.rateMap(c,:) = spkMap(c,:)./(occMap*dt);    
    hist.rateMap_sm(c,:) = circularSmooth(hist.rateMap(c,:),histsmoothWin);
    hist.normrateMap_sm(c,:) = hist.rateMap_sm(c,:)./max(hist.rateMap_sm(c,:));
    [hist.infoMap(c,1), hist.infoMap(c,2)] = infoMeasures(hist.rateMap(c,:),occMap,0);
    
    if doasd
        % ASD estimation
        [asd.rateMap(c,:),~] = runASD_1d(bin_phi,z',Nbins);
        asd.normrateMap(c,:) = asd.rateMap(c,:)./max(asd.rateMap(c,:));
        [asd.infoMap(c,1), asd.infoMap(c,2)] = infoMeasures(asd.rateMap(c,:)',ones(Nbins,1),0);
    end
end

%% find location preference and field size, active time within putative place field
[ hist.pfLoc, hist.fieldSize, hist.pfBins ] = prefLoc_fieldSize_1d( hist.rateMap_sm );

%% calculate fraction of dwell time within place field in which cell was active
for c = 1:Ncells
   % get list of spikes for this cell
   z = activespk(c,:);
   pfBins_activet = 0;
   pfBins_t = 0;
   for k = 1:length(hist.pfBins{c})
       bin = hist.pfBins{c}(k);
       bin_idx = find(bin_phi == bin);
       pfBins_activet = pfBins_activet + length(find(z(bin_idx)>0));
       pfBins_t = pfBins_t + length(bin_idx);
   end
   hist.pf_activet(c) = pfBins_activet/pfBins_t;
end
if doasd
    [ asd.pfLoc, asd.fieldSize, asd.pfBins ] = prefLoc_fieldSize_1d( asd.rateMap );
    for c = 1:Ncells
       % get list of spikes for this cell
       z = activespk(c,:);
       pfBins_activet = 0;
       pfBins_t = 0;
       for k = 1:length(hist.pfBins{c})
           bin = asd.pfBins{c}(k);
           bin_idx = find(bin_phi == bin);
           pfBins_activet = pfBins_activet + length(find(z(bin_idx)>0));
           pfBins_t = pfBins_t + length(bin_idx);
       end
       asd.pf_activet(c) = pfBins_activet/pfBins_t;
    end
end

%% output
PFdata.files_Ntrials = files_Ntrials;
PFdata.phi_trials = phi_trials;
PFdata.spk_trials = spk_trials;
PFdata.idx_trials = idx_trials;
PFdata.bintime_trials = bintime_trials;
PFdata.bintime = bintime;
PFdata.spkRaster = spkRaster;
PFdata.normspkRaster = normspkRaster;
PFdata.ytick_files = ytick_files;
PFdata.activetrials = activetrials;
PFdata.spk_eventrate = spk_eventrate;
PFdata.spk_rate = spk_rate;
PFdata.bin_activet = bin_activet;
PFdata.occMap = occMap;
PFdata.spkMap = spkMap;
PFdata.normspkMap = normspkMap;

if ~doasd, asd = []; end