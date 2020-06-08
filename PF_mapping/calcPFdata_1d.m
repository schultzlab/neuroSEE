function [PFdata, hist, asd] = calcPFdata_1d(bin_phi, activephi, activespk, activet, Nbins, fr)
if nargin<6, fr = 30.9; end

%% separate phi and spike data into trials (laps)
bp = bin_phi;
p = activephi;    
ntrial = 1;
ytick_files = 1;

% find the delineations for diff videos: find t = 0
idx_file = find(diff(activet) < 0);
idx_file = [0; idx_file; numel(activet)] +1; 
    
for f = 1:numel(idx_file)-1
    % bin_phi for video jj
    bp_file = bp(idx_file(f):idx_file(f+1)-1);
    
    % find the delineations for diff trials (laps)
    idx_tr = find( bp_file == bp_file(1) );
    % Eliminate false delineations - if indices are less than 2*Nbins
    % apart, they are probably due to animal moving back and forth or any
    % artefacts in tracking (if animal moves 2*Nbins (204 cm) in 2*Nbins
    % time points (30.9 Hz recording), this corresponds to >300 mm/s!
    % Animal doesn't run this fast, so 2*Nbins is a safe threshold.
    for l = numel(idx_tr):-1:2
        if (idx_tr(l) - idx_tr(l-1)) <= 2*Nbins 
            idx_tr(l) = 0;
        end
    end
    idx_tr = idx_tr( idx_tr > 0 );
    Ntrials(f) = numel(idx_tr)-1;
    
    if Ntrials(f) > 0
        p_file = p(idx_file(f):idx_file(f+1)-1);
        for l = 1:numel(idx_tr)-1
            bp_trials{ntrial} = bp_file(idx_tr(l)+1:idx_tr(l+1));
            phi_trials{ntrial} = p_file(idx_tr(l)+1:idx_tr(l+1));
            for c = 1:Ncells
                s = activespk(c,:);
                s_file = s(idx_file(f):idx_file(f+1)-1);
                spk_trials{c}{ntrial} = s_file(idx_tr(l)+1:idx_tr(l+1));
            end
            ntrial = ntrial + 1;
        end

        if f == numel(idx_file)-1
            ytick_files = [ytick_files; sum(Ntrials(1:f))];
        else
            ytick_files = [ytick_files; sum(Ntrials(1:f))+1];
        end
    end
end

%% find time spent in bin for different trials
bintime_trials = zeros(Ntrials,Nbins);
for tr = 1:Ntrials
    for bin = 1:Nbins
        bintime_trials(tr,bin) = length(find(bp_trials{tr} == bin))/fr;
    end
end
bintime = sum(bintime_trials,1);

%% generate spike rastergrams and find active laps
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
    activetrials = numel(find(meanspkRaster))/numel(bp_trials);
end

%% calculate event counts per bin, event and activity rates
L = size(activespk,2);
T = L / fr;
spk_rate = zeros(Ncells);
spk_amplrate = zeros(Ncells);
bin_activity = zeros(Ncells, Nbins);
spkMap = zeros(Ncells, Nbins);              % spike map
normspkMap = zeros(Ncells, Nbins);          % normalised spike map

for c = 1:Ncells
   % get list of spikes for this cell
   z = activespk(c,:);
   spk_idx = z(c,:)>0;
   spk_ampl = z(c,spk_idx);
   spk_rate(c) = length(spk_ampl)/T;
   spk_amplrate(c) = sum(spk_ampl)/T;
   for bin = 1:Nbins
       bin_idx = find(bin_phi == bin);
       bin_activity(c,bin) = length(find(activespk(c,bin_idx)>0))/length(bin_idx);
       spkMap(c,bin) = sum(z(bin_phi == bin));
   end
    normspkMap(c,:) = spkMap(c,:)./max(spkMap(c,:));
end


%% calculate rate maps for entire session 
hist.rMap = zeros(Ncells, Nbins);          % place field map
hist.rMap_sm = zeros(Ncells, Nbins);       % smoothened place field map
hist.normrMap_sm = zeros(Ncells, Nbins);   % smoothened normalised place field map
hist.infoMap = zeros(Ncells,2);             % info values
hist.pfLoc = zeros(Ncells,1);
hist.fieldSize = zeros(Ncells,1);
hist.pfBins = zeros(Ncells,1);

if doasd
    asd.rMap = zeros(Ncells, Nbins);           % place field map for asd
    asd.rMap = zeros(Ncells, Nbins);       % normalised place field map for asd
    asd.infoMap = zeros(Ncells,2);              % info values
    asd.pfLoc = zeros(Ncells,1);
    asd.fieldSize = zeros(Ncells,1);
    asd.pfBins = zeros(Ncells,1);
end

occMap = histcounts(bin_phi,Nbins);
for c = 1:Ncells
    % histogram estimation
    hist.rMap(c,:) = spkMap(c,:)./(occMap*dt);    
    hist.rMap_sm(c,:) = circularSmooth(hist.rMap(c,:),histsmoothWin);
    hist.normrMap_sm(c,:) = hist.rMap_sm(c,:)./max(hist.rMap_sm(c,:));
    [hist.infoMap(c,1), hist.infoMap(c,2)] = infoMeasures(hist.rMap(c,:),occMap,0);
    
    if doasd
        % ASD estimation
        [asd.rMap(c,:),~] = runASD_1d(bin_phi,z',Nbins);
        asd.normrMap(c,:) = asd.rMap(c,:)./max(asd.rMap(c,:));
        [asd.infoMap(c,1), asd.infoMap(c,2)] = infoMeasures(asd.rMap(c,:)',ones(Nbins,1),0);
    end
end

%% find location preference and field size
[ hist.pfLoc, hist.fieldSize, hist.pfBins ] = pfLoc_fieldSize_1d( hist.rMap_sm );

if doasd
    [ asd.pfLoc, asd.fieldSize, asd.pfBins ] = pfLoc_fieldSize_1d( asd.rMap );
end

%% output
PFdata.phi_trials = phi_trials;
PFdata.spk_trials = spk_trials;
PFdata.bintime_trials = bintime_trials;
PFdata.bintime = bintime;
PFdata.spkRaster = spkRaster;
PFdata.normspkRaster = normspkRaster;
PFdata.ytick_files = ytick_files;
PFdata.activetrials = activetrials;
PFdata.spk_rate = spk_rate;
PFdata.spk_amplrate = spk_amplrate;
PFdata.bin_activity = bin_activity;
PFdata.occMap = occMap;
PFdata.spkMap = spkMap;
PFdata.normspkMap = normspkMap;