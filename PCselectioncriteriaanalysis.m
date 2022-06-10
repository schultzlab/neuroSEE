% criteria parameters
Nbins = 50;
Nrand = 1000;
Vthr = 20;
prctile_thr = 99; %1
rate_thr = 3; %2
activetrials_thr = 0.4; %3
pfactivet_thr = 0.03; %4

% regenerate some variables
phi = downTrackdata.phi;
speed = downTrackdata.speed;
activeind = find(speed > Vthr);
activephi   = phi(activeind);
activespk   = spikes(:,activeind);
[bin_phi,~] = discretize(activephi,Nbins);

infoMap = hist.infoMap;
pf_activet = hist.pf_activet; 
activetrials = PFdata.activetrials;

dt = 1/30.9;
Nbins = max(bin_phi);
Ncells = size(activespk,1); % number of cells
spikeMap = zeros(1,Nbins);
SIsec = zeros(1,Nrand); 
SIspk = zeros(1,Nrand); 
occMap = histcounts(bin_phi,Nbins);

include_SIsec = []; 
include_SIspk = []; 

% remove cells 
% 1) with info < info of 99th percentile of shuffled distribution
exclude_SIsec1 = [];
exclude_SIspk1 = [];
for c = 1:Ncells
    z = activespk(c,:);

    % initialisation for shuffling method 
    a = 309; % 10s
    b = numel(bin_phi) - a; % length of recording - 10 s
    r = round((b-a)*rand(Nrand,1) + a);
    zs = zeros(size(activespk(1,:)));

    for j = 1:Nrand
        % offset spike timeseries by a random time between 10s and length
        % of recording-10 s
        for k = r(j)+1:numel(bin_phi)
            zs(k) = z(k-r(j)); 
        end
        for k = 1:r(j)
            zs(k) = z(numel(bin_phi)-r(j)+k);
        end
        for k = 1:Nbins
            spikeMap(k) = sum(zs(bin_phi == k));
        end

        pcMap = spikeMap./(occMap*dt);
        pcMap(isnan(pcMap)) = 0; pcMap(isinf(pcMap)) = 0; 
        pcMap_sm = circularSmooth(pcMap,5);
        [SIsec(j),SIspk(j)] = infoMeasures(pcMap_sm, occMap, 0);
    end

    if infoMap(c,1) < prctile(SIsec,prctile_thr)
        exclude_SIsec1 = [exclude_SIsec1; c];
    end

    if infoMap(c,2) < prctile(SIspk,prctile_thr)
        exclude_SIspk1 = [exclude_SIspk1; c];
    end
end
% plotSpikeRasterTrials( PFdata.normspkRaster(exclude_SIspk1), PFdata.ytick_files, 'NonPC', false, [], false );

% remove cells 
% 2) for which the mean in-field event rate < 3x mean out-of-field event
% rate
exclude_SIsec2 = [];
exclude_SIspk2 = [];
for bin = 1:Nbins
   t_bin(bin) = length(find(bin_phi == bin));
end
for c = 1:Ncells
    z = activespk(c,:);
    for bin = 1:Nbins
       spikes_bin(c,bin) = sum(z(bin_phi == bin));
       rate_bin(c,bin) = spikes_bin(c,bin)/t_bin(bin);
    end
    infield = sum(rate_bin(c,hist.pfBins{c}))/length(hist.pfBins{c});
    out = setdiff(1:Nbins,hist.pfBins{c});
    outfield = sum(rate_bin(c,out))/length(out);
    
    if infield < rate_thr*outfield
        exclude_SIspk2 = [exclude_SIspk2; c];
    end
end
% plotSpikeRasterTrials( PFdata.normspkRaster(exclude_SIspk2), PFdata.ytick_files, 'NonPC', false, [], false );

% remove cells
% 3) that are not active for more than activetrials_thr of total number of
% trials (laps)
exclude_SIsec3 = [];
exclude_SIspk3 = [];
for c = 1:Ncells
    if activetrials(c) < activetrials_thr
        exclude_SIsec3 = unique([exclude_SIsec3; c]);
        exclude_SIspk3 = unique([exclude_SIspk3; c]);
    end 
end
% plotSpikeRasterTrials( PFdata.normspkRaster(exclude_SIspk3), PFdata.ytick_files, 'NonPC', false, [], false );

exclude_SIspk = unique([exclude_SIspk1; exclude_SIspk2; exclude_SIspk3]);
include_SIspk = setdiff(1:Ncells,exclude_SIspk);
[~,sort_exclude_SIspk] = sort(infoMap(exclude_SIspk,1),'descend');
[~,sort_include_SIspk] = sort(infoMap(include_SIspk,1),'descend');
plotSpikeRasterTrials( PFdata.normspkRaster(exclude_SIspk), PFdata.ytick_files, 'NonPC', false, [], false );
plotSpikeRasterTrials( PFdata.normspkRaster(include_SIspk), PFdata.ytick_files, 'PC', false, [], false );
title('A');

% remove cells 
% 4) that are not active for more than pfactivet_thr of total dwell time inside place field
exclude_SIsec4 = [];
exclude_SIspk4 = [];
for c = 1:Ncells
    if pf_activet(c) < pfactivet_thr
        exclude_SIsec4 = unique([exclude_SIsec4; c]);
        exclude_SIspk4 = unique([exclude_SIspk4; c]);
    end
end
exclude_SIspk = unique([exclude_SIspk1; exclude_SIspk2; exclude_SIspk3; exclude_SIspk4]);
include_SIspk = setdiff(1:Ncells,exclude_SIspk);
[~,sort_exclude_SIspk] = sort(infoMap(exclude_SIspk,1),'descend');
[~,sort_include_SIspk] = sort(infoMap(include_SIspk,1),'descend');
plotSpikeRasterTrials( PFdata.normspkRaster(include_SIspk), PFdata.ytick_files, 'PC', false, [], false );
title('B');

exclude_SIspk = unique([exclude_SIspk1; exclude_SIspk2; exclude_SIspk4]);
include_SIspk = setdiff(1:Ncells,exclude_SIspk);
[~,sort_exclude_SIspk] = sort(infoMap(exclude_SIspk,1),'descend');
[~,sort_include_SIspk] = sort(infoMap(include_SIspk,1),'descend');
plotSpikeRasterTrials( PFdata.normspkRaster(include_SIspk), PFdata.ytick_files, 'PC', false, [], false );
title('C');

pfactivet_thr = 0.04; %4
exclude_SIsec4 = [];
exclude_SIspk4 = [];
for c = 1:Ncells
    if pf_activet(c) < pfactivet_thr
        exclude_SIsec4 = unique([exclude_SIsec4; c]);
        exclude_SIspk4 = unique([exclude_SIspk4; c]);
    end 
end
exclude_SIspk = unique([exclude_SIspk1; exclude_SIspk2; exclude_SIspk4]);
include_SIspk = setdiff(1:Ncells,exclude_SIspk);
[~,sort_exclude_SIspk] = sort(infoMap(exclude_SIspk,1),'descend');
[~,sort_include_SIspk] = sort(infoMap(include_SIspk,1),'descend');
plotSpikeRasterTrials( PFdata.normspkRaster(include_SIspk), PFdata.ytick_files, 'PC', false, [], false );
title('D');

rate_thr = 4;
exclude_SIsec2 = [];
exclude_SIspk2 = [];
for bin = 1:Nbins
   t_bin(bin) = length(find(bin_phi == bin));
end
for c = 1:Ncells
    z = activespk(c,:);
    for bin = 1:Nbins
       spikes_bin(c,bin) = sum(z(bin_phi == bin));
       rate_bin(c,bin) = spikes_bin(c,bin)/t_bin(bin);
    end
    infield = sum(rate_bin(c,hist.pfBins{c}))/length(hist.pfBins{c});
    out = setdiff(1:Nbins,hist.pfBins{c});
    outfield = sum(rate_bin(c,out))/length(out);
    
    if infield < rate_thr*outfield
        exclude_SIspk2 = [exclude_SIspk2; c];
    end
end
exclude_SIspk = unique([exclude_SIspk1; exclude_SIspk2; exclude_SIspk3]);
include_SIspk = setdiff(1:Ncells,exclude_SIspk);
[~,sort_exclude_SIspk] = sort(infoMap(exclude_SIspk,1),'descend');
[~,sort_include_SIspk] = sort(infoMap(include_SIspk,1),'descend');
plotSpikeRasterTrials( PFdata.normspkRaster(include_SIspk), PFdata.ytick_files, 'PC', false, [], false );
title('E');
