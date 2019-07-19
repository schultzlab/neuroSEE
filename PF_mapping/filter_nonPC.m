% Written by Ann Go
% Function for filtering non place cells 
% Based on Indersmitten et al 2019, Front Neurosci

function [ pcIdx, pcIdx_asd ] = filter_nonPC( bin_phi, activespikes, infoMap, infoMap_asd, Nbins )

N = size(activespikes,1);
info_type = 1; % 1 for info/sec, 2 for info/spike
spikeMap = zeros(N,Nbins);
MI = zeros(1,1000);
pcIdx = []; pcIdx_asd = [];

occ = histcounts(bin_phi,Nbins);
        
for i = 1:N
    for k = 1:Nbins
        spikeMap(i,k) = sum(activespikes(bin_phi == k));
    end
    for j = 1:1000
        randind = randperm(length(occ));
        occMap = occ(randind);
        
        pcMap = spikeMap(i,:)./occMap;
        [MI(i,j),~] = infoMeasures(pcMap, occMap, 0);
    end
    if infoMap(i,info_type) > prctile(MI(i,:),99)
        pcIdx = [pcIdx; i];
    end
    if infoMap(i,info_type) > prctile(MI(i,:),99)
        pcIdx_asd = [pcIdx_asd; i];
    end
end

