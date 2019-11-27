% Written by Ann Go
% Function for filtering non place cells 
% Based on Indersmitten et al 2019, Front Neurosci

function [ pcIdx, pcIdx_asd ] = filter_nonPC( bin_phi, activespikes, infoMap, infoMap_asd, Nbins, prctile_thr )

Ncells = size(activespikes,1); % number of cells
info_type = 1; % 1 for info/sec, 2 for info/spike
spikeMap = zeros(Ncells,Nbins);
MI = zeros(1,1000);
pcIdx = []; pcIdx_asd = [];
% prctile_thr = 70;

if iscell(bin_phi)
    for id = 1:Ncells
        occ = histcounts(bin_phi{id},Nbins);
        z = activespikes{id};

        for k = 1:Nbins
            spikeMap(id,k) = sum(z(bin_phi{id} == k));
        end
        for j = 1:1000
            randind = randperm(length(occ));
            occMap = occ(randind);

            pcMap = spikeMap(id,:)./occMap;
            [MI(id,j),~] = infoMeasures(pcMap, occMap, 0);
        end
        if infoMap(id,info_type) > prctile(MI(id,:),prctile_thr)
            pcIdx = [pcIdx; id];
        end
        if infoMap_asd(id,info_type) > prctile(MI(id,:),prctile_thr)
            pcIdx_asd = [pcIdx_asd; id];
        end
    end
else
    occ = histcounts(bin_phi,Nbins);

    for id = 1:Ncells
        for k = 1:Nbins
            spikeMap(id,k) = sum(activespikes(bin_phi == k));
        end
        for j = 1:1000
            randind = randperm(length(occ));
            occMap = occ(randind);

            pcMap = spikeMap(id,:)./occMap;
            [MI(id,j),~] = infoMeasures(pcMap, occMap, 0);
        end
        if infoMap(id,info_type) > prctile(MI(id,:),prctile_thr)
            pcIdx = [pcIdx; id];
        end
        if infoMap_asd(id,info_type) > prctile(MI(id,:),prctile_thr)
            pcIdx_asd = [pcIdx_asd; id];
        end
    end
end

