% Written by Ann Go
% Function for filtering non place cells out of sorted place map
% calculation.
% Based on Indersmitten et al 2019, Front Neurosci

% function [placeMap_PC, infoMap_PC, id_PC] = filter_nonPC(activephi,spikeMap,placeMap,infoMap,Nbins)
activephi = activeData.phi; activespikes = activeData.spikes; Nbins = params.Nbins;
info_type = 1; % 1 for info/sec, 2 for info/spike
MI = zeros(1,1000);
id_PC = [];
for i = 1:size(spikeMap,1)
    for j = 1:1000
        [bin_phi,~] = discretize(activephi,Nbins);
        occ = histcounts(bin_phi,Nbins);
        randind = randperm(length(occ));
        occMap = occ(randind);
        for k = 1:Nbins
            spikeMap(i,k) = sum(activespikes(bin_phi == k));
        end

        pcMap = spikeMap(i,:)./occMap;
        [MI(i,j),~] = infoMeasures(pcMap, occMap, 0);
    end
    if infoMap(i,info_type) > prctile(MI(i,:),99)
        id_PC = [id_PC; i];
    end
end

% for i = 1:size(spikeMap,1)
%     for j = 1:1000
%         randind = randperm(length(activephi));
%         randphi = activephi(randind);
%         [bin_phi,~] = discretize(randphi,Nbins);
%         occMap = histcounts(bin_phi,Nbins);
%         for k = 1:Nbins
%             spikeMap(i,k) = sum(activespikes(bin_phi == k));
%         end
% 
%         pcMap = spikeMap(i,:)./occMap;
%         [MI(i,j),~] = infoMeasures(pcMap, occMap, 0);
%     end
%     if infoMap(i,info_type) > prctile(MI(i,:),99)
%         id_PC = [id_PC; i];
%     end
% end

placeMap_PC = placeMap(id_PC,:);
infoMap_PC = infoMap(id_PC,:);