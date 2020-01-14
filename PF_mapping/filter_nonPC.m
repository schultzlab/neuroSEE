% Written by Ann Go
% Function for filtering non place cells 
% Based on Indersmitten et al 2019, Front Neurosci

function [ pcIdx, pcIdx_asd ] = filter_nonPC( bin_phi, activespikes, infoMap, infoMap_asd, Nbins, prctile_thr, randN, info, info_type )

if nargin<9, info_type = 1; end % 1 for info/sec, 2 for info/spike
if nargin<8, info = 'MI'; end
if nargin<7, randN = 1000; end
if nargin<6, prctile_thr = 99; end

Ncells = size(activespikes,1); % number of cells
spikeMap = zeros(Ncells,Nbins);
I = zeros(1,randN);
pcIdx = []; pcIdx_asd = [];

if iscell(bin_phi)
    for id = 1:Ncells
        occ = histcounts(bin_phi{id},Nbins);
        z = activespikes{id};

        for k = 1:Nbins
            spikeMap(id,k) = sum(z(bin_phi{id} == k));
        end
        for j = 1:randN
            randind = randperm(length(occ));
            occMap = occ(randind);

            pcMap = spikeMap(id,:)./occMap;
            [I(id,j),~] = infoMeasures(pcMap, occMap, 0, info);
        end
        if infoMap(id,info_type) > prctile(I(id,:),prctile_thr)
            pcIdx = [pcIdx; id];
        end
        if infoMap_asd(id,info_type) > prctile(I(id,:),prctile_thr)
            pcIdx_asd = [pcIdx_asd; id];
        end
    end
else
    occ = histcounts(bin_phi,Nbins);

    for id = 1:Ncells
        for k = 1:Nbins
            spikeMap(id,k) = sum(activespikes(bin_phi == k));
        end
        for j = 1:randN
            randind = randperm(length(occ));
            spikeMap(id,k) = occ(randind);

            pcMap = spikeMap(id,:)./occMap;
            [I(id,j),~] = infoMeasures(pcMap, occMap, 0, info);
        end
        if infoMap(id,info_type) > prctile(I(id,:),prctile_thr)
            pcIdx = [pcIdx; id];
        end
        if infoMap_asd(id,info_type) > prctile(I(id,:),prctile_thr)
            pcIdx_asd = [pcIdx_asd; id];
        end
    end
end

