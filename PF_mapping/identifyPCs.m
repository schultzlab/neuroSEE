% Written by Ann Go
% Function for identifying place cells 
% pcIdx_MI, pcIdx_SIsec, pcIdx_SIspk are sorted in descending order of
% infr content
% Based on Indersmitten et al 2019, Front Neurosci

function [ pcIdx_MI, pcIdx_SIsec, pcIdx_SIspk, nonpcIdx_MI, nonpcIdx_SIsec, nonpcIdx_SIspk ] ...
    = identifyPCs( spkRaster, spkPeak, bin_phi, activespikes, infoMap, Nbins, prctile_thr, randN )

if nargin<8, randN = 1000; end
if nargin<7, prctile_thr = 95; end

Ncells = size(activespikes,1); % number of cells
spikeMap = zeros(1,Nbins);
MI = zeros(1,randN);     
SIsec = zeros(1,randN); 
SIspk = zeros(1,randN); 

occMap = histcounts(bin_phi,Nbins);

% remove cells 
% 1) with peak PSD > 10 
% 2) that do not fire for half of the trials
% 3) whose information content does not exceed 99th percentile of
% shuffled distribution
include_MI = []; exclude_MI = [];
include_SIsec = []; exclude_SIsec = [];
include_SIspk = []; exclude_SIspk = [];

for id = 1:Ncells
    if spkPeak(id) < 10
        spkTr = mean(spkRaster{id},2);
        activeTr = numel(find(spkTr));
        if activeTr >= 0.5*size(spkRaster{id},1)
            z = activespikes(id,:);
            
            for k = 1:Nbins
                spikeMap(k) = sum(z(bin_phi == k));
            end
            
            for j = 1:randN
                randind = randperm(length(spikeMap));
                spikeMap = spikeMap(randind);
%                 randind = randperm(length(occMap));
%                 occMap = occMap(randind);
                pcMap = spikeMap./occMap;
                [MI(j),SIsec(j),SIspk(j)] = infoMeasures(pcMap, occMap, 0);
            end
            
            if infoMap(id,1) > prctile(MI,prctile_thr)
                include_MI = [include_MI; id];
            else
                exclude_MI = [exclude_MI; id];
            end
            
            if infoMap(id,2) > prctile(SIsec,prctile_thr)
                include_SIsec = [include_SIsec; id];
            else
                exclude_SIsec = [exclude_SIsec; id];
            end
            
            if infoMap(id,3) > prctile(SIspk,prctile_thr)
                include_SIspk = [include_SIspk; id];
            else
                exclude_SIspk = [exclude_SIspk; id];
            end
        else
            exclude_MI = [exclude_MI; id];
            exclude_SIsec = [exclude_SIsec; id];
            exclude_SIspk = [exclude_SIspk; id];
        end
    else
        exclude_MI = [exclude_MI; id];
        exclude_SIsec = [exclude_SIsec; id];
        exclude_SIspk = [exclude_SIspk; id];
    end
end

% sort cells in order of decreasing info
[~,sort_incMI] = sort(infoMap(include_MI,1),'descend');
[~,sort_excMI] = sort(infoMap(exclude_MI,1),'descend');
pcIdx_MI = include_MI(sort_incMI);
nonpcIdx_MI = exclude_MI(sort_excMI);

[~,sort_incSIsec] = sort(infoMap(include_SIsec,1),'descend');
[~,sort_excSIsec] = sort(infoMap(exclude_SIsec,1),'descend');
pcIdx_SIsec = include_SIsec(sort_incSIsec);
nonpcIdx_SIsec = exclude_SIsec(sort_excSIsec);

[~,sort_incSIspk] = sort(infoMap(include_SIspk,1),'descend');
[~,sort_excSIspk] = sort(infoMap(exclude_SIspk,1),'descend');
pcIdx_SIspk = include_SIspk(sort_incSIspk);
nonpcIdx_SIspk = exclude_SIspk(sort_excSIspk);
