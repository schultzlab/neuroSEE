% Written by Ann Go
% Function for identifying place cells 
% pcIdx_MI, pcIdx_SIsec, pcIdx_SIspk are sorted in descending order of
% infr content
% Based on Indersmitten et al 2019, Front Neurosci

function [ pcIdx_SIsec, pcIdx_SIspk, nonpcIdx_SIsec, nonpcIdx_SIspk ] ...
    = identifyPCs_2d( bin_pos, activespikes, infoMap, Nbins, prctile_thr, randN, method )

if nargin<7, method = 'hist'; end
if nargin<6, randN = 1000; end
if nargin<5, prctile_thr = 99; end

dt = 1/30.9;
Ncells = size(activespikes,1); % number of cells
SIsec = zeros(1,randN); 
SIspk = zeros(1,randN); 

include_SIsec = []; exclude_SIsec = [];
include_SIspk = []; exclude_SIspk = [];

for id = 1:Ncells
    z = activespikes(id,:);

    for j = 1:randN
        % shuffling method 4
        randind = randperm(length(bin_pos));
        bin_pos = bin_pos(randind);
        
        if strcmpi(method, 'hist')
            [xh, yh] = ind2sub(Nbins, bin_pos);
            spikeMap = full(sparse( xh, yh, z, Nbins(1), Nbins(2) ));
            occMap = full(sparse( xh , yh, 1, Nbins(1), Nbins(2) ));
            mode = 0; % the mask is obtained by imfill only
            envMask = getEnvEdgePrior(occMap,mode); % hist

            pcMap = spikeMap./(occMap*dt);
            pcMap(isnan(pcMap)) = 0; pcMap(isinf(pcMap)) = 0; 
            pcMap = imgaussfilt(pcMap, 2); 
            pcMap(~envMask) = 0;
            [SIsec(j),SIspk(j)] = infoMeasures(pcMap, occMap, 0);
        else
            [xa,ya] = ind2sub(Nbins, bin_pos); 
            occMap = full(sparse(xa, ya, 1, Nbins(1) , Nbins(2)));
            mode = 2; % the mask is obtained by dilation and imfill
            envMask_asd = getEnvEdgePrior(occMap,mode); % ASD
            pcMap = runASD_2d(bin_pos, z', Nbins, envMask_asd);
            [SIsec(j),SIspk(j)] = infoMeasures(pcMap', ones(Nbins), 0);
        end
    end

    if infoMap(id,1) > prctile(SIsec,prctile_thr)
        include_SIsec = [include_SIsec; id];
    else
        exclude_SIsec = [exclude_SIsec; id];
    end

    if infoMap(id,2) > prctile(SIspk,prctile_thr)
        include_SIspk = [include_SIspk; id];
    else
        exclude_SIspk = [exclude_SIspk; id];
    end
end

% sort cells in order of decreasing info
[~,sort_incSIsec] = sort(infoMap(include_SIsec,1),'descend');
[~,sort_excSIsec] = sort(infoMap(exclude_SIsec,1),'descend');
pcIdx_SIsec = include_SIsec(sort_incSIsec);
nonpcIdx_SIsec = exclude_SIsec(sort_excSIsec);

[~,sort_incSIspk] = sort(infoMap(include_SIspk,1),'descend');
[~,sort_excSIspk] = sort(infoMap(exclude_SIspk,1),'descend');
pcIdx_SIspk = include_SIspk(sort_incSIspk);
nonpcIdx_SIspk = exclude_SIspk(sort_excSIspk);
