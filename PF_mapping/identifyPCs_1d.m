% Written by Ann Go
% Function for identifying place cells 
% pcIdx_MI, pcIdx_SIsec, pcIdx_SIspk are sorted in descending order of
% infr content
% Based on Indersmitten et al 2019, Front Neurosci

function [ pcIdx_SIsec, pcIdx_SIspk, nonpcIdx_SIsec, nonpcIdx_SIspk ] ...
    = identifyPCs_1d( spkRaster, spkPeak, bin_phi, activespikes, infoMap, Nbins, prctile_thr, randN, mode )

if nargin<9, mode = 'hist'; end
if nargin<8, randN = 1000; end
if nargin<7, prctile_thr = 99; end

dt = 1/30.9;
Ncells = size(activespikes,1); % number of cells
spikeMap = zeros(1,Nbins);
SIsec = zeros(1,randN); 
SIspk = zeros(1,randN); 

% remove cells 
% 1) with peak PSD > 10 
% 2) that do not fire for half of the trials
% 3) whose information content does not exceed 99th percentile of
% shuffled distribution
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
                % shuffling method 4
                randind = randperm(length(bin_phi));
                bin_phi = bin_phi(randind);
                
                if strcmpi(mode, 'hist')
                    occMap = histcounts(bin_phi,Nbins);
                    pcMap = spikeMap./(occMap*dt);
                    pcMap(isnan(pcMap)) = 0; pcMap(isinf(pcMap)) = 0; 
                    pcMap_sm = circularSmooth(pcMap,5);
                    [SIsec(j),SIspk(j)] = infoMeasures(pcMap_sm, occMap, 0);
                else
                    pcMap = runASD_1d(bin_phi,z',Nbins);
                    [SIsec(j),SIspk(j)] = infoMeasures(pcMap', ones(Nbins,1), 0);
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
        else
            exclude_SIsec = [exclude_SIsec; id];
            exclude_SIspk = [exclude_SIspk; id];
        end
    else
        exclude_SIsec = [exclude_SIsec; id];
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

% Other shuffling methods
% shuffling method 1
% a = 618;
% b = numel(bin_phi) - a;
% r = round((b-a)*rand(randN,1) + a);
% zs = zeros(size(z));
% z = activespikes(id,:);
% for k = r(j)+1:numel(bin_phi)
%     zs(k) = z(k-r(j)); 
% end
% for k = 1:r(j)
%     zs(k) = z(numel(bin_phi)-r(j)+k);
% end
% for k = 1:Nbins
%     spikeMap(k) = sum(zs(bin_phi == k));
% end

% shuffling method 2    
% z = activespikes(id,:);
% randind = randperm(length(z));
% z = z(randind);
% for k = 1:Nbins
%     spikeMap(k) = sum(z(bin_phi == k));
% end

% shuffling method 3
% randind = randperm(length(spikeMap));
% spikeMap = spikeMap(randind);


            

