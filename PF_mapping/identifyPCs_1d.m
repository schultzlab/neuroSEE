% Written by Ann Go
% Function for identifying place cells 
% pcIdx_MI, pcIdx_SIsec, pcIdx_SIspk are sorted in descending order of
% info content
% Based on Indersmitten et al 2019, Front Neurosci and Gonzales et al 2019,
% Science

function [ pcIdx_SIsec, pcIdx_SIspk, nonpcIdx_SIsec, nonpcIdx_SIspk ] = identifyPCs_1d( ...
    bin_phi, activespk, infoMap, pf_activet, activetrials, pfBins, prctile_thr, pfactivet_thr, activetrials_thr, fieldrate_thr, Nrand, mode, shuffle_method )

if nargin<13, shuffle_method = 2; end
if nargin<12, mode = 'hist'; end
if nargin<11, Nrand = 1000; end
if nargin<10, fieldrate_thr = 3; end
if nargin<9, activetrials_thr = 0.35; end
if nargin<8, pfactivet_thr = 0.04; end
if nargin<7, prctile_thr = 99; end

dt = 1/30.9;
Nbins = max(bin_phi);
Ncells = size(activespk,1); % number of cells
spikeMap = zeros(1,Nbins);
SIsec = zeros(1,Nrand); 
SIspk = zeros(1,Nrand); 
occMap = histcounts(bin_phi,Nbins);
t_bin = zeros(1,Nbins);
spikes_bin = zeros(Ncells,Nbins);
rate_bin = zeros(Ncells,Nbins);
                    
for bin = 1:Nbins
   t_bin(bin) = length(find(bin_phi == bin));
end

% remove cells 
% 1) with info < info of 99th percentile of shuffled distribution
% 2) that are not active for more than pfactivet_thr of total dwell time inside place field
% 3) that are not active for more than activetrials_thr of total number of
% trials (laps)
% 4) with mean in-field event rate < 3x mean out-of-field event rate

include_SIsec = []; exclude_SIsec = [];
include_SIspk = []; exclude_SIspk = [];

for c = 1:Ncells
    if pf_activet(c) >= pfactivet_thr % criteria 2
        if activetrials(c) >= activetrials_thr % criteria 3
            z = activespk(c,:);
            for bin = 1:Nbins
               spikes_bin(c,bin) = sum(z(bin_phi == bin));
               rate_bin(c,bin) = spikes_bin(c,bin)/t_bin(bin);
            end
            infield = sum(rate_bin(c,pfBins{c}))/length(pfBins{c});
            out = setdiff(1:Nbins,pfBins{c});
            outfield = sum(rate_bin(c,out))/length(out);

            if infield >= fieldrate_thr*outfield % criteria 4
                if shuffle_method == 2
                    % initialisation for shuffling method 2 inside subloop
                    a = 309; % 10s
                    b = numel(bin_phi) - a; % length of recording - 10 s
                    r = round((b-a)*rand(Nrand,1) + a);
                    zs = zeros(size(activespk(1,:)));
                end

                for j = 1:Nrand
                    if shuffle_method == 1
                        % shuffling method 1 
                        % shuffle spike event times
                        randind = randperm(length(z));
                        z = z(randind);
                        for k = 1:Nbins
                            spikeMap(k) = sum(z(bin_phi == k));
                        end
                    end

                    if shuffle_method == 2
                        % shuffling method 2
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
                    end

                    if strcmpi(mode, 'hist')
                        pcMap = spikeMap./(occMap*dt);
                        pcMap(isnan(pcMap)) = 0; pcMap(isinf(pcMap)) = 0; 
                        pcMap_sm = circularSmooth(pcMap,5);
                        [SIsec(j),SIspk(j)] = infoMeasures(pcMap_sm, occMap, 0);
                    else
                        pcMap = runASD_1d(bin_phi,z',Nbins);
                        [SIsec(j),SIspk(j)] = infoMeasures(pcMap', ones(Nbins,1), 0);
                    end
                end

                % criteria 1
                if infoMap(c,1) > prctile(SIsec,prctile_thr)
                    include_SIsec = [include_SIsec; c];
                else
                    exclude_SIsec = [exclude_SIsec; c];
                end

                if infoMap(c,2) > prctile(SIspk,prctile_thr)
                    include_SIspk = [include_SIspk; c];
                else
                    exclude_SIspk = [exclude_SIspk; c];
                end
            else
                exclude_SIsec = [exclude_SIsec; c];
                exclude_SIspk = [exclude_SIspk; c];
            end % if infield >= fieldrate_thr*outfield
        else
            exclude_SIsec = [exclude_SIsec; c];
            exclude_SIspk = [exclude_SIspk; c];
        end % if activetrials(c) < activetrials_thr
    else
        exclude_SIsec = [exclude_SIsec; c];
        exclude_SIspk = [exclude_SIspk; c];
    end % if pf_activet(c) >= pfactivet_thr
end

% sort cells in order of decreasing info
[~,sort_incSIsec] = sort(infoMap(include_SIsec,1),'descend');
[~,sort_excSIsec] = sort(infoMap(exclude_SIsec,1),'descend');
pcIdx_SIsec = include_SIsec(sort_incSIsec);
nonpcIdx_SIsec = exclude_SIsec(sort_excSIsec);

[~,sort_incSIspk] = sort(infoMap(include_SIspk,2),'descend');
[~,sort_excSIspk] = sort(infoMap(exclude_SIspk,2),'descend');
pcIdx_SIspk = include_SIspk(sort_incSIspk);
nonpcIdx_SIspk = exclude_SIspk(sort_excSIspk);
