% Written by Ann Go
% Function for identifying place cells 
% pcIdx_MI, pcIdx_SIsec, pcIdx_SIspk are sorted in descending order of
% infr content
% Based on Indersmitten et al 2019, Front Neurosci

function [ pcIdx_SIsec, pcIdx_SIspk, nonpcIdx_SIsec, nonpcIdx_SIspk ] ...
    = identifyPCs_2d( activespk, xh, yh, infoMap, pf_activet, Nbins, prctile_thr, pfactivet_thr, Nrand, mode, shuffle_method, gaussfiltSigma )

if nargin<12, gaussfiltSigma = 1.5; end
if nargin<11, shuffle_method = 2; end
if nargin<10, mode = 'hist'; end
if nargin<9, Nrand = 1000; end
if nargin<8, pfactivet_thr = 0.05; end
if nargin<7, prctile_thr = 99; end

params = neuroSEE_setparams;
fr = params.PFmap.fr;
dt = 1/fr;
Ncells = size(activespk,1); % number of cells
SIsec = zeros(1,Nrand); 
SIspk = zeros(1,Nrand); 
occMap = full(sparse( yh , xh, 1, Nbins(1), Nbins(2) ));

include_SIsec = []; exclude_SIsec = [];
include_SIspk = []; exclude_SIspk = [];

for c = 1:Ncells
    if pf_activet(c) >= pfactivet_thr
        z = activespk(c,:);
        
        if shuffle_method == 2
            % initialisation for shuffling method 2 inside subloop
            a = 927; % 30s
            b = size(activespk,2) - a; % length of recording - 30 s
            r = round((b-a)*rand(Nrand,1) + a);
            zs = zeros(size(activespk(1,:)));
        end

        for j = 1:Nrand
            if shuffle_method == 1
                % shuffling method 1 
                % shuffle spike event times
                randind = randperm(size(activespk,2));
                xh = xh(randind);
                yh = yh(randind);
                spikeMap = full(sparse( yh, xh, z, Nbins(1), Nbins(2) ));
            end
            
            if shuffle_method == 2
                % shuffling method 2
                % offset spike timeseries by a random time between 10s and length
                % of recording-10 s
                for k = r(j)+1:size(activespk,2)
                    zs(k) = z(k-r(j)); 
                end
                for k = 1:r(j)
                    zs(k) = z(size(activespk,2)-r(j)+k);
                end
                spikeMap = full(sparse( yh, xh, zs, Nbins(1), Nbins(2) ));
            end
            
            if strcmpi(mode, 'hist')
                envmode = 0; % the mask is obtained by imfill only
                envMask = getEnvEdgePrior(occMap,envmode); % hist
                pcMap = spikeMap./(occMap*dt);
                pcMap(isnan(pcMap)) = 0; pcMap(isinf(pcMap)) = 0; 
                pcMap = imgaussfilt(pcMap, gaussfiltSigma); 
                pcMap(~envMask) = 0;
                [SIsec(j),SIspk(j)] = infoMeasures(pcMap, occMap, 0);
            else
                envmode = 2; % the mask is obtained by dilation and imfill
                envMask_asd = getEnvEdgePrior(occMap,envmode); % ASD
                pcMap = runASD_2d(bin_pos, z', Nbins, envMask_asd);
                [SIsec(j),SIspk(j)] = infoMeasures(pcMap', ones(Nbins), 0);
            end
        end

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
    end
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
