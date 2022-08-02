% Written by Ann Go
% Script for spike estimation

function spikes = extractSpikes( ts, options )
if nargin<2
    params = neuroSEE_setparams();
    options = params.spkExtract;
end

bl_prctile = options.bl_prctile;
spk_SNR = options.spk_SNR;
fr = options.fr;
decay_time = options.decay_time; 
lam_pr = options.lam_pr;
N = size(ts,1); T = size(ts,2);

for i = 1:N
    fo = ones(1,T) * prctile(ts(i,:), bl_prctile);
    ts(i,:) = (ts(i,:) - fo); % ./ fo;
end

spikes = zeros(N,T);
for i = 1:N
    ts_up = upsample(ts(i,:),5);
    spkmin = spk_SNR*GetSn(ts(i,:));
    lam = choose_lambda(exp(-1/(fr*decay_time)),GetSn(ts(i,:)),lam_pr);

    [~,spk,~] = deconvolveCa(ts(i,:),'ar2','method','thresholded','optimize_pars',true,'maxIter',20,...
                            'window',150,'lambda',lam,'smin',spkmin);
                        
    spikes(i,:) = downsample(spk(:),5);
end