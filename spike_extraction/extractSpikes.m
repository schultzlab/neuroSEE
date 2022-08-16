% Written by Ann Go
% Script for spike estimation

function [spikes, df_f_deconv, bl, neuron_sn, g] = extractSpikes( ts, options )
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

spikes = zeros(N,T);        % deconvolved neural activity
df_f_deconv = zeros(N,T);   % deconvolved DF/F traces
bl = zeros(N,1);            % baseline for each trace (should be close to zero since traces are DF/F)
neuron_sn = zeros(N,1);     % noise level at each trace
g = cell(N,1);              % discrete time constants for each trace
for i = 1:N
    spkmin = spk_SNR*GetSn(ts(i,:));
    lam = choose_lambda(exp(-1/(fr*decay_time)),GetSn(ts(i,:)),lam_pr);

    [cc,spk,opts_oasis] = deconvolveCa(ts(i,:),'ar2','method','thresholded','optimize_pars',true,'maxIter',20,...
                            'window',150,'lambda',lam,'smin',spkmin);
                        
    bl(i) = opts_oasis.b;
    df_f_deconv(i,:) = cc(:)' + bl(i);
    spikes(i,:) = spk;
    neuron_sn(i) = opts_oasis.sn;
    g{i} = opts_oasis.pars(:)';
end