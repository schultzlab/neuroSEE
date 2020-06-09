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
    spkmin = spk_SNR*GetSn(ts(i,:));
    lam = choose_lambda(exp(-1/(fr*decay_time)),GetSn(ts(i,:)),lam_pr);

    [~,spk,~] = deconvolveCa(ts(i,:),'ar2','method','thresholded','optimize_pars',true,'maxIter',20,...
                            'window',150,'lambda',lam,'smin',spkmin);
    spikes(i,:) = spk(:);
end

%% extract DF/F and deconvolve DF/F traces

[F_dff,F0] = detrend_df_f(A_keep,[b,ones(prod(FOV),1)],C_full,[f_full;-double(F_dark)*ones(1,T)],R_full,options);

C_dec = zeros(N,T);         % deconvolved DF/F traces
S_dec = zeros(N,T);         % deconvolved neural activity
bl = zeros(N,1);            % baseline for each trace (should be close to zero since traces are DF/F)
neuron_sn = zeros(N,1);     % noise level at each trace
g = cell(N,1);              % discrete time constants for each trace
if p == 1; model_ar = 'ar1'; elseif p == 2; model_ar = 'ar2'; else; error('This order of dynamics is not supported'); end

for i = 1:N
    spkmin = options.spk_SNR*GetSn(F_dff(i,:));
    lam = choose_lambda(exp(-1/(options.fr*options.decay_time)),GetSn(F_dff(i,:)),options.lam_pr);
    [cc,spk,opts_oasis] = deconvolveCa(F_dff(i,:),model_ar,'method','thresholded','optimize_pars',true,'maxIter',20,...
                                'window',150,'lambda',lam,'smin',spkmin);
    bl(i) = opts_oasis.b;
    C_dec(i,:) = cc(:)' + bl(i);
    S_dec(i,:) = spk(:);
    neuron_sn(i) = opts_oasis.sn;
    g{i} = opts_oasis.pars(:)';
    disp(['Performing deconvolution. Trace ',num2str(i),' out of ',num2str(N),' finished processing.'])
end