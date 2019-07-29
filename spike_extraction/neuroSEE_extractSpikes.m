% Written by Ann Go

% This function extracts spikes from df_f

% INPUTS
%   df_f
%   ddf_f       : decontaminated (fissa-corrected) df_f
%   data_locn   : GCaMP data repository
%   file        : part of filename of 2P image in the format
%                   yyyymmdd_HH_MM_SS
%   params.
%       RsmoothFac   : width of smoothing window
%       g            : OASIS fluorescence impulse factor
%       lambda       : OASIS sparsity penalty
%   force       : if =1, existing R and spike data is overwritten
%
% OUTPUS
%   spikes      : spikes extracted from R
%   params

function [spikes, params, fname_mat] = neuroSEE_extractSpikes( df_f, ddf_f, data_locn, file, params, force )
   if nargin<6, force = 0; end

    mcorr_method = params.methods.mcorr_method;
    segment_method = params.methods.segment_method;

    if params.methods.dofissa
        str_fissa = 'FISSA';
    else
        str_fissa = 'noFISSA';
    end
    filedir = [data_locn,'Data/',file(1:8),'/Processed/',file,'/mcorr_',mcorr_method,'/',segment_method,'/',str_fissa,'/'];
        if ~exist(filedir,'dir'), mkdir(filedir); end
    fname_mat = [filedir file '_spikes_output.mat'];
    fname_fig = [filedir file '_spikes.fig'];

    if ~isempty(ddf_f)
        C = ddf_f;
    else
        C = df_f;
    end
    N = size(C,1); T = size(C,2);

    if force || ~exist(fname_mat,'file')
        str = sprintf( '%s: Doing spike extraction\n', file );
        cprintf(str)

        for i = 1:N
            fo = ones(1,T) * prctile(C(i,:),params.spkExtract.bl_prctile);
            C(i,:) = (C(i,:) - fo); % ./ fo;
        end
        spikes = zeros(N,T);
        for i = 1:N
            spkmin = params.spkExtract.spk_SNR*GetSn(C(i,:));
            lam = choose_lambda(exp(-1/(params.fr*params.spkExtract.decay_time)),GetSn(C(i,:)),params.spkExtract.lam_pr);

            [~,spk,~] = deconvolveCa(C(i,:),'ar2','method','thresholded','optimize_pars',true,'maxIter',20,...
                                    'window',150,'lambda',lam,'smin',spkmin);
            spikes(i,:) = spk(:);

        end

        makeplot(spikes);
        
        currstr = sprintf( '%s: Spike extraction done\n', file );
        refreshdisp(currstr,str)
    else
        spike_output = load(fname_mat);
        spikes = spike_output.spikes;
        params.spkExtract = spike_output.params;

        if ~exist(fname_fig,'file')
            makeplot(spikes);
        end

        fprintf( '%s: Spike extraction data found and loaded\n', file );
    end

    function makeplot(spikes)
        fig = figure;
        iosr.figures.multiwaveplot(1:size(spikes,2),1:size(spikes,1),spikes,'gain',5); yticks([]); xticks([]);
        title('dF/F','Fontweight','normal','Fontsize',12);
        savefig(fig,[filedir file '_spikes']);
        saveas(fig,[filedir file '_spikes'],'jpg');
        close(fig);
    end

end % function
