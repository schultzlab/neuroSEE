% Written by Ann Go

% This function extracts spikes from df_f or ddf_f if it exists

% INPUTS
%   df_f        : deltaf/f time series from green channel (ROI segmentation output)
%   ddf_f       : decontaminated (fissa-corrected) df_f
%   params.
%       RsmoothFac   : width of smoothing window
%       g            : OASIS fluorescence impulse factor
%       lambda       : OASIS sparsity penalty
%   fileorlist  : file - part of file name of image stacks in the format
%               yyyymmdd_hh_mm_ss
%               * Can also be a list - name of text file containing filenames of files 
%               to be concatenated for ROI segmentation. Typically in the format 
%               'list_m##_expname.txt'. When a list is given, a reffile
%               should be provided.%   conc_runs  : flag if rois were segmented from concatenated files from
%               different runs e.g. fam1fam2fam1-fam1 but rois are
%               for fam1fam2fam1. DO NOT flag for fam1fam2fam1 files since in
%               this case it is understood that the rois are from the
%               concatenated runs.
%   reffile     : (optional) file to be used as registration template. Required if 
%               fileorlist above is a list.This file is usually part of 'list'
%               but does not have to be. 
%   force       : (optional, default: false) if true, existing R and spike data is overwritten
%
% OUTPUS
%   spikes      : spikes extracted from R
%   params

function [spikes, df_f_deconv, bl, neuron_sn, g, params] = neuroSEE_extractSpikes( df_f, ddf_f, data_locn, params, fileorlist, reffile, grp_sdir, force )
    if nargin<8, force = 0; end
    if nargin<7, list = []; end
    if nargin<6, reffile = []; end

    mcorr_method = params.methods.mcorr_method;
    segment_method = params.methods.segment_method;
    if params.methods.dofissa
        str_fissa = 'FISSA';
    else
        str_fissa = 'noFISSA';
    end
    
    % determine whether fileorlist is a file or list
    if ~strncmp(fileorlist,'list',4)
        file = fileorlist; list = [];
    else
        list = fileorlist; file = [];
    end

    if isempty(list) % file
        savedir = [data_locn 'Data/' file(1:8) '/Processed/' file '/mcorr_' mcorr_method '/'...
                   segment_method '/' str_fissa '/'];
        if ~exist([savedir file '_spikes_output.mat'],'file')
            fname_mat = [savedir file '_spikes.mat'];
        else
            fname_mat = [savedir file '_spikes_output.mat'];
        end
        fname_fig = [savedir file '_spikes.fig'];
    else
        [ mouseid, expname ] = find_mouseIDexpname(list);
        savedir = [grp_sdir '/' str_fissa '/'];
        fname_mat = [savedir mouseid '_' expname '_ref' reffile '_spikes.mat'];
        fname_fig = [savedir mouseid '_' expname '_ref' reffile '_spikes.fig'];
    end    

    if ~isempty(ddf_f)
        C = ddf_f;
    else
        C = df_f;
    end

    prevstr = [];
    if force || ~exist(fname_mat,'file')
        if isempty(list)
            str = sprintf( '%s: Doing spike estimation\n', file );
        else
            str = sprintf( '%s: Doing spike estimation\n', [mouseid '_' expname] );
        end
        refreshdisp(str, prevstr);
        prevstr = str;

        [spikes, df_f_deconv, bl, neuron_sn, g] = extractSpikes( C, params.spkExtract );
        
        % Save output
        spike_output.spikes = spikes;
        spike_output.df_f_deconv = df_f_deconv;
        spike_output.bl = bl;
        spike_output.neuron_sn = neuron_sn;
        spike_output.g = g;
        spike_output.params = params.spkExtract;
        
        if isempty(list)
            str = sprintf( '%s: Saving spike data\n', file );
        else
            str = sprintf( '%s: Saving spike data\n', [mouseid '_' expname] );
        end
        refreshdisp(str, prevstr);
        prevstr = str;
        
        if ~exist(savedir,'dir'), mkdir(savedir); fileattrib(savedir,'+w','g','s'); end
        save(fname_mat,'-struct','spike_output');

        % Make and save plot
        plotSpikes(spikes, fname_fig(1:end-4));
        
        if isempty(list)
            str = sprintf( '%s: Spike estimation done\n', file );
        else
            str = sprintf( '%s: Spike estimation done\n', [mouseid '_' expname] );
        end
        refreshdisp(str, prevstr);
    else
        if isempty(list)
            str = sprintf('%s: Loading spike data\n', file);
        else
            str = sprintf('%s: Loading spike data\n', [mouseid '_' expname]);
        end
        refreshdisp(str, prevstr);
        prevstr = str;
        spike_output = load(fname_mat);
        spikes = spike_output.spikes;
        df_f_deconv = spike_output.df_f_deconv;
        bl = spike_output.bl;
        neuron_sn = spike_output.neuron_sn;
        g = spike_output.g;
        params.spkExtract = spike_output.params;
        
        if isempty(list)
            str = sprintf('%s: Spike data loaded\n', file);
        else
            str = sprintf('%s: Spike data loaded\n', [mouseid '_' expname]);
        end
        refreshdisp(str, prevstr);
        prevstr = str;
        
        if ~exist(fname_fig,'file')
            plotSpikes(spikes, fname_fig(1:end-4));
        end

        if isempty(list)
            str = sprintf( '%s: Spike estimation data found and loaded\n', file );
        else
            str = sprintf( '%s: Spike estimation data found and loaded\n', [mouseid '_' expname] );
        end
        refreshdisp(str, prevstr);
    end

end % function
