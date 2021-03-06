% Written by Ann Go

% This function extracts spikes from df_f or ddf_f if it exists

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

function [spikes, params, fname_mat] = neuroSEE_extractSpikes( df_f, ddf_f, data_locn, file, params, force, list, reffile, fsave )
    if nargin<9, fsave = true; end
    if nargin<8, reffile = []; end
    if nargin<7, list = []; end
    if nargin<6, force = 0; end

    mcorr_method = params.methods.mcorr_method;
    segment_method = params.methods.segment_method;
    bl_prctile = params.spkExtract.bl_prctile;
    if params.methods.dofissa
        str_fissa = 'FISSA';
    else
        str_fissa = 'noFISSA';
    end
    
    if isempty(list)
        filedir = [data_locn 'Data/' file(1:8) '/Processed/' file '/mcorr_' mcorr_method '/'...
                   segment_method '/' str_fissa '/bl_prctile' num2str(bl_prctile) '/'];
        if ~exist([filedir file '_spikes_output.mat'],'file')
            fname_mat = [filedir file '_spikes.mat'];
        else
            fname_mat = [filedir file '_spikes_output.mat'];
        end
        fname_fig = [filedir file '_spikes.fig'];
    else
        [ mouseid, expname ] = find_mouseIDexpname(list);
        if fsave
            filedir = [data_locn 'Data/' file(1:8) '/Processed/' file '/mcorr_' mcorr_method '/' ...
                        segment_method '_' mouseid '_' expname '/' str_fissa '/bl_prctile' num2str(bl_prctile) '/' ];
            fname_mat = [filedir file '_' mouseid '_' expname '_ref' reffile '_spikes.mat'];
            fname_fig = [filedir file '_' mouseid '_' expname '_ref' reffile '_spikes.fig'];
        else
            fname_mat = [];
        end
    end    

    if ~isempty(ddf_f)
        C = ddf_f;
    else
        C = df_f;
    end

    if force || ~exist(fname_mat,'file')
        if isempty(list)
            str = sprintf( '%s: Doing spike extraction\n', file );
        else
            if fsave
                str = sprintf( '%s: Doing spike extraction\n', [mouseid '_' expname '_' file] );
            else
                str = sprintf( '%s: Doing spike extraction\n', [mouseid '_' expname] );
            end
        end
        cprintf(str)

        spikes = extractSpikes( C, params.spkExtract );
        
        if fsave
            % Save output
            spike_output.spikes = spikes;
            if ~exist(filedir,'dir'), mkdir(filedir); end
            spike_output.params = params.spkExtract;
            save(fname_mat,'-struct','spike_output');

            % Make and save plot
            plotSpikes(spikes, fname_fig(1:end-4));
        end
        
        if isempty(list)
            currstr = sprintf( '%s: Spike extraction done\n', file );
        else
            if fsave
                currstr = sprintf( '%s: Spike extraction done\n', [mouseid '_' expname '_' file] );
            else
                currstr = sprintf( '%s: Spike extraction done\n', [mouseid '_' expname] );
            end
        end
        refreshdisp(currstr,str)
    else
        spike_output = load(fname_mat);
        spikes = spike_output.spikes;
        params.spkExtract = spike_output.params;
        
        if fsave && ~exist(fname_fig,'file')
            plotSpikes(spikes, fname_fig(1:end-4));
        end

        if isempty(list)
            fprintf( '%s: Spike extraction data found and loaded\n', file );
        else
            if ~isempty(file)
                fprintf( '%s: Spike extraction data found and loaded\n', [mouseid '_' expname '_' file] );
            else
                fprintf( '%s: Spike extraction data found and loaded\n', [mouseid '_' expname] );
            end
        end
    end

end % function
