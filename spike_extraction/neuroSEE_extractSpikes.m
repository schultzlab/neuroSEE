% Written by Ann Go
%
% This function calculates the ratiometric Ca time series (R) and extracts
% the spikes
% 
% INPUTS
%   tsG         : time series in green channel
%   tsR         : time series in red channel
%   data_locn   : directory of GCaMP6 data
%   file        : part of filename of 2P image in the format
%                   yyyymmdd_HH_MM_SS
%   params.
%       RsmoothFac   : width of smoothing window
%       g            : OASIS fluorescence impulse factor
%       lambda       : OASIS sparsity penalty
%   force       : if =1, existing R and spike data is overwritten
%
% OUTPUS
%   R           : ratiometric Ca time series
%   spikes      : spikes extracted from R
%   params

function [R, spikes, params] = neuroSEE_extractSpikes( tsG, tsR, data_locn, file, params, force )
    if nargin<6, force = 0; end
        
    filedir = fullfile(data_locn,'Data/',file(1:8),'/Processed/',file,'/');
    fname_mat = [filedir file '_spikes_output.mat'];
    
    % If asked to force overwrite
    if force
        str = sprintf( '%s: Extracting spikes\n', file );
        cprintf(str)

        % Calculate ratiometric Ca time series
        R = ratiometric_Ca( tsG, tsR, params(1).RsmoothFac );
    
        spikes = nndORoasis(R, 2, params(1).g, params(1).lambda); % always use 2 (oasis), 1 is for nnd
%         save(fname_mat,'R','spikes','params');
        
        currstr = sprintf( '%s: Spikes extracted\n', file );
        refreshdisp(currstr,str)
    else
        yn_fname_mat = exist(fname_mat,'file');
        % Find out if spike extraction data already exists. If not,
        % generate data
        if yn_fname_mat
            m = load(fname_mat);
            R = m.R;
            spikes = m.spikes;
            params = m.params;
            
            str = sprintf( '%s: Spike extraction data loaded\n', file );
            cprintf(str)
        else
            str = sprintf( '%s: Extracting spikes\n', file );
            cprintf(str)
            
            % Calculate ratiometric Ca time series
            R = ratiometric_Ca( tsG, tsR, params(1).RsmoothFac );

            spikes = nndORoasis(R, 2, params(1).g, params(1).lambda); % always use 2 (oasis), 1 is for nnd
%             save(fname_mat,'R','spikes','params');
            
            currstr = sprintf( '%s: Spikes extracted\n', file );
            refreshdisp(currstr,str)
        end
    end
end
