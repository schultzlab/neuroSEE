% Written by Ann Go
%
% This function generates place field maps if they don't yet exist and then
% sorts them according to location of maximum Skagg's information (or mutual 
% information if this is the chosen metric). It then
% saves the data into a mat file and generates summary plots.
%
%% INPUTS
%   spikes      : ncell rows 
%   downTrackdata   : cell of downsampled tracking data with fields x, y, phi, speed, time
%   data_locn   : repository of GCaMP data
%   file        : part of filename of video to be processed in the format
%                   yyyymmdd_hh_MM_ss
%   params      : settings
%
%% OUTPUTS
%   occMap      : occupancy map
%   spikeMap    : spike map (Ncells rows x Nbins columns)
%   infoMap     : information map 
%   activeTrackdata  : downsampled tracking data for when animal was moving; fields are
%                       x, y, phi, speed, time, spikes, spikes_pc

% place field maps
%   pfMap       : place field map obtained with histogram estimation (same size as spikeMap)
%   pfMap_sm    : smoothed version of pfMap 
%   pfMap_asd   : place field map obtained with ASD
%   normpfMap       : normalised place field map obtained with histogram estimation (same size as spikeMap)
%   normpfMap_sm    : normalised smoothed version of pfMap 
%   normpfMap_asd   : normalised place field map obtained with ASD

%   params

%% OUTPUTS for '1D' only 
% sorted according to location of maximum information
%   sorted_pfMap          
%   sorted_pfMap_sm  
%   sorted_pfMap_asd    

% normalised and sorted according to location of maximum information
%   sorted_normpfMap          
%   sorted_normpMap_sm  
%   sorted_normpfMap_asd    

% per trial spike maps
%   spikeMap_pertrial
%   normspikeMap_pertrial

%   pcIdx   : row indices of spikes corresponding to place cells
%   sortIdx : sorted row indices corresponding to sorted_pfMap

function [ hist, asd, pfData, activeData, params ] = neuroSEE_mapPF( spikes, downTrackdata, data_locn, file, params, force, list, reffile)
    if nargin<8, reffile = []; end
    if nargin<7, list = []; end
    if nargin<6, force = 0; end

    mcorr_method = params.methods.mcorr_method;
    segment_method = params.methods.segment_method;
    
    if params.methods.dofissa
        str_fissa = 'FISSA';
    else
        str_fissa = 'noFISSA';
    end
    
    if isempty(list)
        sdir = [data_locn,'Data/',file(1:8),'/Processed/',file,'/mcorr_',mcorr_method,'/',segment_method,'/',str_fissa,'/PFmaps/'];
        fname_pref = file;
        fname_mat = [sdir fname_pref '_PFmap_output.mat'];
    else
        [ mouseid, expname ] = find_mouseIDexpname(list);
        filedir = [ data_locn 'Analysis/' mouseid '/' mouseid '_' expname '/group_proc/' mcorr_method '_' segment_method '_'...
                    str_fissa '/' mouseid '_' expname '_imreg_ref' reffile '/'];
        fname_pref = [mouseid '_' expname '_ref' reffile];
        fname_mat = [filedir fname_pref '_PFmap_output.mat'];
        sdir = [filedir '/PFdata/'];
    end
    
    if force || ~exist(fname_mat,'file')
        str = sprintf( '%s: Generating place field maps\n', file );
        cprintf(str)

        % If imaging timestamps exist, use them. If not, generate timestamps from
        % known scanning frame rate.
%         dir_timestamps = [data_locn 'Data/' file(1:8) '/Timestamps/'];
%         if exist(dir_timestamps,'dir')
%             imtime = extractImtime(dir_timestamps);
%         else
%            imtime = [];
%         end
        
        Nepochs = params.PFmap.Nepochs;
        if strcmpi(params.mode_dim,'1D')
            % Generate place field maps
            [hist, asd, activeData, pfData] = generatePFmap_1d( spikes, downTrackdata, params );
           
            % Make plots
            plotPF_1d(hist, asd, PFdata, true, true, sdir, fname_pref)
        
            % Save output
            output.hist = hist;
            output.asd = asd;
            output.pfData = pfData;
            output.activeData = activeData;
            output.params = params.PFmap;
            save(fname_mat,'-struct','output');
        else % '2D'
            [occMap, spkMap, spkIdx, hist, asd, downData, activeData] = generatePFmap_2d(spikes, imtime, downTrackdata, params);
            
             % Make plots
            makeplot_2d(spkMap, activeData, hist, asd);
        
            % Save output
            output.occMap = occMap;
            output.spkMap = spkMap;
            output.spkIdx = spkIdx;
            output.hist = hist;
            output.asd = asd;
            output.downData = downData;
            output.activeData = activeData;
            output.params = params.PFmap;
            save(fname_mat,'-struct','output');
        end
              
        currstr = sprintf( '%s: Place field maps generated\n', file );
        refreshdisp(currstr,str)
    else
        m = load(fname_mat);
        hist = m.hist;
        asd = m.asd;
        pfData = m.pfData;
        activeData = m.activeData;
        params.PFmap = m.params;
        Nepochs = params.PFmap.Nepochs;
        
        str = sprintf( '%s: Place field map data loaded\n', file );
        cprintf(str)
    end
    

    function makeplot_2d(spkMap, activeData, hist, asd)
        Nspk = size(spkMap,3);
        nPlot = 4;
        for e = 1:Nepochs
            for ii=0:(Nspk/nPlot)-1 
                fh = figure; 
                map = viridisMap; colormap(map);
                ha = tight_subplot(nPlot,4,[.01 .005],[.01 .07],[.01 .01]);
                for jj=0:3
                    if (ii*nPlot+jj) <= Nspk
                        axes(ha(jj*nPlot+1));
                        z = activeData.spikes(spkIdx(ii*nPlot+jj+1),:);
                        hold on; axis off;
                        plot(activeData.x,-activeData.y); plot(activeData.x(z>0),-activeData.y(z>0),'r.','markersize',10);
                        title(['Cell ',num2str(ii*nPlot+jj+1)],'fontsize',15);
                        axes(ha(jj*nPlot+2));
                        imagesc(squeeze(hist.pfMap(:,:,ii*nPlot+jj+1,e))');
                        axis off; colorbar; % caxis([0 0.06]);
                        if Nepochs >1 
                            title(['Epoch ',num2str(e)],'fontsize',15);
                        end
                        axes(ha(jj*nPlot+3)); 
                        imagesc(squeeze(hist.pfMap_sm(:,:,ii*nPlot+jj+1,e))');
                        axis off; colorbar; % caxis([0 0.005]);
                        axes(ha(jj*nPlot+4));
                        imagesc(squeeze(asd.pfMap(:,:,ii*nPlot+jj+1,e))');
                        axis off; colorbar; % caxis([0 0.003]);
                    end
                end
                if Nspk/nPlot <= 1
                    if Nepochs == 1
                        fname_fig = [filedir file '_PFmaps'];
                    else
                        fname_fig = [filedir file '_PFmaps_' num2str(e) 'of' num2str(Nepochs) 'ep' ];
                    end
                else
                    if Nepochs == 1
                        fname_fig = [filedir file '_PFmaps_' num2str(ii+1)];
                    else
                        fname_fig = [filedir file '_PFmaps_' num2str(ii+1) '_' num2str(e) 'of' num2str(Nepochs) 'ep' ];
                    end
                end
                savefig( fh, fname_fig );
                saveas( fh, fname_fig, 'png' );
                close( fh );
            end 
        end
    end
end % function


