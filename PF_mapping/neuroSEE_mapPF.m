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

function [ occMap, hist, asd, activeTrackdata, params, spkMap, spkIdx ] = neuroSEE_mapPF( spikes, downTrackdata, data_locn, file, params, force, list, reffile)
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
    
    filedir = [data_locn,'Data/',file(1:8),'/Processed/',file,'/mcorr_',mcorr_method,'/',segment_method,'/',str_fissa,'/PFmaps/'];
    if ~exist(filedir,'dir'), mkdir(filedir); end
    fname_mat = [filedir file '_PFmap_output.mat'];
    
    if force || ~exist(fname_mat,'file')
        str = sprintf( '%s: Generating place field maps\n', file );
        cprintf(str)

        % If imaging timestamps exist, use them. If not, generate timestamps from
        % known scanning frame rate.
%         dir_timestamps = [data_locn 'Data/' file(1:8) '/Timestamps/'];
%         if exist(dir_timestamps,'dir')
%             imtime = extractImtime(dir_timestamps);
%         else
            imtime = [];
%         end
        
        Nepochs = params.PFmap.Nepochs;
        if strcmpi(params.mode_dim,'1D')
            % Generate place field maps
            [occMap, hist, asd, downData, activeData] = generatePFmap_1d(spikes, imtime, downTrackdata, params);
           
            % If 1D, sort place field maps 
            [ hist.sort_pfMap, hist.sortIdx ] = sortPFmap_1d( hist.pfMap, hist.infoMap, Nepochs );
            [ asd.sort_pfMap, asd.sortIdx ] = sortPFmap_1d( asd.pfMap, asd.infoMap, Nepochs );
            for en = 1:Nepochs
                hist.sort_pfMap_sm(:,:,en) = hist.pfMap_sm(hist.sortIdx(:,en),:,en);
                hist.sort_normpfMap(:,:,en) = hist.normpfMap(hist.sortIdx(:,en),:,en);
                hist.sort_normpfMap_sm(:,:,en) = hist.normpfMap_sm(hist.sortIdx(:,en),:,en);
            
                asd.sort_pfMap(:,:,en) = asd.pfMap(asd.sortIdx(:,en),:,en);
                asd.sort_normpfMap(:,:,en) = asd.normpfMap(asd.sortIdx(:,en),:,en);
            end
            
            % Make plots
            makeplot_1d(occMap, hist, asd);
        
            % Save output
            output.occMap = occMap;
            output.hist = hist;
            output.asd = asd;
            output.downData = downData;
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
        occMap = m.occMap;
        hist = m.hist;
        asd = m.asd;
        downData = m.downData;
        activeData = m.activeData;
        params.PFmap = m.params;
        Nepochs = params.PFmap.Nepochs;
        if isfield(m,'spkMap')
            spkMap = m.spkMap;
        end
        if isfield(m,'spkIdx')
            spkIdx = m.spkIdx;
        end
        
        str = sprintf( '%s: Place field map data loaded\n', file );
        cprintf(str)
    end
    
    function makeplot_1d(occMap, hist, asd)
        Npcs = length(hist.pcIdx);
        Npcs_asd = length(asd.pcIdx);
        
        % summary of occMap, spkMaps, pfMaps
        for e = 1:Nepochs
            fh = figure('Position',[1087 648 800 800]);
            map = viridisMap; colormap(map);
            subplot(10,8,2:4); imagesc(occMap(e,:));
                xticks([]); yticks([]); ylabel('Occ');
                title('Histogram estimation'); colorbar;
            subplot(10,8,6:8); imagesc(occMap(e,:));
                xticks([]); yticks([]); ylabel('Occ');
                title('ASD'); colorbar;

            subplot(10,8,[10:12,18:20,26:28]);
                imagesc(hist.spkMap(hist.sortIdx(:,e),:,e));
                xticks([]);
                yticks([1 Npcs]); ylabel('Cell #'); 
                title('Spike map'); colorbar;
            subplot(10,8,[14:16,22:24,30:32]);
                imagesc(asd.spkMap(asd.sortIdx(:,e),:,e));
                xticks([]);
                yticks([1 Npcs_asd]);  ylabel('Cell #'); 
                title('Spike map'); colorbar;
            
            subplot(10,8,[33,41,49]); imagesc(hist.infoMap(hist.sortIdx,1,e));
                xticks([]);
                yticks([1 Npcs]); ylabel('Cell #'); 
                title('Max MI'); colorbar;
            subplot(10,8,[34:36,42:44,50:52]);    
                imagesc(hist.sort_pfMap(:,:,e)); 
                xticks([]); yticks([1 Npcs]);
                title('Place field map'); colorbar;
            subplot(10,8,[37,45,53]); imagesc(asd.infoMap(asd.sortIdx,1,e));
                xticks([]);
                yticks([1 Npcs_asd]); 
                title('Max MI'); colorbar;
            subplot(10,8,[38:40,46:48,54:56]);    
                imagesc(asd.sort_pfMap(:,:,e)); 
                yticks([1 Npcs_asd]);
                xticks([1 15 30]); xticklabels([1 50 100]);
                xlabel('Position (cm)');
                title('Place field map'); colorbar;

            subplot(10,8,[58:60,66:68,74:76]);    
                imagesc(hist.sort_pfMap_sm(:,:,e)); 
                yticks([1 Npcs]); ylabel('Cell #');
                xticks([1 15 30]); xticklabels([1 50 100]);
                xlabel('Position (cm)');
                title('Smoothened pf map'); colorbar; 

            if Nepochs == 1
                fname_fig = [filedir file '_PFmaps'];
            else
                fname_fig = [filedir file '_PFmaps_' num2str(e) 'of' num2str(Nepochs) 'ep'];
            end
            savefig( fh, fname_fig );
            saveas( fh, fname_fig, 'png' );
            close( fh );
        end
        
        % summary of occMap, spkMaps, normpfMaps
        for e = 1:Nepochs
            fh = figure('Position',[1087 648 800 800]);
            map = viridisMap; colormap(map);
            subplot(10,8,2:4); imagesc(occMap(e,:));
                xticks([]); yticks([]); ylabel('Occ');
                title('Histogram estimation'); colorbar;
            subplot(10,8,6:8); imagesc(occMap(e,:));
                xticks([]); yticks([]); ylabel('Occ');
                title('ASD'); colorbar;

            subplot(10,8,[10:12,18:20,26:28]);
                imagesc(hist.normspkMap(hist.sortIdx(:,e),:,e));
                xticks([]);
                yticks([1 Npcs]); ylabel('Cell #'); 
                title('Normalised spk map'); colorbar;
            subplot(10,8,[14:16,22:24,30:32]);
                imagesc(asd.normspkMap(asd.sortIdx(:,e),:,e));
                xticks([]);
                yticks([1 Npcs_asd]);  ylabel('Cell #'); 
                title('Normalised spk map'); colorbar;
            
            subplot(10,8,[33,41,49]); imagesc(hist.infoMap(hist.sortIdx,1,e));
                xticks([]);
                yticks([1 Npcs]); ylabel('Cell #'); 
                title('Max MI'); colorbar;
            subplot(10,8,[34:36,42:44,50:52]);    
                imagesc(hist.sort_normpfMap(:,:,e)); 
                xticks([]); yticks([1 Npcs]);
                title('Normalised pf map'); colorbar;
            subplot(10,8,[37,45,53]); imagesc(asd.infoMap(asd.sortIdx,1,e));
                xticks([]);
                yticks([1 Npcs_asd]); 
                title('Max MI'); colorbar;
            subplot(10,8,[38:40,46:48,54:56]);    
                imagesc(asd.sort_normpfMap(:,:,e)); 
                yticks([1 Npcs_asd]);
                xticks([1 15 30]); xticklabels([1 50 100]);
                xlabel('Position (cm)');
                title('Normalised pf map'); colorbar;

            subplot(10,8,[58:60,66:68,74:76]);    
                imagesc(hist.sort_normpfMap_sm(:,:,e)); 
                yticks([1 Npcs]); ylabel('Cell #');
                xticks([1 15 30]); xticklabels([1 50 100]);
                xlabel('Position (cm)');
                title('Norm smooth pf map'); colorbar; 

            if Nepochs == 1
                fname_fig = [filedir file '_normPFmaps'];
            else
                fname_fig = [filedir file '_normPFmaps_' num2str(e) 'of' num2str(Nepochs) 'ep'];
            end
            savefig( fh, fname_fig );
            saveas( fh, fname_fig, 'png' );
            close( fh );
        end
        
        % per trial spike maps
        [nRow, nCol] = getnRownCol(Npcs);
        nPlot = nRow*nCol;
        
        % histogram
        Ntrials = size(hist.spkMap_pertrial,1);
        for ii=0:Npcs/nPlot
            fh = figure;
            map = viridisMap; colormap(map);
            ha = tight_subplot(nRow,nCol,[.01 .01],[.01 .05],[.01 .01]);
            for jj=0:nPlot-1
                if (ii*nPlot+jj+1) <= Npcs
                    axes(ha(+jj+1));
                    imagesc(hist.spkMap_pertrial(:,:,ii*nPlot+jj+1)); 
                    yticks(1:Ntrials:Ntrials); yticklabels([1,Ntrials]); ylabel('Trial #');
                    xticks([1 15 30]); xticklabels([1 50 100]); xlabel('Position (cm)');
                    axis off; title(['Cell ' num2str(ii*nPlot+jj+1)],'fontsize',15);
                end
            end
            if Npcs/nPlot <= 1
                fname_fig = [filedir file '_spk_pertrial_hist'];
            else
                fname_fig = [filedir file '_spk_pertrial_hist_' num2str(ii+1)];
            end
            savefig( fh, fname_fig );
            saveas( fh, fname_fig, 'png' );
            close( fh );
        end 

        for ii=0:Npcs/nPlot
            fh = figure;
            map = viridisMap; colormap(map);
            ha = tight_subplot(nRow,nCol,[.01 .01],[.01 .05],[.01 .01]);
            for jj=0:nPlot-1
                if (ii*nPlot+jj+1) <= Npcs
                    axes(ha(+jj+1));
                    imagesc(hist.normspkMap_pertrial(:,:,ii*nPlot+jj+1)); 
                    yticks(1:Ntrials:Ntrials); yticklabels([1,Ntrials]); ylabel('Trial #');
                    xticks([1 15 30]); xticklabels([1 50 100]); xlabel('Position (cm)');
                    axis off; title(['Cell ' num2str(ii*nPlot+jj+1)],'fontsize',15);
                end
            end
            if Npcs/nPlot <= 1
                fname_fig = [filedir file '_normspk_pertrial_hist'];
            else
                fname_fig = [filedir file '_normspk_pertrial_hist_' num2str(ii+1)];
            end
            savefig( fh, fname_fig );
            saveas( fh, fname_fig, 'png' );
            close( fh );
        end 
        
        % asd
        [nRow, nCol] = getnRownCol(Npcs_asd);
        nPlot = nRow*nCol;

        Ntrials = size(asd.spkMap_pertrial,1);
        for ii=0:Npcs_asd/nPlot
            fh = figure;
            map = viridisMap; colormap(map);
            ha = tight_subplot(nRow,nCol,[.01 .01],[.01 .05],[.01 .01]);
            for jj=0:nPlot-1
                if (ii*nPlot+jj+1) <= Npcs_asd
                    axes(ha(+jj+1));
                    imagesc(asd.spkMap_pertrial(:,:,ii*nPlot+jj+1)); 
                    yticks(1:Ntrials:Ntrials); yticklabels([1,Ntrials]); ylabel('Trial #');
                    xticks([1 15 30]); xticklabels([1 50 100]); xlabel('Position (cm)');
                    axis off; title(['Cell ' num2str(ii*nPlot+jj+1)],'fontsize',15);
                end
            end
            if Npcs_asd/nPlot <= 1
                fname_fig = [filedir file '_spk_pertrial_asd'];
            else
                fname_fig = [filedir file '_spk_pertrial_asd_' num2str(ii+1)];
            end
            savefig( fh, fname_fig );
            saveas( fh, fname_fig, 'png' );
            close( fh );
        end 

        for ii=0:Npcs_asd/nPlot
            fh = figure;
            map = viridisMap; colormap(map);
            ha = tight_subplot(nRow,nCol,[.01 .01],[.01 .05],[.01 .01]);
            for jj=0:nPlot-1
                if (ii*nPlot+jj+1) <= Npcs_asd
                    axes(ha(+jj+1));
                    imagesc(asd.normspkMap_pertrial(:,:,ii*nPlot+jj+1)); 
                    yticks(1:Ntrials:Ntrials); yticklabels([1,Ntrials]); ylabel('Trial #');
                    xticks([1 15 30]); xticklabels([1 50 100]); xlabel('Position (cm)');
                    axis off; title(['Cell ' num2str(ii*nPlot+jj+1)],'fontsize',15);
                end
            end
            if Npcs_asd/nPlot <= 1
                fname_fig = [filedir file '_normspk_pertrial_asd'];
            else
                fname_fig = [filedir file '_normspk_pertrial_asd_' num2str(ii+1)];
            end
            savefig( fh, fname_fig );
            saveas( fh, fname_fig, 'png' );
            close( fh );
        end 
        
        % remapping within a session
        if Nepochs > 1
            fh = figure;
            map = viridisMap; colormap(map);
            for ei = 1:Nepochs % rows: sorting
                for ej = 1:Nepochs % cols: epochs 
                    subplot(Nepochs, Nepochs, (ei-1)*Nepochs + ej); imagesc(hist.normpfMap(hist.sortIdx(:,ei),:,ej)); 
                    title(['Epoch ' num2str(ej)]); ylabel(['Epoch' num2str(ei) ' sorting']);
                end
            end
            fname_fig = [filedir file '_remapping_hist'];
            savefig( fh, fname_fig );
            saveas( fh, fname_fig, 'png' );
            close( fh );
            
            fh = figure;
            map = viridisMap; colormap(map);
            for ei = 1:Nepochs % rows: sorting
                for ej = 1:Nepochs % cols: epochs 
                    subplot(Nepochs, Nepochs, (ei-1)*Nepochs + ej); imagesc(asd.normpfMap(asd.sortIdx(:,ei),:,ej)); 
                    title(['Epoch ' num2str(ej)]); ylabel(['Epoch' num2str(ei) ' sorting']);
                end
            end
            fname_fig = [filedir file '_remapping_asd'];
            savefig( fh, fname_fig );
            saveas( fh, fname_fig, 'png' );
            close( fh );
        end
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


