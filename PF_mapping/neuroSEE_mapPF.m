% Written by Ann Go
%
% This function generates place field maps if they don't yet exist and then
% sorts them according to location of maximum Skagg's information. It then
% saves the data into a mat file and generates summary plots.
%
% INPUTS
%   spikes      : ncell rows 
%   trackdata   : cell of tracking data with fields x, y, phi, speed, time
%   data_locn   : repository of GCaMP data
%   file        : part of filename of video to be processed in the format
%                   yyyymmdd_hh_MM_ss
%   params.
%       mode_dim    : 1 for 1D, 2 for 2D 
%       mode_method : 1 for ASD, 2 for histogram estimation
%       Nbins       : number of location bins
%       Nepochs     : number of epochs for each 4 min video
%       histsmoothFac : Gaussian smoothing window for histogram estimation
%       Vthr        : speed threshold (mm/s) Note: David Dupret uses 20 mm/s, Neurotar uses 8 mm/s
%
% OUTPUTS
%   occMap      : occupancy map, 1D: Ncells rows x Nbins columns
%                                2D: Nbins rows, Nbins, columns, Ncells stack
%   spikesMap   : spike map (same size as occMap)
%   infoMap     : information map (same size as occMap)
%   placeMap    : place field map (same size as occMap)
%   downData    : tracking data downsampled to imaging frequency, fields are
%                   x, y, phi, speed, time
%   activaData  : downsampled tracking data for when animal was moving, fields are
%                   x, y, phi, speed, time, spikes
%   placeMap_smooth : smoothed version of placeMap (for when histogram estimation is used)
%   sorted_placeMap : placeMap_smooth (or placeMap if no smoothed version)
%                       sorted by max location of Skagg's information
%   sortIdx         : sorted indices corresponding to sorted_placeMap
%   params


function [occMap, spikeMap, infoMap, placeMap, downData, activeData, placeMap_smooth, sorted_placeMap, ...
            normsorted_placeMap, sortIdx, params]= neuroSEE_mapPF( spikes, trackdata, data_locn, file, params, force)
    
    if nargin<6, force = 0; end
    
    filedir = fullfile(data_locn,'Data/',file(1:8),'/Processed/',file,'/');
    fname_mat = [filedir file '_PFmapping_output.mat'];
    fname_fig1 = [filedir file '_2P_PFmap_summary_norm.fig'];
    fname_fig2 = [filedir file '_2P_PFmap_summary.fig'];
    
    % If asked to force overwrite
    if force
        str = sprintf( '%s: Generating place field maps\n', file );
        cprintf(str)

        % If imaging timestamps exist, use them. If not, generate timestamps from
        % known scanning frame rate.
        dir_timestamps = [data_locn 'Data/' file(1:8) '/Timestamps/'];
        yn_dir_timestamps = exist(dir_timestamps,'dir');

%         if yn_dir_timestamps
%             imtime = extractImtime(dir_timestamps);
%         else
            imtime = [];
%         end
    
        % Generate place field maps
        imrate = params.imrate;
        Vthr = params.Vthr;
        mode_dim = params.mode_dim;
        mode_method = params.mode_method;
        Nbins = params.Nbins;
        Nepochs = params.Nepochs;
        histsmoothFac = params.histsmoothFac;
        [ occMap, spikeMap, infoMap, placeMap, downData, activeData, placeMap_smooth ] = ...
            generatePFmap( spikes, imtime, trackdata, imrate, Vthr, mode_dim, mode_method, Nbins, Nepochs, histsmoothFac );
        
        % If 1D, sort place field maps 
        if mode_dim == 1
            [ sorted_placeMap, normsorted_placeMap, sortIdx ] = sortPFmap( placeMap_smooth, infoMap, Nepochs);
            if Nepochs == 1
                % create summary figure and save it
                fh = figure('Position',[1087 648 500 800]);
                subplot(9,5,2:5); imagesc(occMap);
                    % xticks(Nbins/6:Nbins/6:Nbins); yticks([]);
                    xticks([]); yticks([]);
                    title('Occupancy map'); % colorbar;
                subplot(9,5,[6,11,16,21]); imagesc(infoMap(sortIdx,2));
                    xticks([]); yticks([]); ylabel('Cells'); 
                    title('Info map'); colorbar;
                subplot(9,5,[7:10,12:15,17:20,22:25]);
                    nspikeMap = spikeMap./max(max(spikeMap));
                    imagesc(nspikeMap(sortIdx,:)); xticks([]); yticks([]);
                    title('Spike maps');
                subplot(9,5,[27:30,32:35,37:40,42:45]);    
                    imagesc(normsorted_placeMap); 
                    xticks(Nbins/6:Nbins/6:Nbins); yticks([]);
                    degperbin = 360/Nbins; xticklabels(degperbin*(Nbins/6:Nbins/6:Nbins));
                    % colorbar; % caxis([0,0.005]);
                    %xticks([1 30 60 90 120 150 180]); xticklabels([1 60 120 180 240 300 360]); 
                    % yticklabels(sortIdx);
                    title('Normalised place field maps');
                    xlabel('Position (degrees)'); ylabel('Cells');
                savefig( fh, fname_fig1 );
                saveas( fh, fname_fig1(1:end-4), 'pdf' );
                close( fh );
                
                fh2 = figure('Position',[1087 648 500 800]);
                subplot(9,5,2:5); imagesc(occMap);
                    % xticks(Nbins/6:Nbins/6:Nbins); yticks([]);
                    xticks([]); yticks([]);
                    title('Occupancy map'); % colorbar;
                subplot(9,5,[6,11,16,21]); imagesc(infoMap(sortIdx,2));
                    xticks([]); yticks([]); ylabel('Cells'); 
                    title('Info map'); colorbar;
                subplot(9,5,[7:10,12:15,17:20,22:25]);
                    nspikeMap = spikeMap./max(max(spikeMap));
                    imagesc(nspikeMap(sortIdx,:)); xticks([]); yticks([]);
                    title('Spike maps');
                subplot(9,5,[27:30,32:35,37:40,42:45]);    
                    imagesc(sorted_placeMap); 
                    xticks(Nbins/6:Nbins/6:Nbins); yticks([]);
                    degperbin = 360/Nbins; xticklabels(degperbin*(Nbins/6:Nbins/6:Nbins));
                    % colorbar; % caxis([0,0.005]);
                    %xticks([1 30 60 90 120 150 180]); xticklabels([1 60 120 180 240 300 360]); 
                    % yticklabels(sortIdx);
                    title('Place field maps');
                    xlabel('Position (degrees)'); ylabel('Cells');
                    colorbar('Location','southoutside');
                    caxis([0, 0.012]);
                savefig( fh2, fname_fig2 );
                saveas( fh2, fname_fig2(1:end-4), 'pdf' );
                close( fh2 );
            end
        else
            sorted_placeMap = [];
            normsorted_placeMap = [];
            sortIdx = [];
        end
        
        % Save
        save(fname_mat,'occMap','spikeMap','infoMap','placeMap','downData', 'activeData', ...
                        'placeMap_smooth','sorted_placeMap','normsorted_placeMap','sortIdx','params');
                    
        currstr = sprintf( '%s: Place field maps generated\n', file );
        refreshdisp(currstr,str)
    else
        yn_fname_mat = exist(fname_mat,'file');
        if yn_fname_mat
            m = load(fname_mat);
            occMap = m.occMap;
            spikeMap = m.spikeMap;
            infoMap = m.infoMap;
            placeMap = m.placeMap;
            downData = m.downData;
            activeData = m.activeData;
            placeMap_smooth = m.placeMap_smooth;
            sorted_placeMap = m.sorted_placeMap;
            sortIdx = m.sortIdx;
            params = m.params;
            
            str = sprintf( '%s: Place field map data loaded\n', file );
            cprintf(str)
        else
            str = sprintf( '%s: Generating place field maps\n', file );
            cprintf(str)
        
            % If imaging timestamps exist, use them. If not, generate timestamps from
            % known scanning frame rate.
            dir_timestamps = [data_locn 'Data/' file(1:8) '/Timestamps/'];
            yn_dir_timestamps = exist(dir_timestamps,'dir');

            if yn_dir_timestamps
                imtime = extractImtime(dir_timestamps);
            else
                imtime = [];
            end

            % Generate place field maps
            imrate = params.imrate;
            Vthr = params.Vthr;
            mode_dim = params.mode_dim;
            mode_method = params.mode_method;
            Nbins = params.Nbins;
            Nepochs = params.Nepochs;
            histsmoothFac = params.histsmoothFac;
            [ occMap, spikeMap, infoMap, placeMap, downData, activeData, placeMap_smooth ] = ...
                generatePFmap( spikes, imtime, trackdata, imrate, Vthr, mode_dim, mode_method, Nbins, Nepochs, histsmoothFac );
            
            % If 1D, sort place field maps 
            if mode_dim == 1
                [ sorted_placeMap, normsorted_placeMap, sortIdx ] = sortPFmap( placeMap_smooth, infoMap, Nepochs);
                if Nepochs == 1
                    % create summary figure and save it
                    fh = figure('Position',[1087 648 500 800]);
                    subplot(9,5,2:5); imagesc(occMap);
                        % xticks(Nbins/6:Nbins/6:Nbins); yticks([]);
                        xticks([]); yticks([]);
                        title('Occupancy map'); % colorbar;
                    subplot(9,5,[6,11,16,21]); imagesc(infoMap(sortIdx,2));
                        xticks([]); yticks([]); ylabel('Cells'); 
                        title('Info map'); colorbar;
                    subplot(9,5,[7:10,12:15,17:20,22:25]);
                        nspikeMap = spikeMap./max(max(spikeMap));
                        imagesc(nspikeMap(sortIdx,:)); xticks([]); yticks([]);
                        title('Spike maps');
                    subplot(9,5,[27:30,32:35,37:40,42:45]);    
                        imagesc(normsorted_placeMap); 
                        xticks(Nbins/6:Nbins/6:Nbins); yticks([]);
                        degperbin = 360/Nbins; xticklabels(degperbin*(Nbins/6:Nbins/6:Nbins));
                        % colorbar; % caxis([0,0.005]);
                        %xticks([1 30 60 90 120 150 180]); xticklabels([1 60 120 180 240 300 360]); 
                        % yticklabels(sortIdx);
                        title('Normalised place field maps');
                        xlabel('Position (degrees)'); ylabel('Cells');
                    savefig( fh, fname_fig1 );
                    saveas( fh, fname_fig1(1:end-4), 'pdf' );
                    close( fh );

                    fh2 = figure('Position',[1087 648 500 800]);
                        subplot(9,5,2:5); imagesc(occMap);
                        % xticks(Nbins/6:Nbins/6:Nbins); yticks([]);
                        xticks([]); yticks([]);
                        title('Occupancy map'); % colorbar;
                    subplot(9,5,[6,11,16,21]); imagesc(infoMap(sortIdx,2));
                        xticks([]); yticks([]); ylabel('Cells'); 
                        title('Info map'); colorbar;
                    subplot(9,5,[7:10,12:15,17:20,22:25]);
                        nspikeMap = spikeMap./max(max(spikeMap));
                        imagesc(nspikeMap(sortIdx,:)); xticks([]); yticks([]);
                        title('Spike maps');
                    subplot(9,5,[27:30,32:35,37:40,42:45]);    
                        imagesc(sorted_placeMap); 
                        xticks(Nbins/6:Nbins/6:Nbins); yticks([]);
                        degperbin = 360/Nbins; xticklabels(degperbin*(Nbins/6:Nbins/6:Nbins));
                        % colorbar; % caxis([0,0.005]);
                        %xticks([1 30 60 90 120 150 180]); xticklabels([1 60 120 180 240 300 360]); 
                        % yticklabels(sortIdx);
                        title('Place field maps');
                        xlabel('Position (degrees)'); ylabel('Cells');
                        colorbar('Location','southoutside');
                    savefig( fh2, fname_fig2 );
                    saveas( fh2, fname_fig2(1:end-4), 'pdf' );
                    close( fh2 );
                end
            else
                sorted_placeMap = [];
                normsorted_placeMap = [];
                sortIdx = [];
            end

            % Save
            save(fname_mat,'occMap','spikeMap','infoMap','placeMap','downData', 'activeData', ...
                            'placeMap_smooth','sorted_placeMap','normsorted_placeMap','sortIdx','params');
                    
            currstr = sprintf( '%s: Place field maps generated\n', file );
            refreshdisp(currstr,str)
        end
    end
end


