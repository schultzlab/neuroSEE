% Written by Ann Go
%
% This function generates place field maps if they don't yet exist and then
% sorts them according to location of maximum Skagg's information. It then
% saves the data into a mat file and generates summary plots.
%
%% INPUTS
%   spikes      : ncell rows 
%   trackData   : cell of tracking data with fields x, y, phi, speed, time
%   data_locn   : repository of GCaMP data
%   file        : part of filename of video to be processed in the format
%                   yyyymmdd_hh_MM_ss
%   params      : settings
%
%% OUTPUTS
%   occMap      : occupancy map
%   spikeMap    : spike map (Ncells rows x Nbins columns)
%   infoMap     : information map 
%   downData    : tracking data downsampled to imaging frequency, fields are
%                   x, y, phi, speed, time
%   activeData  : downsampled tracking data for when animal was moving, fields are
%                   x, y, phi, speed, time, spikes

% place field maps
%   placeMap    : place field map obtained with histogram estimation (same size as spikeMap)
%   placeMap_smooth : smoothed version of placeMap 
%   placeMap_asd    : place field map obtained with ASD

%   params

%% OUTPUTS for '1D' only 
% sorted according to location of maximum mutual information
%   sorted_placeMap          
%   sorted_placeMap_smooth  
%   sorted_placeMap_asd    

% normalised and sorted according to location of maximum mutual information
%   normsorted_placeMap          
%   normsorted_placeMap_smooth  
%   normsorted_placeMap_asd    

% per trial place field maps
%   placeMap_pertrial
%   normplaceMap_pertrial

%   PCidx   : row indices of spikes corresponding to place cells
%   sortIdx : sorted row indices corresponding to sorted_placeMap

function varargout = neuroSEE_mapPF( spikes, trackData, data_locn, file, params, force)
    
    if nargin<6, force = 0; end

    mcorr_method = params.methods.mcorr_method;
    segment_method = params.methods.segment_method;
    dofissa = params.methods.dofissa;
    
    if params.methods.dofissa
        str_fissa = 'FISSA';
    else
        str_fissa = 'noFISSA';
    end
    filedir = [data_locn,'Data/',file(1:8),'/Processed/',file,'/mcorr_',mcorr_method,'/',segment_method,'/',str_fissa,'/'];
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

        if strcmpi(params.mode_dim,'1d')
            % Generate place field maps
            [occMap, spikeMap, infoMap, infoMap_asd,...
             placeMap, placeMap_smooth, placeMap_asd, ...
             placeMap_pertrial, normplaceMap_pertrial, PCidx,...
             downData, activeData]...
               = generatePFmap_1D( spikes, imtime, trackData, params );
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
    end
end


