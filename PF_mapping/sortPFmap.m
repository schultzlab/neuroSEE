% Written by Ann Go (adapted from Giuseppe's PF_ASD_1d.m)
%
% This function sorts the rows of placeMap (if 1D, the planes if 2D)
% according to location of max Skagg's info contained in infoMap
%
% INPUTS
%   placeMap    : place field map (same size as occMap)
%   infoMap     : information map (same size as occMap)
%   Nepochs     : number of epochs for each 4 min video
%
% OUTPUTS
%   sorted_placeMap : placeMap sorted by max location of Skagg's information
%   sortIdx         : sorted indices corresponding to sorted_placeMap

function [ sorted_placeMap, normsorted_placeMap, allIdx ] = sortPFmap( placeMap, infoMap, Nepochs)

    if nargin<3, Nepochs = 1; end
    
    n_info = 1/10; % ratio of units to consider as not informative
    info_type = 2; % 1 is info/sec, 2 is info/spk
    Ncells = size(placeMap,1);
    Nbins = size(placeMap,2);
    hdl = figure;

    % Initialise arrays
    sorted_placeMap = zeros(Ncells,Nbins,Nepochs);
    normsorted_placeMap = zeros(Ncells,Nbins,Nepochs);
    allIdx = zeros(Ncells,Nepochs);
    
    for e = 1:Nepochs
        mat = squeeze(placeMap(:,:,e));                                 % Ncells x Nbins
        % sort by information
        [~,idx] = sort( squeeze(infoMap(:,:,e)) );                      % Ncells x 2
        idx_info = idx( round(n_info*Ncells)+1:Ncells, info_type );     % info idx (<Ncells x 1)
        idx_no = idx( 1:round(n_info*Ncells), info_type );              % no info idx (<Ncells x 1)
        pf_info = mat(idx_info,:);                                      % subset of placeMap with info
        maxloc = zeros(length(idx_info),1);
        for i = 1:length(idx_info)
            [~,maxloc(i)] = max(pf_info(i,:));
        end
        [~,sortIdx(:,e)] = sort(maxloc);                                % sort according to phi of max info
        rawsorted_placeMap = [ pf_info(sortIdx(:,e),:); mat(flip(idx_no),:) ]; % cat 2 sorted subsamples
        
        % normalise for visualisation
        for i = 1:Ncells
            sorted_placeMap(i,:,e) = rawsorted_placeMap(i,:);
            normsorted_placeMap(i,:,e) = rawsorted_placeMap(i,:)/max(rawsorted_placeMap(i,:));
        end
        
        allIdx(:,e) = [idx_info(sortIdx(:,e)); flip(idx_no)];
        
        % plot
        if Nepochs > 1
            figure(hdl); 
            subplot(1,Nepochs,e); imagesc(sorted_placeMap(allIdx,:,e)); 
            % colorbar; % caxis([0,0.005]);
            xticks([1 30 60 90 120 150 180]); xticklabels([1 60 120 180 240 300 360]); 
            yticks(1:2:Ncells); 
            yticklabels(allIdx(1:2:Ncells,e)); 
            str = sprintf('Epoch %g', e);
            title(str); 
            xlabel('Position (degrees)'); ylabel('Cell number');
        else
            figure(hdl); subplot(121);
            imagesc(normsorted_placeMap(:,:)); 
            % colorbar; % caxis([0,0.005]);
            % xticks([1 30 60 90 120 150 180]); xticklabels([1 60 120 180 240 300 360]); 
            % yticks(1:2:Ncells); 
            %yticklabels(allIdx); %(1:2:Ncells,e)); 
            title('Normalised place field maps');
            xlabel('Position (degrees)'); ylabel('Cell number');
            subplot(122);
            imagesc(sorted_placeMap(:,:)); 
            % colorbar; % caxis([0,0.005]);
            % xticks([1 30 60 90 120 150 180]); xticklabels([1 60 120 180 240 300 360]); 
            % yticks(1:2:Ncells); 
            %yticklabels(allIdx); %(1:2:Ncells,e)); 
            title('Place field maps');
            xlabel('Position (degrees)'); ylabel('Cell number');
        end
    end

end
