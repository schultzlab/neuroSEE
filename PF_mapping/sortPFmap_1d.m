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

function [ sorted_pfMap, sortIdx ] = sortPFmap_1d( pfMap, infoMap, Nepochs)

    if nargin<3, Nepochs = 1; end
    
    n_info = 0; % ratio of units to consider as not informative
    info_type = 1; % 1 is info/sec, 2 is info/spk
    Npcs = size(pfMap,1);
    Nbins = size(pfMap,2);

    % Initialise arrays
    sorted_pfMap = zeros(Npcs,Nbins,Nepochs);
    sortIdx = zeros(Npcs,Nepochs);
    
    for e = 1:Nepochs
        mat = squeeze(pfMap(:,:,e));                                 % Ncells x Nbins
        % sort by information
        [~,idx] = sort( squeeze(infoMap(:,:,e)) );                      % Ncells x 2
        idx_info = idx( round(n_info*Npcs)+1:Npcs, info_type );     % info idx (<Ncells x 1)
        idx_no = idx( 1:round(n_info*Npcs), info_type );              % no info idx (<Ncells x 1)
        pf_info = mat(idx_info,:);                                      % subset of placeMap with info
        maxloc = zeros(length(idx_info),1);
        for i = 1:length(idx_info)
            [~,maxloc(i)] = max(pf_info(i,:));
        end
        [~,sIdx] = sort(maxloc);                                % sort according to phi of max info
        sorted_pfMap(:,:,e) = [ pf_info(sIdx,:); mat(flip(idx_no),:) ]; % cat 2 sorted subsamples
        
        sortIdx(:,e) = [idx_info(sIdx); flip(idx_no)];
    end

end
