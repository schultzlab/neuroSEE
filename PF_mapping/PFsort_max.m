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

function [ sorted_pfMap, sortIdx ] = PFsort_max( pfMap )

    Npcs = size(pfMap,1);
    Nbins = size(pfMap,2);
    Nepochs = size(pfMap,3);

    % Initialise arrays
    sorted_pfMap = zeros(Npcs,Nbins,Nepochs);
    sortIdx = zeros(Npcs,Nepochs);
    
    % Sort
    for e = 1:Nepochs
        mat = squeeze(pfMap(:,:,e));                                    % Ncells x Nbins
        maxloc = zeros(Npcs,1);
        for i = 1:Npcs
            [~,maxloc(i)] = max(mat(i,:));
        end
        [~,sortIdx(:,e)] = sort(maxloc);                                % sort according to phi of max value
        sorted_pfMap(:,:,e) = mat(sortIdx(:,e),:);          
    end

end
