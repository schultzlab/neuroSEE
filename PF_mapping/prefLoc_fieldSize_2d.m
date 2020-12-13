% Written by Ann Go
% prefLoc: centroid of place field
% fieldSize: full width at half max 
% pfBins: position bin indices encompassing place field

function [ prefLoc, fieldSize, pfBins ] = prefLoc_fieldSize_2d( pfMap )

Ncells = size( pfMap, 3 );
Nbins = size( pfMap, 1 );

prefLoc = zeros( Ncells, 1 );
fieldSize = zeros( Ncells, 1 );
pfBins = cell( Ncells, 1 );

area_2d = pi*(32.5/2)^2;

for c = 1:Ncells
    rMap = pfMap(:,:,c);
    rMap(rMap < max(max(rMap))/2) = 0;
    
    % field size
    fieldSize(c,1) = numel(find(rMap > 0)) * (32.5/Nbins)^2; % in cm^2
    fieldSize(c,2) = fieldSize(c,1)/area_2d;                 % in % of arena
    
    % centroid location
    rMap_BW = rMap > 0; 
    se = strel('disk', 1);
    eroded_rMap_BW = imerode(rMap_BW,se); 
    if isempty(find(eroded_rMap_BW, 1))
        rMap_BW2 = rMap_BW;
    else
        rMap_BW2 = imdilate(eroded_rMap_BW,se);
    end
    
    CC = bwconncomp(rMap_BW2);
    
    if CC.NumObjects > 1
        % find blob where maximum value lies
        for i = 1:CC.NumObjects
            if numel(intersect(find(rMap == max(max(rMap))), sub2ind(size(rMap), CC.PixelIdxList{i}))) > 0
                useInd = i;
            end
        end
        stats = regionprops(CC,rMap,'WeightedCentroid');
        prefLoc(c) = sub2ind(size(rMap), round(stats(useInd).WeightedCentroid(1)), round(stats(useInd).WeightedCentroid(2)));
        % pfBins{c} = CC.PixelIdxList{useInd};
        pfBins{c} = find(rMap_BW);
    elseif CC.NumObjects == 1
        stats = regionprops(CC,rMap,'WeightedCentroid');
        prefLoc(c) = sub2ind(size(rMap), round(stats.WeightedCentroid(1)), round(stats.WeightedCentroid(2)));
        % pfBins{c} = CC.PixelIdxList{1};
        pfBins{c} = find(rMap_BW);
    else
        prefLoc(c) = [];
        pfBins{c} = [];
    end
    
    
end
