% Written by Ann Go

function [ prefLoc, fieldSize ] = prefLoc_fieldSize_2d( pfMap )

Ncells = size( pfMap, 3 );
Nbins = size( pfMap, 1 );

prefLoc = zeros( Ncells, 1 );
fieldSize = zeros( Ncells, 1 );

area_2d = pi*(32.5/2)^2;

for id = 1:Ncells
    rMap = pfMap(:,:,id);
    rMap(rMap < max(max(rMap))/2) = 0;
    
    % field size
    fieldSize(id,1) = numel(find(rMap > 0)) * (32.5/Nbins)^2;
    fieldSize(id,2) = fieldSize(id,1)/area_2d;
    
    % centroid location
    rMap_BW = rMap > 0; 
    se = strel('disk', 1);
    eroded_rMap_BW = imerode(rMap_BW,se); 
    rMap_BW2 = imdilate(eroded_rMap_BW,se); 
    
    CC = bwconncomp(rMap_BW2);
    
    if CC.NumObjects > 1
        for i = 1:CC.NumObjects
            if numel(intersect(find(rMap == max(max(rMap))), sub2ind(size(rMap), CC.PixelIdxList{i}))) > 0
                useInd = i;
            end
        end
        stats = regionprops(CC,rMap,'WeightedCentroid');
        prefLoc(id) = sub2ind(size(rMap), round(stats(useInd).WeightedCentroid(1)), round(stats(useInd).WeightedCentroid(2)));
    else
        stats = regionprops(CC,rMap,'WeightedCentroid');
        prefLoc(id) = sub2ind(size(rMap), round(stats.WeightedCentroid(1)), round(stats.WeightedCentroid(2)));
    end
    
end
