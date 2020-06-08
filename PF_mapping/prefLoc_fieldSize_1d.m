% Adapted from prefLoc.m by Simon Schultz
% This script finds the preferred location of a place cell by
%   phi_pref1 : calculating vector average place preference
%   phi_pref2 : finding peak of place field map (firing rate/occupancy)

function [ prefLoc, fieldSize, fieldLoc ] = prefLoc_fieldSize_1d( pfMap )

Ncells = size( pfMap, 1 );
Nbins = size( pfMap, 2 );

fieldSize = zeros( Ncells, 1 );
prefLoc = zeros( Ncells, 1 );

for i = 1:Ncells
    % preferred location
    [maxval,maxloc] = max( pfMap(i,:) );
    prefLoc(i) = maxloc;
    
    % field size
    fieldSize(i) = 0;
    % circularly shift array to center the peak
    shift = round(Nbins/2)-maxloc;
    pf_sh = circshift(pfMap(i,:),shift);
    
    % go through array values to the right of the peak
    for j = round(Nbins/2)+1:Nbins
        if pf_sh(j) >= maxval/2
            fieldSize(i) = fieldSize(i) + 1;
        else
            rInd = (j-1)-shift;
            break
        end
    end
    
    % go through array values to the left of the peak
    for j = round(Nbins/2)-1:-1:1
        if pf_sh(j) >= maxval/2
            fieldSize(i) = fieldSize(i) + 1;
        else
            lInd = (j+1)-shift;
            break
        end
    end   
    
    % exact field location
    if lInd >= 1 && rInd <= Nbins
        fieldLoc = (lInd:rInd);
    elseif lInd < 1 && rInd <= Nbins
        fieldLoc = [1:rInd, lInd+Nbins:Nbins];
    elseif lInd >= 1 && rInd > Nbins
        fieldLoc = [1:rInd-Nbins, lInd:Nbins];
    end
end


