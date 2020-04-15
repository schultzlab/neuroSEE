% Written by Ann Go
% This script calculates the field size for a place field map

function pfSize = fieldSize( pfMap )

Ncells = size( pfMap, 1 );
Nbins = size( pfMap, 2 );
pfSize = zeros( Ncells, 1 );

for i = 1:Ncells
    pfSize(i) = 0;
    [maxval,maxloc] = max( pfMap(i,:) );
    
    % circularly shift array to center the peak
    pf_sh = circshift(pfMap(i,:),round(Nbins/2)-maxloc);
    
    % go through array values to the right of the peak
    for j = round(Nbins/2)+1:Nbins
        if pf_sh(j) >= maxval/2
            pfSize(i) = pfSize(i) + 1;
        else
            y3 = pf_sh(j-1);
            y2 = maxval/2;
            y1 = pf_sh(j);
            x3 = j-1;
            x1 = j;
            x2 = (x3*(y2-y1) + x1*(y3-y2))/(y3-y1);
            pfSize(i) = pfSize(i) + x2-(j-1);
            break
        end
    end
    
    % go through array values to the left of the peak
    for j = round(Nbins/2)-1:-1:1
        if pf_sh(j) >= maxval/2
            pfSize(i) = pfSize(i) + 1;
        else
            y3 = pf_sh(j+1);
            y2 = maxval/2;
            y1 = pf_sh(j);
            x3 = j+1;
            x1 = j;
            x2 = (x3*(y2-y1) + x1*(y3-y2))/(y3-y1);
            pfSize(i) = pfSize(i) + (j+1)-x2;
            break
        end
    end
end

pfSize = pfSize*103/Nbins;

