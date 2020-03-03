% Adapted from neurosee_PFsort.m by Seigfred Prado
% This script calculates the vector average place preference.

function [ sorted_pfMap, sortIdx, sorted_meanphi, sorted_cv ] = PFsort_vecAve( pfMap )

Nbins = size( pfMap, 2 );
if Nbins > 360
    disp('Error: cannot perform circular statistics.');
else
    i = 1:Nbins:360;
    sin_i = sind(i);
    cos_i = cosd(i);
    
    % Initialise matrices
    alpha = zeros(size(pfMap,1), 1); % mean direction
    cv = zeros(size(pfMap,1), 1); % circular variance
    % y = zeros(size(placeMap,1), 1);
    
    % Sort
    for k = 1:size(pfMap,1)
        alpha(k) = atan2d(sum(pfMap(k,:).*sin_i),sum(pfMap(k,:).*cos_i));
        % y(k,:) = abs(sum(placeMap(k,:).*sin_i)/sum(placeMap(k,:).*cos_i));
        if alpha(k)<0
            alpha(k) = alpha(k)+360;
        end
        cv(k) = 1 - (sqrt(sum(pfMap(k,:).*sin_i)^2 + sum(pfMap(k,:).*cos_i)^2))/360; % circular variance = 1 - ResultantLength
    end
    
%     alpha_placey = find(y>15);
%     [~, sortIdx] = sortrows(alpha_placey);
%     pfMap_plot = pfMap(sortIdx,:);
%     sorted_pfMap = pfMap_plot(:,any(pfMap_plot));
    
    [sorted_meanphi, sortIdx] = sortrows(alpha);   
    sorted_pfMap = pfMap(sortIdx,:);
    sorted_cv = cv(sortIdx);
    
end

end

