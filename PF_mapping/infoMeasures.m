function [ SIsec, SIspk, MI ] = infoMeasures( placemap, occMap, mask )

if mask % if an env mask or occupational map is used
    mask = occMap>0;
    meanRate = nanmean(nanmean(placemap(mask))); % obtain mean spatial FR
else; meanRate = nanmean(nanmean(placemap));
end
OccMapProb = occMap ./ nansum(nansum(occMap)); %obtain spatial prob density
information = 0; % initialise information estimation

% Skaggs spatial info estimation
for bini=1:size(OccMapProb, 1)
    for binj=1:size(OccMapProb, 2)
        if (placemap(bini,binj) > 0) % loop over every 
            information = information + placemap(bini,binj)*...
                log2(placemap(bini,binj)/meanRate)*OccMapProb(bini,binj);       
        end
    end
end
SIsec = information;
SIspk = SIsec/meanRate;

% mutual info calculations
% pR = 0; % p(r)
% for bini=1:size(OccMapProb, 1)
%     for binj=1:size(OccMapProb, 2)
%         if (placemap(bini,binj) > 0) % loop over every non-0 bin
%             pR = pR + placemap(bini,binj)*OccMapProb(bini,binj);
%         end
%     end
% end
% MI = 0;
% for bini=1:size(OccMapProb, 1)
%     for binj=1:size(OccMapProb, 2)
%         if (placemap(bini, binj) > 0) % loop over every non-0 bin
%             MI = MI + OccMapProb(bini,binj)*placemap(bini,binj)*...
%                 log2(placemap(bini,binj)/pR);
%         end
%     end
% end

%spkRateMap = discretize(placemap,4);
%MI = mutInfo(spkRateMap, occMap);

MI = SIsec;
end


