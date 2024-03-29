function hSIstruct = sortPFmaps(rateMap, rateMap_sm, normrateMap_sm, pfLoc, hSIstruct)
    % Sort place field maps
    if ~iscell(rateMap)
        if ~isempty(hSIstruct.pcIdx)
            [ ~, sortIdx ] = sort( pfLoc(hSIstruct.pcIdx) );
            hSIstruct.sortpcIdx = hSIstruct.pcIdx(sortIdx);
            hSIstruct.sort_pfMap = rateMap(hSIstruct.sortpcIdx,:);
            hSIstruct.sort_normpfMap_sm = normrateMap_sm(hSIstruct.sortpcIdx,:);
            if ~isempty(rateMap_sm)
                hSIstruct.sort_pfMap_sm = rateMap_sm(hSIstruct.sortpcIdx,:);
            end
        end
    else
        for e = 1:length(rateMap)
            if ~isempty(hSIstruct.pcIdx{e})
                [ ~, sortIdx ] = sort( pfLoc{e}(hSIstruct.pcIdx{e}) );
                hSIstruct.sortpcIdx{e} = hSIstruct.pcIdx{e}(sortIdx);
                hSIstruct.sort_pfMap{e} = rateMap{e}(hSIstruct.sortpcIdx{e},:);
                hSIstruct.sort_normpfMap_sm{e} = normrateMap_sm{e}(hSIstruct.sortpcIdx{e},:);
                if ~isempty(rateMap_sm)
                    hSIstruct.sort_pfMap_sm{e} = rateMap_sm{e}(hSIstruct.sortpcIdx{e},:);
                end
            end
        end
    end
end

