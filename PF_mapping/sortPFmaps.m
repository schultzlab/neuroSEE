function hSIstruct = sortPFmaps(rateMap, rateMap_sm, normrateMap_sm, pfLoc, hSIstruct)
    % Sort place field maps
    if ~iscell(rateMap)
        if ~isempty(hSIstruct.pcIdx)
            [ ~, sortIdx ] = sort( pfLoc(hSIstruct.pcIdx) );
            hSIstruct.sortpcIdx = hSIstruct.pcIdx(sortIdx);
            hSIstruct.sort_pfMap = rateMap(sortIdx,:);
            hSIstruct.sort_normpfMap_sm = normrateMap_sm(sortIdx,:);
            if ~isempty(rateMap_sm)
                hSIstruct.sort_pfMap_sm = rateMap_sm(sortIdx,:);
            end
        end
    else
        for e = 1:length(rateMap)
            if ~isempty(hSIstruct.pcIdx{e})
                [ ~, sortIdx ] = sort( pfLoc{e}(hSIstruct.pcIdx{e}) );
                hSIstruct.sortpcIdx{e} = hSIstruct.pcIdx{e}(sortIdx);
                hSIstruct.sort_pfMap{e} = rateMap{e}(sortIdx,:);
                hSIstruct.sort_normpfMap_sm{e} = normrateMap_sm{e}(sortIdx,:);
                if ~isempty(rateMap_sm)
                    hSIstruct.sort_pfMap_sm{e} = rateMap_sm{e}(sortIdx,:);
                end
            end
        end
    end
end

