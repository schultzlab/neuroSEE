pfactivet_thr = 0.05; 
activetrials_thr = 0.4; 

inc_idx = []; exc_idx = hist.SIspk.nonpcIdx;
Npcs = numel(hist.SIspk.pcIdx);
for n = 1:Npcs
    c = hist.SIspk.pcIdx(n);
    if hist.pf_activet(c) >= pfactivet_thr
        %if pfData.activetrials(c) >= activetrials_thr
           inc_idx = [inc_idx; c];
%         else
%            exc_idx = [exc_idx; c];
%         end
    else
        exc_idx = [exc_idx; c];
    end
end
hist.SIspk.pcIdx = inc_idx;
hist.SIspk.nonpcIdx = exc_idx;

hist.SIspk = sortPFmaps(hist.rateMap, hist.rateMap_sm, hist.normrateMap_sm, hist.pfLoc, hist.SIspk);
plotPF_1d(hist, [], pfData, false)