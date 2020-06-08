y = ddf_f(1,:); 
clear pks locs
[pks,locs]= findpeaks(y,'MinPeakHeight',0.5*max(y));
[sortpks, sortIdx ] = sort(pks,'descend');
if length(pks) >= 20
    tpks = sortpks(1:20);
    sortlocs = locs(sortIdx);
    tlocs = sortlocs(1:20);
else
    tpks = sortpks;
    sortlocs = locs(sortIdx);
    tlocs = sortlocs;
end
maxval = 0.05;
for i = 1:size(tpks,2)
    startInd = tlocs(i);
    % go through array values to the left of the peak
    for j = startInd-1:-1:1
        if y(j) <= maxval*max(y)
            leftInd = j;
            break
        end
    end
    % go through array values to the right of the peak
    for j = startInd+1:numel(y)
        if y(j) <= maxval*max(y)
            rightInd = j;
            break
        end
    end
    % find the skewness of the calcium transient
    ts{i} = y(leftInd:rightInd);
    skew(i) = skewness(ts{i},0);
end

% figure; plot(y); hold on; plot(tlocs,0.05*max(y)*ones(size(tpks)),'.'); hold off
figure; 
for i=1:size(tpks,2)
    subplot(4,5,i);
    plot(ts{i}); 
end
