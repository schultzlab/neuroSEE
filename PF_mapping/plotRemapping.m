function plotRemapping(normpfMap, sortIdx, fname)
    % remapping within a session
    fh = figure;
    map = viridisMap; 
    for ei = 1:Nepochs % rows: sorting
        for ej = 1:Nepochs % cols: epochs 
            subplot(Nepochs, Nepochs, (ei-1)*Nepochs + ej); imagesc(normpfMap(sortIdx(:,ei),:,ej)); colormap(map);
            title(['Epoch ' num2str(ej)]); ylabel(['Epoch' num2str(ei) ' sorting']);
        end
    end
    if fsave
        savefig( fh, fname );
        saveas( fh, fname, 'png' );
        if fclose, close( fh ); end
    end
end
