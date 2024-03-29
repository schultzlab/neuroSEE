function plotpfMaps_2d(activeData_all, pfMap_h_all, normpfMap_h_all, normpfMap_sm_h_all, infoMap_h_all, ...
                        normpfMap_a_all, infoMap_a_all, pcIdx, cell_str, info_str, fname, fclose)
    if nargin<12, fclose = true; end
    if nargin<11, fname = []; end

    Nepochs = size( pfMap_h_all, 4 );
    Ncells = length(pcIdx);
    nRow = 4;
    if isempty(normpfMap_a_all)
        nCol = 3;
    else
        nCol = 4;
    end
    
    % find the delineations for the video: find t = 0
    idx_file = find(diff(activeData_all.t) < 0);
    idx_file = [0; idx_file; numel(activeData_all.t)] +1;
    max_spksize = 18;
    spk_col = 'b';

    for e = 1:Nepochs
        for ii=0:ceil(Ncells/nRow)-1 
            if isempty(normpfMap_a_all)
                fh = figure('Position',[680 678 500 550]); 
            end
            ha = tight_subplot(nRow,nCol,[.04 .02],[.02 .05],[.02 .02]);
            for jj=0:nRow-1
                if (ii*nRow+jj+1) <= Ncells
                    axes(ha(jj*nCol+1));
                    hold on; axis([-140 140 -140 140]); axis off;
                    % no mistake here in the order of plotting x&y, this is
                    % to match the image pixel indexing to matrix row and
                    % column indexing
                    for kk = 1: numel(idx_file) - 1
                        plot(activeData_all.y(idx_file(kk):idx_file(kk+1)-1),-activeData_all.x(idx_file(kk):idx_file(kk+1)-1),...
                            'Color',[0.8 0.8 0.8], 'LineWidth',1); axis square; 
                    end
                    z = activeData_all.spikes(pcIdx(ii*nRow+jj+1),:);
                    ind = find(z>0);
                    x = activeData_all.x(ind);
                    y = activeData_all.y(ind);
                    spikes = z(ind);
                    [spkampl_sorted,sort_ind] = sort(spikes);
                    spk_size = spkampl_sorted/max(spkampl_sorted)*max_spksize;
                    scatter(y(sort_ind),-x(sort_ind),spk_size,spk_col,'filled');
                    title_str = sprintf('%s %g', cell_str, (ii*nRow+jj+1)); 
                    title(title_str,'fontsize',12);
                    hold off; 

                    axes(ha(jj*nCol+2));
                    cmap = viridisMap_whitelowest;
                    colormap(cmap);
                    imagesc(squeeze(normpfMap_h_all(:,:,pcIdx(ii*nRow+jj+1),e))');
                    axis off; axis square; % colorbar; 
                    title_str = sprintf('Max %.2f events/s', max(max(pfMap_h_all(:,:,pcIdx(ii*nRow+jj+1))))); 
                    title(title_str,'fontsize',12);
                    
                    axes(ha(jj*nCol+3));
                    imagesc(squeeze(normpfMap_sm_h_all(:,:,pcIdx(ii*nRow+jj+1),e))');
                    axis off; axis square; % colorbar; % caxis([0 0.005]);
                    title_str = sprintf('%.2f %s', infoMap_h_all(pcIdx(ii*nRow+jj+1)), info_str); 
                    title(title_str,'fontsize',12);
                    
                    if ~isempty(normpfMap_a_all)
                        axes(ha(jj*nRow+4));
                        imagesc(squeeze(normpfMap_a_all(:,:,pcIdx(ii*nRow+jj+1),e))');
                        axis off; % colorbar; % caxis([0 0.003]);
                        title_str = sprintf('%.2f %s', infoMap_a_all(pcIdx(ii*nRow+jj+1)), info_str); 
                        title(title_str,'fontsize',12);
                    end
                end
            end

            if ~isempty(fname)
                if Ncells/nRow <= 1
                    if Nepochs == 1
                        fname_fig = fname;
                    else
                        fname_fig = [fname '_' num2str(e) 'of' num2str(Nepochs) 'ep' ];
                    end
                else
                    if Nepochs == 1
                        fname_fig = [fname '_' num2str(ii+1)];
                    else
                        fname_fig = [fname '_' num2str(ii+1) '_' num2str(e) 'of' num2str(Nepochs) 'ep' ];
                    end
                end
                savefig( fh, fname_fig );
                saveas( fh, fname_fig, 'png' );
                if fclose, close( fh ); end
            end
        end 
    end
end