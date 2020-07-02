function plotpfMaps_2d(activeData, pfMap_h, pfMap_sm_h, pfMap_a, pcIdx, cell_str, fname, fclose)
    if nargin<8, fclose = true; end
    if nargin<7, fname = []; end

    Nepochs = size( pfMap_h, 4 );
    Ncells = length(pcIdx);
    nRow = 4;
    if isempty(pfMap_a)
        nCol = 3;
    else
        nCol = 4;
    end

    for e = 1:Nepochs
        for ii=0:0%ceil(Ncells/nRow)-1 
            if isempty(pfMap_a)
                fh = figure('Position',[680 678 500 550]); 
            end
            ha = tight_subplot(nRow,nCol,[.03 .02],[.02 .05],[.02 .02]);
            for jj=0:nRow-1
                if (ii*nRow+jj+1) <= Ncells
                    axes(ha(jj*nCol+1));
                    z = activeData.spikes(pcIdx(ii*nRow+jj+1),:);
                    hold on; axis off;
                    % no mistake here in the order of plotting x&y, this is
                    % to match the image pixel indexing to matrix row and
                    % column indexing
                    plot(activeData.y,-activeData.x,'Color',[0.6 0.6 0.6],'LineWidth',1.2); 
                    ind = find(z>0);
                    x = activeData.x(ind);
                    y = activeData.y(ind);
                    spikes = z(ind);
                    [spikes_sorted,sort_ind] = sort(spikes);
                    scatter(y(sort_ind),-x(sort_ind),20,spikes_sorted,'filled');
                    title_str = sprintf('%s %g', cell_str, (ii*nRow+jj+1)); %, max(spikes));
                    title(title_str,'fontsize',15);
                    hold off; axis square;

                    axes(ha(jj*nCol+2));
                    imagesc(squeeze(pfMap_h(:,:,ii*nRow+jj+1,e))');
                    axis off; colorbar; % caxis([0 0.06]);
                    if Nepochs >1 
                        title(['Epoch ',num2str(e)],'fontsize',11);
                    end
                    axes(ha(jj*nCol+3)); 
                    imagesc(squeeze(pfMap_sm_h(:,:,ii*nRow+jj+1,e))');
                    axis off; colorbar; % caxis([0 0.005]);
                    if ~isempty(pfMap_a)
                        axes(ha(jj*nRow+4));
                        imagesc(squeeze(pfMap_a(:,:,ii*nRow+jj+1,e))');
                        axis off; colorbar; % caxis([0 0.003]);
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