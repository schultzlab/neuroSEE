function plotPF_2d(spkMap, activeData, hist, asd, fsave, fclose)
    Nspk = size(spkMap,3);
    
    nPlot = 4;
    for e = 1:Nepochs
        for ii=0:(Nspk/nPlot)-1 
            fh = figure; 
            ha = tight_subplot(nPlot,4,[.01 .005],[.01 .07],[.01 .01]);
            for jj=0:3
                if (ii*nPlot+jj) <= Nspk
                    axes(ha(jj*nPlot+1));
                    z = activeData.spikes(spkIdx(ii*nPlot+jj+1),:);
                    hold on; axis off;
                    plot(activeData.x,-activeData.y); plot(activeData.x(z>0),-activeData.y(z>0),'r.','markersize',10);
                    title(['Cell ',num2str(ii*nPlot+jj+1)],'fontsize',15);
                    axes(ha(jj*nPlot+2));
                    imagesc(squeeze(hist.pfMap(:,:,ii*nPlot+jj+1,e))');
                    axis off; colorbar; % caxis([0 0.06]);
                    if Nepochs >1 
                        title(['Epoch ',num2str(e)],'fontsize',15);
                    end
                    axes(ha(jj*nPlot+3)); 
                    imagesc(squeeze(hist.pfMap_sm(:,:,ii*nPlot+jj+1,e))');
                    axis off; colorbar; % caxis([0 0.005]);
                    axes(ha(jj*nPlot+4));
                    imagesc(squeeze(asd.pfMap(:,:,ii*nPlot+jj+1,e))');
                    axis off; colorbar; % caxis([0 0.003]);
                end
            end
            
            if fsave
                if ~exist([fdir 'PFmaps/'],'dir'), mkdir([fdir 'PFmaps/']); end
                if Nspk/nPlot <= 1
                    if Nepochs == 1
                        fname_fig = [fdir 'PFmaps/' mouseid '_' expname '_ref_' files(1,:) '_PFmaps'];
                    else
                        fname_fig = [fdir 'PFmaps/' mouseid '_' expname '_ref_' files(1,:) '_PFmaps_' num2str(e) 'of' num2str(Nepochs) 'ep' ];
                    end
                else
                    if Nepochs == 1
                        fname_fig = [fdir 'PFmaps/' mouseid '_' expname '_ref_' files(1,:) '_PFmaps_' num2str(ii+1)];
                    else
                        fname_fig = [fdir 'PFmaps/' mouseid '_' expname '_ref_' files(1,:) '_PFmaps_' num2str(ii+1) '_' num2str(e) 'of' num2str(Nepochs) 'ep' ];
                    end
                end
                savefig( fh, fname_fig );
                saveas( fh, fname_fig, 'png' );
                if fclose, close( fh ); end
            end
        end 
    end
end