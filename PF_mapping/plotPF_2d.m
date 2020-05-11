function plotPF_2d(hist, asd, activeData, fclose, fsave, sdir, fname_pref)
    if nargin < 7, fsave = false; fname_pref = ''; sdir = ''; end
    if nargin < 6, fsave = false; end
    if nargin < 5, fsave = false; end
    if nargin < 4, fclose = false; fsave = false; end

    Nepochs = size(hist.occMap,3);
    Nbins = size(hist.occMap(:,:,1));

%% SI (bits/sec)
    % histogram estimation
    % plots for PLACE cells 
    if fsave
        if ~exist([sdir 'hist_SI_bitspersec/'],'dir'), mkdir([sdir 'hist_SI_bitspersec/']); end
    end
    if ~isempty(hist.SIsec.pcIdx)
        Npcs = numel(hist.SIsec.pcIdx);
        fname = [sdir 'hist_SI_bitspersec/' fname_pref '_pfMaps_PCs'];
        if ~isempty(asd)
            pfMaps_2d(activeData, hist.SIsec.pfMap, hist.SIsec.pfMap_sm, asd.rMap(:,:,hist.SIsec.pcIdx,:), hist.SIsec.pcIdx, 'PC');
        else
            pfMaps_2d(activeData, hist.SIsec.pfMap, hist.SIsec.pfMap_sm, [], hist.SIsec.pcIdx, 'PC');
        end
        
        fname = [sdir 'hist_SI_bitspersec/' fname_pref '_centroid_PCs'];
        plot_pfChar( Npcs, hist.SIsec.normpfMap_sm, hist.SIsec.centroid, 'PC', fsave, fname, fclose );

        fname = [sdir 'hist_SI_bitspersec/' fname_pref  '_populSummary_PC'];
        plot_populSummary( hist.SIsec.spkPeak, hist.SIsec.spkMean, hist.SIsec.infoMap,...
            hist.SIsec.fieldSize, hist.SIsec.centroid, hist.bin_pos, activeData.r, fsave, fname, fclose );
            
        if Nepochs > 1
            fname = [sdir 'hist_SI_bitspersec/' fname_pref '_remapping_hist_SIsec'];
            plotRemapping( hist.SIsec.normpfMap_sm, fname );
        end
     end
    
    % spike raster plots for NON-place cells 
    if ~isempty(hist.SIsec.nonpcIdx)
        Nnonpcs = numel(hist.SIsec.nonpcIdx);
        fname = [sdir 'hist_SI_bitspersec/' fname_pref '_pfMaps_nonPCs'];
        if ~isempty(asd)
            pfMaps_2d(activeData, hist.rMap(:,:,hist.SIsec.nonpcIdx), hist.rMap_sm(:,:,hist.SIsec.nonpcIdx), ...
                        asd.rMap(:,:,hist.SIsec.nonpcIdx), hist.SIsec.nonpcIdx, 'nonPC');
        else
            pfMaps_2d(activeData, hist.rMap(:,:,hist.SIsec.nonpcIdx), hist.rMap_sm(:,:,hist.SIsec.nonpcIdx), ...
                        [], hist.SIsec.nonpcIdx, 'nonPC');
        end
        
        fname = [sdir 'hist_SI_bitspersec/' fname_pref '_centroid_nonPCs'];
        plot_pfChar( Nnonpcs, hist.normrMap_sm(:,:,hist.SIsec.nonpcIdx), hist.centroid(hist.SIsec.nonpcIdx), 'nonPC', fsave, fname, fclose );

    end
    
%% SI (bits/spike)
    % histogram estimation
    % plots for PLACE cells 
    if fsave
        if ~exist([sdir 'hist_SI_bitsperspk/'],'dir'), mkdir([sdir 'hist_SI_bitsperspk/']); end
    end
    if ~isempty(hist.SIspk.pcIdx)
        Npcs = numel(hist.SIspk.pcIdx);
        fname = [sdir 'hist_SI_bitsperspk/' fname_pref '_pfMaps_PCs'];
        if ~isempty(asd)
            pfMaps_2d(activeData, hist.SIspk.pfMap, hist.SIspk.pfMap_sm, asd.rMap(:,:,hist.SIspk.pcIdx,:), hist.SIspk.pcIdx, 'PC');
        else
            pfMaps_2d(activeData, hist.SIspk.pfMap, hist.SIspk.pfMap_sm, [], hist.SIspk.pcIdx, 'PC');
        end
        
        fname = [sdir 'hist_SI_bitsperspk/' fname_pref '_centroid_PCs'];
        plot_pfChar( Npcs, hist.SIspk.normpfMap_sm, hist.SIspk.centroid, 'PC', fsave, fname, fclose );

        fname = [sdir 'hist_SI_bitsperspk/' fname_pref  '_populSummary_PC'];
        plot_populSummary( hist.SIspk.spkPeak, hist.SIspk.spkMean, hist.SIspk.infoMap,...
            hist.SIspk.fieldSize, hist.SIspk.centroid, hist.bin_pos, activeData.r, fsave, fname, fclose );
            
        if Nepochs > 1
            fname = [sdir 'hist_SI_bitsperspk/' fname_pref '_remapping_hist_SIsec'];
            plotRemapping( hist.SIspk.normpfMap_sm, fname );
        end
    end
    
    % spike raster plots for NON-place cells 
    if ~isempty(hist.SIspk.nonpcIdx)
        Nnonpcs = numel(hist.SIspk.nonpcIdx);
        fname = [sdir 'hist_SI_bitsperspk/' fname_pref '_pfMaps_nonPCs'];
        if ~isempty(asd)
            pfMaps_2d(activeData, hist.rMap(:,:,hist.SIspk.nonpcIdx), hist.rMap_sm(:,:,hist.SIspk.nonpcIdx), ...
                        asd.rMap(:,:,hist.SIspk.nonpcIdx), hist.SIspk.nonpcIdx, 'nonPC');
        else
            pfMaps_2d(activeData, hist.rMap(:,:,hist.SIspk.nonpcIdx), hist.rMap_sm(:,:,hist.SIspk.nonpcIdx), ...
                        [], hist.SIspk.nonpcIdx, 'nonPC');
        end
        
        fname = [sdir 'hist_SI_bitsperspk/' fname_pref '_centroid_nonPCs'];
        plot_pfChar( Nnonpcs, hist.normrMap_sm(:,:,hist.SIspk.nonpcIdx), hist.centroid(hist.SIspk.nonpcIdx), 'nonPC', fsave, fname, fclose );

    end
    
    
%% subfunctions
    function pfMaps_2d(activeData, pfMap_h, pfMap_sm_h, pfMap_a, pcIdx, cell_str)
        Nepochs = size( pfMap_h, 4 );
        Ncells = length(pcIdx);
        nRow = 4;
        if isempty(asd)
            nCol = 3;
        else
            nCol = 4;
        end
        
        for e = 1:Nepochs
            for ii=0:ceil(Ncells/nRow)-1 
                fh = figure; 
                ha = tight_subplot(nRow,nCol,[.01 .005],[.01 .07],[.01 .01]);
                for jj=0:nRow-1
                    if (ii*nRow+jj+1) <= Ncells
                        axes(ha(jj*nCol+1));
                        z = activeData.spikes(pcIdx(ii*nRow+jj+1),:);
                        hold on; axis off;
                        plot(activeData.x,-activeData.y,'Color',[0.6 0.6 0.6],'LineWidth',1.2); 
                        ind = find(z>0);
                        x = activeData.x(ind);
                        y = activeData.y(ind);
                        spikes = z(ind);
                        [spikes_sorted,sort_ind] = sort(spikes);
                        scatter(x(sort_ind),-y(sort_ind),20,spikes_sorted,'filled');
                        title_str = sprintf('%s %d (%0.2g)', cell_str, (ii*nRow+jj+1), max(spikes));
                        title(title_str,'fontsize',15);
                        hold off

                        axes(ha(jj*nCol+2));
                        imagesc(squeeze(pfMap_h(:,:,ii*nRow+jj+1,e))');
                        axis off; colorbar; % caxis([0 0.06]);
                        if Nepochs >1 
                            title(['Epoch ',num2str(e)],'fontsize',15);
                        end
                        axes(ha(jj*nCol+3)); 
                        imagesc(squeeze(pfMap_sm_h(:,:,ii*nRow+jj+1,e))');
                        axis off; colorbar; % caxis([0 0.005]);
                        if ~isempty(asd)
                            axes(ha(jj*nRow+4));
                            imagesc(squeeze(pfMap_a(:,:,ii*nRow+jj+1,e))');
                            axis off; colorbar; % caxis([0 0.003]);
                        end
                    end
                end

                if fsave
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

    function plot_pfChar( Ncells, normpfMap_sm, centroid, title_str, fsave, fname, fclose )
        [nRow, nCol] = getnRownCol(Ncells);
        nPlot = nRow*nCol;
        Nfig = ceil(Ncells/nPlot)-1;
        if Nfig<0, Nfig = 0; end

        for ii=0:Nfig
            fh = figure;
            ha = tight_subplot(nRow,nCol,[.02 .015],[.01 .05],[.01 .01]);
            for jj=0:nPlot-1
                if (ii*nPlot+jj+1) <= Ncells
                    axes(ha(jj+1));
                    map = viridisMap; 
                    imagesc(normpfMap_sm(:,:,ii*nPlot+jj+1)); colormap(map);
                    yticks([]); xticks([]); 
                    title([title_str ' ' num2str(ii*nPlot+jj+1)],'fontsize',12);
                    if ~isempty(centroid)
                        hold on
                        [row, col] = ind2sub(Nbins,centroid(ii*nPlot+jj+1));
                        plot(row,col,'ko','markerfacecolor','k','markersize',8); 
                        hold off; 
                    end
                end
            end
            if fsave
                if Ncells/nPlot <= 1
                    fname_fig = fname;
                else
                    fname_fig = [fname '_' num2str(ii+1)];
                end
                savefig( fh, fname_fig );
                saveas( fh, fname_fig, 'png' );
                if fclose, close( fh ); end
            end
        end 
    end

    function plot_populSummary(spkPeak, spkMean, infoMap, fieldSize, centroid, bin_pos, activer, fsave, fname, fclose )
        % summary of population data
        for e = 1:Nepochs
            fh = figure;
            subplot(241); histogram(spkPeak,'Normalization','probability','FaceColor',[0.5 0.5 0.5],...
                'LineStyle','none','BinMethod','fd'); 
                hold on;
                plot(mean(spkPeak),0.2,'kv','markerfacecolor','k','markersize',6); hold off;
                xlabel('PSD peak'); 
                ylabel('Fraction of Cells');
            subplot(242); histogram(spkMean,'Normalization','probability','FaceColor',[0.5 0.5 0.5],...
                'LineStyle','none','BinMethod','fd'); 
                hold on;
                plot(mean(spkMean),0.3,'kv','markerfacecolor','k','markersize',6); hold off;
                xlabel('PSD mean'); 
                ylabel('Fraction of Cells');
            subplot(243); histogram(infoMap,'Normalization','probability','FaceColor',[0.5 0.5 0.5],...
                'LineStyle','none','BinMethod','fd'); 
                hold on;
                plot(mean(infoMap),0.3,'kv','markerfacecolor','k','markersize',6); hold off;
                xlabel('Info content');
                ylabel('Fraction of Cells');
            subplot(244); histogram(fieldSize(:,1,e),'Normalization','probability','FaceColor',[0.5 0.5 0.5],...
                'LineStyle','none','BinMethod','sqrt');
                hold on;
                mean_fieldsize = mean(fieldSize(:,1,e));
                plot(mean_fieldsize,0.2,'kv','markerfacecolor','k','markersize',6); hold off;
                xlabel('Field size (cm^2)');
                ylabel('Fraction of Cells');
            subplot(245); histogram(fieldSize(:,2,e),'Normalization','probability','FaceColor',[0.5 0.5 0.5],...
                'LineStyle','none','BinMethod','sqrt');
                hold on;
                mean_fieldsize = mean(fieldSize(:,2,e));
                plot(mean_fieldsize,0.2,'kv','markerfacecolor','k','markersize',6); hold off;
                xlabel('Field size (% environment)');
                ylabel('Fraction of Cells');
            subplot(246); histogram(centroid(:,e),'Normalization','probability','FaceColor',[0.5 0.5 0.5],...
                'LineStyle','none','BinMethod','sqrt');
                hold on;
                mean_pfLoc = mean(centroid(:,e));
                plot(mean_pfLoc,0.5,'kv','markerfacecolor','k','markersize',6); hold off;
                xlabel('Position index');
                ylabel('No. of PCs');
            subplot(247); histogram(bin_pos(:,e),'Normalization','probability','FaceColor',[0.5 0.5 0.5],...
                'LineStyle','none','BinMethod','fd');
                hold on;
                plot(mean(bin_pos(:,e)),0.04,'kv','markerfacecolor','k','markersize',6); hold off;
                xlabel('Position index');
                ylabel('Fraction of Time');
            subplot(248); histogram(activer(:,e),'Normalization','probability','FaceColor',[0.5 0.5 0.5],...
                'LineStyle','none','BinMethod','fd');
                hold on;
                plot(mean(activer(:,e)),0.04,'kv','markerfacecolor','k','markersize',6); hold off;
                xlabel('Radius (mm)');
                ylabel('Fraction of Time');

            if fsave
                if Nepochs == 1
                    fname_fig = fname;
                else
                    fname_fig = [fname '_' num2str(e) 'of' num2str(Nepochs) 'ep'];
                end
                savefig( fh, fname_fig );
                saveas( fh, fname_fig, 'png' );
                if fclose, close( fh ); end
            end
        end
    end

    function plotRemapping(normpfMap_sm, fname)
        % remapping within a session
        Nepochs = size(pfMap_sm,4);
        Ncells = size(pfMap_sm,3);
        nRow = 4;
        nCol = Nepochs;
        
        for e = 1:Nepochs
            for ii=0:(Ncells/nRow)-1 
                fh = figure; 
                ha = tight_subplot(nRow,nCol,[.01 .005],[.01 .07],[.01 .01]);
                for jj=0:nCol
                    if (ii*nRow+jj) <= Ncells
                        axes(ha(jj*nRow+e));
                        imagesc(normpfMap_sm(:,:,ii*nRow+jj,e));
                        title(['ep ' e]);
                        if e == 1
                            ylabel(['PC ',num2str(ii*nRow+jj+1)],'fontsize',15);
                        end
                    end
                end
            end
            if fsave
                savefig( fh, fname );
                saveas( fh, fname, 'png' );
                if fclose, close( fh ); end
            end
        end
    end

end
