function plotPF_2d(hist, asd, PFdata, activeData, fclose, fsave, sdir, fname_pref)
    if nargin < 8, fsave = false; fname_pref = ''; sdir = ''; end
    if nargin < 7, fsave = false; end
    if nargin < 6, fsave = false; end
    if nargin < 5, fclose = false; fsave = false; end

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
        fname = [sdir 'hist_SI_bitspersec/' fname_pref '_PFmaps_PCs'];
        if ~isempty(asd)
            plotpfMaps_2d(activeData, hist.rMap, hist.rMap_sm, asd.rMap, hist.SIsec.pcIdx, 'PC', fname, fclose);
        else
            plotpfMaps_2d(activeData, hist.rMap, hist.rMap_sm, [], hist.SIsec.pcIdx, 'PC', fname, fclose);
        end
        
        fname = [sdir 'hist_SI_bitspersec/' fname_pref '_centroid_PCs'];
        plot_pfChar( Npcs, hist.SIsec.normpfMap_sm, hist.centroid(hist.SIsec.pcIdx), 'PC', fname, fclose );

        fname = [sdir 'hist_SI_bitspersec/' fname_pref  '_populSummary_PC'];
        plot_populSummary( PFdata.spk_eventrate, PFdata.spk_rate, hist.infoMap(:,1), hist.fieldSize, ...
            hist.centroid, hist.SIsec.pcIdx, hist.SIsec.nonpcIdx, hist.bin_pos, activeData.r, 'bits/s',fname, fclose );
            
        if Nepochs > 1
            fname = [sdir 'hist_SI_bitspersec/' fname_pref '_remapping_hist_SIsec'];
            plotRemapping( hist.SIsec.normpfMap_sm, fname, fclose );
        end
     end
    
    % plots for NON-place cells 
    if ~isempty(hist.SIsec.nonpcIdx)
        Nnonpcs = numel(hist.SIsec.nonpcIdx);
        fname = [sdir 'hist_SI_bitspersec/' fname_pref '_pfMaps_nonPCs'];
        if ~isempty(asd)
            plotpfMaps_2d(activeData, hist.rMap, hist.rMap_sm, asd.rMap, hist.SIsec.nonpcIdx, 'nonPC', fname, fclose);
        else
            plotpfMaps_2d(activeData, hist.rMap, hist.rMap_sm, [], hist.SIsec.nonpcIdx, 'nonPC', fname, fclose);
        end
        
        fname = [sdir 'hist_SI_bitspersec/' fname_pref '_centroid_nonPCs'];
        plot_pfChar( Nnonpcs, hist.normrMap_sm(:,:,hist.SIsec.nonpcIdx), hist.centroid(hist.SIsec.nonpcIdx), 'nonPC', fname, fclose );

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
            plotpfMaps_2d(activeData, hist.rMap, hist.rMap_sm, asd.rMap, hist.SIspk.pcIdx, 'PC', fname, fclose);
        else
            plotpfMaps_2d(activeData, hist.rMap, hist.rMap_sm, [], hist.SIspk.pcIdx, 'PC', fname, fclose);
        end
        
        fname = [sdir 'hist_SI_bitsperspk/' fname_pref '_centroid_PCs'];
        plot_pfChar( Npcs, hist.SIspk.normpfMap_sm, hist.centroid(hist.SIspk.pcIdx), 'PC', fname, fclose );

        fname = [sdir 'hist_SI_bitsperspk/' fname_pref  '_populSummary_PC'];
        plot_populSummary( PFdata.spk_eventrate, PFdata.spk_rate, hist.infoMap(:,2), hist.fieldSize, ...
            hist.centroid, hist.SIspk.pcIdx, hist.SIspk.nonpcIdx, hist.bin_pos, activeData.r, 'bits/event', fname, fclose );
            
        if Nepochs > 1
            fname = [sdir 'hist_SI_bitsperspk/' fname_pref '_remapping_hist_SIsec'];
            plotRemapping( hist.SIspk.normpfMap_sm, fname, fclose );
        end
    end
    
    % plots for NON-place cells 
    if ~isempty(hist.SIspk.nonpcIdx)
        Nnonpcs = numel(hist.SIspk.nonpcIdx);
        fname = [sdir 'hist_SI_bitsperspk/' fname_pref '_pfMaps_nonPCs'];
        if ~isempty(asd)
            plotpfMaps_2d( activeData, hist.rMap, hist.rMap_sm, asd.rMap, hist.SIspk.nonpcIdx, 'nonPC', fname, fclose );
        else
            plotpfMaps_2d( activeData, hist.rMap, hist.rMap_sm, [], hist.SIspk.nonpcIdx, 'nonPC', fname, fclose);
        end
        
        fname = [sdir 'hist_SI_bitsperspk/' fname_pref '_centroid_nonPCs'];
        plot_pfChar( Nnonpcs, hist.normrMap_sm(:,:,hist.SIspk.nonpcIdx), hist.centroid(hist.SIspk.nonpcIdx), 'nonPC', fname, fclose );

    end
    
    
%% subfunctions
    function plot_pfChar( Npcs, normpfMap_sm, centroid_pc, title_str, fname, fclose )
        [nRow, nCol] = getnRownCol(Npcs);
        nPlot = nRow*nCol;
        Nfig = ceil(Npcs/nPlot)-1;
        if Nfig<0, Nfig = 0; end

        for ii=0:Nfig
            fh = figure;
            ha = tight_subplot(nRow,nCol,[.02 .015],[.01 .05],[.01 .01]);
            for jj=0:nPlot-1
                if (ii*nPlot+jj+1) <= Npcs
                    axes(ha(jj+1));
                    map = viridisMap; 
                    imagesc(normpfMap_sm(:,:,ii*nPlot+jj+1)); colormap(map);
                    yticks([]); xticks([]); 
                    title([title_str ' ' num2str(ii*nPlot+jj+1)],'fontsize',12);
                    if ~isempty(centroid_pc)
                        hold on
                        [row, col] = ind2sub(Nbins,centroid_pc(ii*nPlot+jj+1));
                        plot(row,col,'ko','markerfacecolor','k','markersize',8); 
                        hold off; 
                    end
                end
            end
            if ~isempty(fname)
                if Npcs/nPlot <= 1
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

    function plot_populSummary(spk_eventrate_all, spk_rate_all, infoMap_all, fieldSize_all, centroid_all, ...
                                pcIdx, nonpcIdx, bin_pos, activer, info_str, fname, fclose )
        % summary of population data
        for e = 1:Nepochs
            fh = figure;
            subplot(241); 
                hold on
                plot(mean(spk_eventrate_all(pcIdx)),0.2,'v','markerfacecolor',[0.4 0.6 1],...
                    'markeredgecolor',[0.4 0.6 1],'markersize',4);
                plot(mean(spk_eventrate_all(nonpcIdx)),0.2,'v','markerfacecolor',[1 0.6 0.6],...
                    'markeredgecolor',[1 0.6 0.6],'markersize',4);
                h = histogram(spk_eventrate_all(pcIdx),'Normalization','probability','DisplayStyle','bar',...
                    'BinMethod','sqrt','FaceColor',[0.4 0.6 1],'LineStyle','none',...
                    'FaceAlpha',0.7);
                histogram(spk_eventrate_all(nonpcIdx),'Normalization','probability','DisplayStyle','bar',...
                    'FaceColor',[1 0.6 0.6],'LineStyle','none',...
                    'FaceAlpha',0.7, 'BinWidth', h.BinWidth);
                legend('pc','nonpc'); legend('boxoff');
                hold off
                xlabel('Event rate (Hz)'); 
                ylabel('Prop. of cells');
            subplot(242); 
                hold on
                plot(mean(spk_rate_all(pcIdx)),0.3,'v','markerfacecolor',[0.4 0.6 1],...
                    'markeredgecolor',[0.4 0.6 1],'markersize',4);
                plot(mean(spk_rate_all(nonpcIdx)),0.3,'v','markerfacecolor',[1 0.6 0.6],...
                    'markeredgecolor',[1 0.6 0.6],'markersize',4);
                h = histogram(spk_rate_all(pcIdx),'Normalization','probability','DisplayStyle','bar',...
                    'BinMethod','sqrt','FaceColor',[0.4 0.6 1],'LineStyle','none',...
                    'FaceAlpha',0.7);
                histogram(spk_rate_all(nonpcIdx),'Normalization','probability','DisplayStyle','bar',...
                    'FaceColor',[1 0.6 0.6],'LineStyle','none',...
                    'FaceAlpha',0.7, 'BinWidth', h.BinWidth);
                % legend('pc','nonpc'); legend('boxoff');
                hold off
                xlabel('Ampl rate (Hz)'); 
                ylabel('Prop. of cells');
            subplot(243); 
                hold on
                plot(mean(infoMap_all(pcIdx)),0.3,'v','markerfacecolor',[0.4 0.6 1],...
                    'markeredgecolor',[0.4 0.6 1],'markersize',4);
                plot(mean(infoMap_all(nonpcIdx)),0.3,'v','markerfacecolor',[1 0.6 0.6],...
                    'markeredgecolor',[1 0.6 0.6],'markersize',4);
                h = histogram(infoMap_all(pcIdx),'Normalization','probability','DisplayStyle','bar',...
                    'BinMethod','sqrt','FaceColor',[0.4 0.6 1],'LineStyle','none',...
                    'FaceAlpha',0.7);
                histogram(infoMap_all(nonpcIdx),'Normalization','probability','DisplayStyle','bar',...
                    'FaceColor',[1 0.6 0.6],'LineStyle','none',...
                    'FaceAlpha',0.7, 'BinWidth', h.BinWidth);
                % legend('pc','nonpc'); legend('boxoff');
                hold off
                xlabel(['Information (' info_str ')']); % xlabel('Information');
                ylabel('Prop. of cells');
            subplot(244); histogram(fieldSize_all(pcIdx,1,e),'Normalization','probability','FaceColor',[0.5 0.5 0.5],...
                'LineStyle','none','BinMethod','sqrt');
                hold on;
                mean_fieldsize = mean(fieldSize_all(pcIdx,1,e));
                plot(mean_fieldsize,0.2,'kv','markerfacecolor','k','markersize',6); hold off;
                xlabel('Field size (cm^2)');
                ylabel('Fraction of Cells');
            subplot(245); histogram(fieldSize_all(pcIdx,2,e),'Normalization','probability','FaceColor',[0.5 0.5 0.5],...
                'LineStyle','none','BinMethod','sqrt');
                hold on;
                mean_fieldsize = mean(fieldSize_all(pcIdx,2,e));
                plot(mean_fieldsize,0.2,'kv','markerfacecolor','k','markersize',6); hold off;
                xlabel('Field size (% environment)');
                ylabel('Fraction of Cells');
            subplot(246); histogram(centroid_all(pcIdx,e),'Normalization','probability','FaceColor',[0.5 0.5 0.5],...
                'LineStyle','none','BinMethod','sqrt');
                hold on;
                mean_pfLoc = mean(centroid_all(pcIdx,e));
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

            if ~isempty(fname)
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

    function plotRemapping(normpfMap_sm, fname, fclose)
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
