function plotPF_1d(hist, asd, PFdata, fclose, fsave, sdir, fname_pref)
    if nargin < 7, fsave = false; fname_pref = '_'; sdir = '_'; end
    if nargin < 6, fsave = false; end
    if nargin < 5, fsave = false; end
    if nargin < 4, fclose = false; fsave = false; end

    plotNon_norm = true;
    ytick_files = PFdata.ytick_files;
    Nbins = size(PFdata.occMap,2);
    if iscell(hist.rateMap)
        Nepochs = length(hist.rateMap);
    else
        Nepochs = 1;
    end
    
%% SI (bits/sec)
    % histogram estimation
    % plots for PLACE cells 
    if fsave
        if Nepochs == 1
            if ~exist([sdir 'hist_SI_bitspersec/'],'dir'), mkdir([sdir 'hist_SI_bitspersec/']); end
        else
            for e = 1:Nepochs
                if ~exist([sdir 'hist_SI_bitspersec/epoch' num2str(e)],'dir'), mkdir([sdir 'hist_SI_bitspersec/epoch' num2str(e)]); end
            end
        end
    end
    if ~isempty(hist.SIsec.pcIdx)
        if Nepochs == 1
            Npcs = numel(hist.SIsec.pcIdx);
            fname = [sdir 'hist_SI_bitspersec/' fname_pref '_normspkRaster_PCs'];
            plotRaster( Npcs, PFdata.spkRaster(), ytick_files, 'PC', fsave, fname, fclose );

            fname = [sdir 'hist_SI_bitspersec/' fname_pref '_rateplot_PCs'];
            plot_pfLocSize( Npcs, hist.SIsec.normpfMap_sm, hist.SIsec.pfLoc, hist.SIsec.pfSize, 'PC', fsave, fname, fclose )

            if plotNon_norm
                fname = [sdir 'hist_SI_bitspersec/' fname_pref '_spkRaster_PCs'];
                plotRaster( Npcs, hist.SIsec.spkRaster_pc, ytick_files, 'PC', fsave, fname, fclose );
            end

            fname = [sdir 'hist_SI_bitspersec/' fname_pref  '_populSummary'];
            plot_populSummary( hist.SIsec.sort_normpfMap_sm, hist.SIsec.sort_pfMap_sm, hist.SIsec.spkPeak, hist.SIsec.spkMean, hist.SIsec.infoMap,...
                hist.SIsec.pfSize, hist.SIsec.pfLoc, PFdata.bin_phi, fsave, fname, fclose );
        end
    end
    
    % spike raster plots for NON-place cells 
    if ~isempty(hist.SIsec.nonpcIdx)
        Nnonpcs = numel(hist.SIsec.nonpcIdx);
        fname = [sdir 'hist_SI_bitspersec/' fname_pref '_normspkRaster_nonPCs'];
        plotRaster( Nnonpcs, hist.SIsec.normspkRaster_nonpc, ytick_files, 'Non-PC', fsave, fname, fclose );
        
        if plotNon_norm
            fname = [sdir 'hist_SI_bitspersec/' fname_pref '_spkRaster_nonPCs'];
            plotRaster( Nnonpcs, hist.SIsec.spkRaster_nonpc, ytick_files, 'Non-PC', fsave, fname, fclose );
        end
        
        fname = [sdir 'hist_SI_bitspersec/' fname_pref '_rateplot_nonPCs'];
        plot_pfLocSize( Nnonpcs, hist.normrMap_sm(hist.SIsec.nonpcIdx,:), [], [], 'Non-PC', fsave, fname, fclose )
    end
    

    % asd
    % plots for PLACE cells 
    if ~isempty(asd)
        if fsave
            if ~exist([sdir 'asd_SI_bitspersec/'],'dir'), mkdir([sdir 'asd_SI_bitspersec/']); end
        end
        if ~isempty(asd.SIsec.pcIdx)
            Npcs = numel(asd.SIsec.pcIdx);
            fname = [sdir 'asd_SI_bitspersec/' fname_pref '_normspkRaster_PCs'];
            plotRaster( Npcs, asd.SIsec.normspkRaster_pc, ytick_files, 'PC', fsave, fname, fclose );

            fname = [sdir 'asd_SI_bitspersec/' fname_pref '_rateplot_PCs'];
            plot_pfLocSize( Npcs, asd.SIsec.normpfMap, asd.SIsec.pfLoc, asd.SIsec.pfSize, 'PC', fsave, fname, fclose )

            if plotNon_norm
                fname = [sdir 'asd_SI_bitspersec/' fname_pref '_spkRaster_PCs'];
                plotRaster( Npcs, asd.SIsec.spkRaster_pc, ytick_files, 'PC', fsave, fname, fclose );
            end

            fname = [sdir 'asd_SI_bitspersec/' fname_pref  '_populSummary'];
            plot_populSummary( asd.SIsec.sort_normpfMap, asd.SIsec.sort_pfMap, asd.SIsec.spkPeak, asd.SIsec.spkMean, asd.SIsec.infoMap,...
                asd.SIsec.pfSize, asd.SIsec.pfLoc, PFdata.bin_phi, fsave, fname, fclose );

            if Nepochs > 1
                fname = [sdir 'asd_SI_bitspersec/' fname_pref '_remapping_asd_SIsec'];
                plotRemapping( asd.SIsec.normpfMap, asd.SIsec.sortIdx, fname );
            end
        end

        % spike raster plots for NON-place cells 
        if ~isempty(asd.SIsec.nonpcIdx)
            Nnonpcs = numel(asd.SIsec.nonpcIdx);
            fname = [sdir 'asd_SI_bitspersec/' fname_pref '_normspkRaster_nonPCs'];
            plotRaster( Nnonpcs, asd.SIsec.normspkRaster_nonpc, ytick_files, 'Non-PC', fsave, fname, fclose );

            if plotNon_norm
                fname = [sdir 'asd_SI_bitspersec/' fname_pref '_spkRaster_nonPCs'];
                plotRaster( Nnonpcs, asd.SIsec.spkRaster_nonpc, ytick_files, 'Non-PC', fsave, fname, fclose );
            end

            fname = [sdir 'asd_SI_bitspersec/' fname_pref '_rateplot_nonPCs'];
            plot_pfLocSize( Nnonpcs, asd.normrMap(asd.SIsec.nonpcIdx,:), [], [], 'Non-PC', fsave, fname, fclose )
        end
    end
    
%% SI (bits/spike)

    
%% PLOTTING FUNCTIONS  
% spike rasters
function plotRaster( Ncells, spkRaster, ytick_files, title_str, fsave, fname, fclose )
    [nRow, nCol] = getnRownCol(Ncells);
    nPlot = nRow*nCol;
    Nfig = ceil(Ncells/nPlot)-1;
    if Nfig<0, Nfig = 0; end

    for ii=0:Nfig
        fh = figure;
        ha = tight_subplot(nRow,nCol,[.05 0.008],[.08 .05],[.06 .02]);
        for jj=0:nPlot-1
            if (ii*nPlot+jj+1) <= Ncells
                axes(ha(jj+1));
                map = viridisMap; 
                imagesc(spkRaster{ii*nPlot+jj+1}); colormap(map);
                % only put ylabels for 1st column plots
                if ii*nPlot+jj+1 < 8
                    if ii*nPlot+jj+1 == 1
                        yticks(ytick_files); yticklabels(ytick_files); 
                        ylabel('Lap #');
                    else
                        yticks([]);
                    end
                else
                    if mod(ii*nPlot+jj+1,nCol) == 1
                        yticks(ytick_files); yticklabels(ytick_files); 
                        ylabel('Lap #');
                    else
                        yticks([]);
                    end
                end
                % only put xlabels for last row plots
                if ii*nPlot+jj+1 > (ii+1)*nPlot - nCol && ii*nPlot+jj+1 < (ii+1)*nPlot + 1
                    xticks([1 Nbins]); xticklabels([1 100]); xlabel('Pos (cm)');
                else
                    xticks([]); 
                end
                title([title_str ' ' num2str(ii*nPlot+jj+1)],'fontsize',12);
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

% pf location and size
function plot_pfLocSize( Ncells, normpfMap_sm, pfLoc, pfBins, title_str, fsave, fname, fclose )
    [nRow, nCol] = getnRownCol(Ncells);
    nPlot = nRow*nCol;
    Nfig = ceil(Ncells/nPlot)-1;
    if Nfig<0, Nfig = 0; end

    for ii=0:Nfig
        fh = figure;
        ha = tight_subplot(nRow,nCol,[.05 .005],[.08 .05],[.06 .02]);
        for jj=0:nPlot-1
            if (ii*nPlot+jj+1) <= Ncells
                axes(ha(jj+1));
                map = viridisMap; 
                plot(normpfMap_sm(ii*nPlot+jj+1,:)); colormap(map);
                yticks([]); xticks([]); 
                title([title_str ' ' num2str(ii*nPlot+jj+1)],'fontsize',12);
                if ~isempty(pfLoc)
                    hold on
                    plot(pfLoc(ii*nPlot+jj+1),1,'kv','markerfacecolor','k','markersize',4); 
                    for i = 1:length(pfBins{ii*nPlot+jj+1})
                        plot(pfBins{ii*nPlot+jj+1}(i),0.5,'r*','markerfacecolor','r','markersize',4)
                    end
                    hold off; 
                end
                % only put ylabels for 1st column plots
                if ii*nPlot+jj+1 < 8
                    if ii*nPlot+jj+1 == 1
                        yticks([0 1]); 
                        ylabel('Ampl rate');
                    else
                        yticks([]);
                    end
                else
                    if mod(ii*nPlot+jj+1,nCol) == 1
                        yticks([0 1]);  
                        ylabel('Ampl rate');
                    else
                        yticks([]);
                    end
                end
                % only put xlabels for last row plots
                if ii*nPlot+jj+1 > (ii+1)*nPlot - nCol && ii*nPlot+jj+1 < (ii+1)*nPlot + 1
                    xticks([1 Nbins]); xticklabels([1 100]); xlabel('Pos (cm)');
                else
                    xticks([]); 
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

% Row 1: bin dwell time per trial, total bin dwell time, % of bin dwell
% time cell is active, % of laps cell is active
% Row 2: sorted pf maps, event rate histogram (pcs & nonpcs), spike
% amplitude rate histogram (pcs & nonpcs), histogram of info (pcs & nonpcs)
% Row 3: sorted norm pf maps, histogram of fieldsize, histogram of field
% loc, frac of pf dwell time cell is active
function plot_populSummary(bintime_trials, bintime, bin_activet, activetrials,...
        sort_pfMap_sm, sort_normpfMap_sm, spk_eventrate, spk_rate, infoMap, ...
        fieldSize, pfLoc, pf_activet, pcIdx, nonpcIdx, info_str, fsave, fname, fclose )
    % summary of population data
    Ncells = length(activetrials);
    Npcs = length(pcIdx);
    fh = figure('Position',[680 678 800 650]);
    subplot(341); 
        imagesc(bintime_trials); 
        xticks([1 Nbins]); xticklabels([1 100]); xlabel('Position (cm)');
        if size(bintime_trials,1)>1 
            yticks([1 size(bintime_trials,1)]); yticklabels([1 size(bintime_trials,1)]);
        else 
            yticks(1); yticklabels(1); 
        end
        ylabel('Lap #'); 
        title('Dwell time (s)'); 
        c = colorbar; c.Ticks = [0 round(max(max(bintime_trials)),1)];
    subplot(342); bar(bintime,'FaceColor',[0.5 0.5 0.5],'LineStyle','none');
        xticks([1 Nbins]); xticklabels([1 100]); xlabel('Position (cm)');
        ylabel('Total dwell time (s)');
        box off;
    subplot(343); 
        imagesc(bin_activet);
        xticks([1 Nbins]); xticklabels([1 100]); xlabel('Position (cm)');
        if Ncells>1 
            yticks([1 Ncells]); yticklabels([1 Ncells]);
        else 
            yticks(1); yticklabels(1); 
        end
        ylabel('Cell #');
        title('Active dwell time (%)'); 
        c = colorbar; c.Ticks = [0 round(max(max(bin_activet)),1)];
    subplot(344); 
        hold on
        plot(mean(activetrials(pcIdx)),0.3,'v','markerfacecolor',[0.4 0.6 1],...
            'markeredgecolor',[0.4 0.6 1],'markersize',4);
        plot(mean(activetrials(nonpcIdx)),0.3,'v','markerfacecolor',[1 0.6 0.6],...
            'markeredgecolor',[1 0.6 0.6],'markersize',4);
        h = histogram(activetrials(pcIdx),'Normalization','probability','DisplayStyle','bar',...
            'BinMethod','sqrt','FaceColor',[0.4 0.6 1],'LineStyle','none',...
            'FaceAlpha',0.7);
        histogram(activetrials(nonpcIdx),'Normalization','probability','DisplayStyle','bar',...
            'FaceColor',[1 0.6 0.6],'LineStyle','none',...
            'FaceAlpha',0.7, 'BinWidth', h.BinWidth);
        legend('pc','nonpc'); legend('boxoff');
        hold off;
        ylabel('Prop. of cells'); 
        xlabel('Active laps (%)');
    subplot(345); 
        map = viridisMap;
        imagesc(sort_pfMap_sm); colormap(map);
        xticks([1 Nbins]); xticklabels([1 100]); xlabel('Position (cm)');
        if Ncells>1 
            yticks([1 Npcs]); yticklabels([1 Npcs ]);
        else 
            yticks(1); yticklabels(1); 
        end
        ylabel('PC #'); 
        title('PF map'); 
        c = colorbar; c.Ticks = [0 round(max(max(sort_pfMap_sm)),1)];
    subplot(349); 
        imagesc(sort_normpfMap_sm); 
        xticks([1 Nbins]); xticklabels([1 100]); xlabel('Position (cm)');
        if Ncells>1 
            yticks([1 Npcs]); yticklabels([1 Npcs]);
        else 
            yticks(1); yticklabels(1); 
        end
        ylabel('PC #'); 
        title('Normalised PF map'); 
        c = colorbar; c.Ticks = [0 1];
    subplot(346); 
        hold on
        plot(mean(spk_eventrate(pcIdx)),0.3,'v','markerfacecolor',[0.4 0.6 1],...
            'markeredgecolor',[0.4 0.6 1],'markersize',4);
        plot(mean(spk_eventrate(nonpcIdx)),0.3,'v','markerfacecolor',[1 0.6 0.6],...
            'markeredgecolor',[1 0.6 0.6],'markersize',4);
        h = histogram(spk_eventrate(pcIdx),'Normalization','probability','DisplayStyle','bar',...
            'BinMethod','sqrt','FaceColor',[0.4 0.6 1],'LineStyle','none',...
            'FaceAlpha',0.7);
        histogram(spk_eventrate(nonpcIdx),'Normalization','probability','DisplayStyle','bar',...
            'FaceColor',[1 0.6 0.6],'LineStyle','none',...
            'FaceAlpha',0.7, 'BinWidth', h.BinWidth);
        legend('pc','nonpc'); legend('boxoff');
        hold off
        xlabel('Event rate (Hz)'); 
        ylabel('Prop. of cells');
    subplot(347); 
        hold on
        plot(mean(spk_rate(pcIdx)),0.3,'v','markerfacecolor',[0.4 0.6 1],...
            'markeredgecolor',[0.4 0.6 1],'markersize',4);
        plot(mean(spk_rate(nonpcIdx)),0.3,'v','markerfacecolor',[1 0.6 0.6],...
            'markeredgecolor',[1 0.6 0.6],'markersize',4);
        h = histogram(spk_rate(pcIdx),'Normalization','probability','DisplayStyle','bar',...
            'BinMethod','sqrt','FaceColor',[0.4 0.6 1],'LineStyle','none',...
            'FaceAlpha',0.7);
        histogram(spk_rate(nonpcIdx),'Normalization','probability','DisplayStyle','bar',...
            'FaceColor',[1 0.6 0.6],'LineStyle','none',...
            'FaceAlpha',0.7, 'BinWidth', h.BinWidth);
        legend('pc','nonpc'); legend('boxoff');
        hold off
        xlabel('Ampl rate (Hz)'); 
        ylabel('Prop. of cells');
    subplot(348); 
        hold on
        plot(mean(infoMap(pcIdx)),0.3,'v','markerfacecolor',[0.4 0.6 1],...
            'markeredgecolor',[0.4 0.6 1],'markersize',4);
        plot(mean(infoMap(nonpcIdx)),0.3,'v','markerfacecolor',[1 0.6 0.6],...
            'markeredgecolor',[1 0.6 0.6],'markersize',4);
        h = histogram(infoMap(pcIdx),'Normalization','probability','DisplayStyle','bar',...
            'BinMethod','sqrt','FaceColor',[0.4 0.6 1],'LineStyle','none',...
            'FaceAlpha',0.7);
        histogram(infoMap(nonpcIdx),'Normalization','probability','DisplayStyle','bar',...
            'FaceColor',[1 0.6 0.6],'LineStyle','none',...
            'FaceAlpha',0.7, 'BinWidth', h.BinWidth);
        legend('pc','nonpc'); legend('boxoff');
        hold off
        xlabel(info_str); % xlabel('Information');
        ylabel('Prop. of cells');
    subplot(3,4,10); 
        hold on;
        histogram(fieldSize(pcIdx),'Normalization','probability','FaceColor',[0.4 0.6 1],...
        'LineStyle','none','BinMethod','sqrt','DisplayStyle','bar');
        plot(mean(fieldSize(pcIdx)),0.3,'kv','markerfacecolor',[0.4 0.6 1],...
            'markeredgecolor',[0.4 0.6 1],'markersize',4); hold off;
        hold off;
        xlabel('Field size (cm)');
        ylabel('Prop. of PCs');
    subplot(3,4,11); 
        hold on;
        histogram(pfLoc(pcIdx),'Normalization','probability','FaceColor',[0.4 0.6 1],...
        'LineStyle','none','BinMethod','sqrt','DisplayStyle','bar');
        plot(mean(pfLoc(pcIdx)),0.3,'kv','markerfacecolor',[0.4 0.6 1],...
            'markeredgecolor',[0.4 0.6 1],'markersize',4); hold off;
        hold off;
        xlabel('Position (cm)');
        ylabel('No. of PCs');
    subplot(3,4,12); 
        hold on;
        histogram(pf_activet(pcIdx)*100,'Normalization','probability','FaceColor',[0.4 0.6 1],...
        'LineStyle','none','BinMethod','sqrt','DisplayStyle','bar');
        plot(mean(pf_activet(pcIdx))*100,0.3,'kv','markerfacecolor',[0.4 0.6 1],...
            'markeredgecolor',[0.4 0.6 1],'markersize',4); hold off;
        hold off;
        xlabel('Active dwell time (%)');
        ylabel('Prop. of PCs');

        
    if fsave
        if ~exist(sdir,'dir'), mkdir(sdir); end

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