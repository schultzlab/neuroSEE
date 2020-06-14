function plotPF_1d(hist, asd, PFdata, fclose, fsave, sdir, fname_pref)
    if nargin < 7, fsave = false; fname_pref = '_'; sdir = '_'; end
    if nargin < 6, fsave = false; end
    if nargin < 5, fsave = false; end
    if nargin < 4, fclose = false; fsave = false; end

    plotNon_norm = true;
    Nbins = size(PFdata.occMap,2);
    ytick_files = PFdata.ytick_files;
    
   
%% SI (bits/sec)
    % histogram estimation
    % plots for PLACE cells 
    if fsave
        if ~exist([sdir 'hist_SI_bitspersec/'],'dir'), mkdir([sdir 'hist_SI_bitspersec/']); end
    end
    if ~isempty(hist.SIsec.pcIdx)
        Npcs = numel(hist.SIsec.pcIdx);
        fname = [sdir 'hist_SI_bitspersec/' fname_pref '_normspkRaster_PCs'];
        plotRaster( Npcs, hist.SIsec.normspkRaster_pc, ytick_files, 'PC', fsave, fname, fclose );
        
        fname = [sdir 'hist_SI_bitspersec/' fname_pref '_rateplot_PCs'];
        plot_pfLoc( Npcs, hist.SIsec.normpfMap_sm, hist.SIsec.pfLoc, hist.SIsec.pfSize, 'PC', fsave, fname, fclose )
        
        if plotNon_norm
            fname = [sdir 'hist_SI_bitspersec/' fname_pref '_spkRaster_PCs'];
            plotRaster( Npcs, hist.SIsec.spkRaster_pc, ytick_files, 'PC', fsave, fname, fclose );
        end
        
        fname = [sdir 'hist_SI_bitspersec/' fname_pref  '_populSummary'];
        plot_populSummary( hist.SIsec.sort_normpfMap_sm, hist.SIsec.sort_pfMap_sm, hist.SIsec.spkPeak, hist.SIsec.spkMean, hist.SIsec.infoMap,...
            hist.SIsec.pfSize, hist.SIsec.pfLoc, PFdata.bin_phi, fsave, fname, fclose );
            
        if Nepochs > 1
            fname = [sdir 'hist_SI_bitspersec/' fname_pref '_remapping_hist_SIsec'];
            plotRemapping( hist.SIsec.normpfMap_sm, hist.SIsec.sortIdx, fname );
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
        plot_pfLoc( Nnonpcs, hist.normrMap_sm(hist.SIsec.nonpcIdx,:), [], [], 'Non-PC', fsave, fname, fclose )
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
            plot_pfLoc( Npcs, asd.SIsec.normpfMap, asd.SIsec.pfLoc, asd.SIsec.pfSize, 'PC', fsave, fname, fclose )

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
            plot_pfLoc( Nnonpcs, asd.normrMap(asd.SIsec.nonpcIdx,:), [], [], 'Non-PC', fsave, fname, fclose )
        end
    end
    
%% SI (bits/spike)
    % histogram estimation
    % plots for PLACE cells 
    if fsave
        if ~exist([sdir 'hist_SI_bitsperspk/'],'dir'), mkdir([sdir 'hist_SI_bitsperspk/']); end
    end
    if ~isempty(hist.SIspk.pcIdx)
        Npcs = numel(hist.SIspk.pcIdx);
        fname = [sdir 'hist_SI_bitsperspk/' fname_pref '_normspkRaster_PCs'];
        plotRaster( Npcs, hist.SIspk.normspkRaster_pc, ytick_files, 'PC', fsave, fname, fclose );
        
        fname = [sdir 'hist_SI_bitsperspk/' fname_pref '_rateplot_PCs'];
        plot_pfLoc( Npcs, hist.SIspk.normpfMap_sm, hist.SIspk.pfLoc, hist.SIspk.pfSize, 'PC', fsave, fname, fclose )
        
        if plotNon_norm
            fname = [sdir 'hist_SI_bitsperspk/' fname_pref '_spkRaster_PCs'];
            plotRaster( Npcs, hist.SIspk.spkRaster_pc, ytick_files, 'PC', fsave, fname, fclose );
        end
        
        fname = [sdir 'hist_SI_bitsperspk/' fname_pref  '_populSummary'];
        plot_populSummary( hist.SIspk.sort_normpfMap_sm, hist.SIspk.sort_pfMap_sm, hist.SIspk.spkPeak, hist.SIspk.spkMean, hist.SIspk.infoMap,...
            hist.SIspk.pfSize, hist.SIspk.pfLoc, PFdata.bin_phi, fsave, fname, fclose );
            
        if Nepochs > 1
            fname = [sdir 'hist_SI_bitsperspk/' fname_pref '_remapping_hist_SIspk'];
            plotRemapping( hist.SIspk.normpfMap, hist.SIspk.sortIdx, fname );
        end
    end
    
    % spike raster plots for NON-place cells 
    if ~isempty(hist.SIspk.nonpcIdx)
        Nnonpcs = numel(hist.SIspk.nonpcIdx);
        fname = [sdir 'hist_SI_bitsperspk/' fname_pref '_normspkRaster_nonPCs'];
        plotRaster( Nnonpcs, hist.SIspk.normspkRaster_nonpc, ytick_files, 'Non-PC', fsave, fname, fclose );
        
        if plotNon_norm
            fname = [sdir 'hist_SI_bitsperspk/' fname_pref '_spkRaster_nonPCs'];
            plotRaster( Nnonpcs, hist.SIspk.spkRaster_nonpc, ytick_files, 'Non-PC', fsave, fname, fclose );
        end
        
        fname = [sdir 'hist_SI_bitsperspk/' fname_pref '_rateplot_nonPCs'];
        plot_pfLoc( Nnonpcs, hist.normrMap_sm(hist.SIspk.nonpcIdx,:), [], [], 'Non-PC', fsave, fname, fclose )
    end
    

    % asd
    % plots for PLACE cells 
    if ~isempty(asd)
        if fsave
            if ~exist([sdir 'asd_SI_bitsperspk/'],'dir'), mkdir([sdir 'asd_SI_bitsperspk/']); end
        end
        if ~isempty(asd.SIspk.pcIdx)
            Npcs = numel(asd.SIspk.pcIdx);
            fname = [sdir 'asd_SI_bitsperspk/' fname_pref '_normspkRaster_PCs'];
            plotRaster( Npcs, asd.SIspk.normspkRaster_pc, ytick_files, 'PC', fsave, fname, fclose );

            fname = [sdir 'asd_SI_bitsperspk/' fname_pref '_rateplot_PCs'];
            plot_pfLoc( Npcs, asd.SIspk.normpfMap, asd.SIspk.pfLoc, asd.SIspk.pfSize, 'PC', fsave, fname, fclose )

            if plotNon_norm
                fname = [sdir 'asd_SI_bitsperspk/' fname_pref '_spkRaster_PCs'];
                plotRaster( Npcs, asd.SIspk.spkRaster_pc, ytick_files, 'PC', fsave, fname, fclose );
            end

            fname = [sdir 'asd_SI_bitsperspk/' fname_pref  '_populSummary'];
            plot_populSummary( asd.SIspk.sort_normpfMap, asd.SIspk.sort_pfMap, asd.SIspk.spkPeak, asd.SIspk.spkMean, asd.SIspk.infoMap,...
                asd.SIspk.pfSize, asd.SIspk.pfLoc, PFdata.bin_phi, fsave, fname, fclose );

            if Nepochs > 1
                fname = [sdir 'asd_SI_bitsperspk/' fname_pref '_remapping_asd_SIspk'];
                plotRemapping( asd.SIspk.normpfMap, asd.SIspk.sortIdx, fname );
            end
        end

        % spike raster plots for NON-place cells 
        if ~isempty(asd.SIspk.nonpcIdx)
            Nnonpcs = numel(asd.SIspk.nonpcIdx);
            fname = [sdir 'asd_SI_bitsperspk/' fname_pref '_normspkRaster_nonPCs'];
            plotRaster( Nnonpcs, asd.SIspk.normspkRaster_nonpc, ytick_files, 'Non-PC', fsave, fname, fclose );

            if plotNon_norm
                fname = [sdir 'asd_SI_bitsperspk/' fname_pref '_spkRaster_nonPCs'];
                plotRaster( Nnonpcs, asd.SIspk.spkRaster_nonpc, ytick_files, 'Non-PC', fsave, fname, fclose );
            end

            fname = [sdir 'asd_SI_bitsperspk/' fname_pref '_rateplot_nonPCs'];
            plot_pfLoc( Nnonpcs, asd.normrMap(asd.SIspk.nonpcIdx,:), [], [], 'Non-PC', fsave, fname, fclose )
        end
    end
    
%% PLOTTING FUNCTIONS  
% spike rasters
function plotRaster( Ncells, spkRaster, ytick_files, title_str, fsave, fname, fclose )
    [nRow, nCol] = getnRownCol(Ncells);
    nPlot = nRow*nCol;
    Nfig = ceil(Ncells/nPlot)-1;
    if Nfig<0, Nfig = 0; end

    for ii=0:Nfig
        fh = figure;
        ha = tight_subplot(nRow,nCol,[.02 .015],[.01 .05],[.01 .01]);
        for jj=0:nPlot-1
            if (ii*nPlot+jj+1) <= Ncells
                axes(ha(+jj+1));
                map = viridisMap; 
                imagesc(spkRaster{ii*nPlot+jj+1}); colormap(map);
                yticks(ytick_files); yticklabels(ytick_files); 
                xticks([]); 
                title([title_str ' ' num2str(ii*nPlot+jj+1)],'fontsize',12);
%                 hold on
%                 yyaxis right; plot(hist.meanspkRaster_MI(ii*nPlot+jj+1,:),'w-'); hold off
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
function plot_pfLoc( Ncells, normpfMap_sm, pfLoc, pfBins, title_str, fsave, fname, fclose )
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
                plot(normpfMap_sm(ii*nPlot+jj+1,:)); colormap(map);
                yticks([]); xticks([]); 
                title([title_str ' ' num2str(ii*nPlot+jj+1)],'fontsize',12);
                if ~isempty(pfLoc)
                    hold on
                    plot(pfLoc(ii*nPlot+jj+1),1,'kv','markerfacecolor','k','markersize',8); 
                    for i = 1:length(pfBins{ii*nPlot+jj+1})
                        plot(pfBins{ii*nPlot+jj+1}(i),0.5,'r*','markerfacecolor','r','markersize',4)
                    end
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

% Row 1: bin dwell time per trial, total bin dwell time, % of bin dwell
% time cell is active, % of laps cell is active
% Row 2: sorted pf maps, event rate histogram (pcs & nonpcs), spike
% amplitude rate histogram (pcs & nonpcs), histogram of info (pcs & nonpcs)
% Row 3: sorted norm pf maps, histogram of fieldsize, histogram of field
% loc, frac of pf dwell time cell is active
function plot_populSummary(bintime_trials, bintime, bin_activet, activetrials,...
        sort_pfMap_sm, sort_normpfMap_sm, spk_eventrate, spk_rate, infoMap, ...
        fieldSize, pfLoc, pf_activet, pcIdx, nonpcIdx, fsave, fname, fclose )
    % summary of population data
    Ncells = size(sort_normpfMap_sm,1);
    fh = figure;
    subplot(341); imagesc(bintime_trials); 
        xticks([1 Nbins]); xticklabels([1 100]); xlabel('Position (cm)');
        if size(bintime_trials,1)>1 
            yticks([1 size(bintime_trials,1)]); yticklabels([1 size(bintime_trials,1)]);
        else 
            yticks(1); yticklabels(1); 
        end
        ylabel('Lap #'); 
        title('Active dwell time (s)'); c = colorbar; c.Ticks = [0 round(max(max(bintime_trials)),1)];
    subplot(342); bar(bintime,'FaceColor',[0.5 0.5 0.5],'LineStyle','none');
        xticks([1 Nbins]); xticklabels([1 100]); xlabel('Position (cm)');
        ylabel('Total dwell time (s)');
    subplot(343); 
        hold on
        histogram(bin_activet(pcIdx)*100,'Normalization','probability','DisplayStyle','bar',...
            'BinMethod','sqrt','FaceColor',[0.4 0.6 1],'LineStyle','none',...
            'FaceAlpha',0.7)
        histogram(bin_activet(nonpcIdx)*100,'Normalization','probability','DisplayStyle','bar',...
            'BinMethod','sqrt','FaceColor',[1 0.6 0.6],'LineStyle','none',...
            'FaceAlpha',0.7)
        plot(mean(bin_activet(pcIdx))*100,0.2,'v','markerfacecolor',[0.4 0.6 1],...
            'markeredgecolor',[0.4 0.6 1],'markersize',4);
        plot(mean(bin_activet(nonpcIdx))*100,0.2,'v','markerfacecolor',[1 0.6 0.6],...
        'markeredgecolor',[1 0.6 0.6],'markersize',4);
        legend('pc','nonpc'); legend('boxoff');
        hold off;
        xlabel('Active bin time (%)'); 
        ylabel('Prop. of cells');
    subplot(344); 
        hold on
        histogram(activetrials(pcIdx),'Normalization','probability','DisplayStyle','bar',...
            'BinMethod','sqrt','FaceColor',[0.4 0.6 1],'LineStyle','none',...
            'FaceAlpha',0.7)
        histogram(activetrials(nonpcIdx),'Normalization','probability','DisplayStyle','bar',...
            'BinMethod','sqrt','FaceColor',[1 0.6 0.6],'LineStyle','none',...
            'FaceAlpha',0.7)
        plot(mean(activetrials(pcIdx)),0.2,'v','markerfacecolor',[0.4 0.6 1],...
            'markeredgecolor',[0.4 0.6 1],'markersize',4);
        plot(mean(activetrials(nonpcIdx)),0.2,'v','markerfacecolor',[1 0.6 0.6],...
            'markeredgecolor',[1 0.6 0.6],'markersize',4);
        hold off;
        ylabel('Prop. of cells'); 
        xlabel('Active trials (%)');
        % legend('pc','nonpc'); legend('boxoff');
    subplot(345); 
        map = viridisMap;
        imagesc(sort_normpfMap_sm); colormap(map);
        xticks([1 Nbins/2 Nbins]); xticklabels([1 50 100]); xlabel('Position (cm)');
        if Ncells>1 
            yticks([1 Ncells]); yticklabels([1 Ncells]);
        else 
            yticks(1); yticklabels(1); 
        end
        ylabel('Cell #'); 
        title('Normalised PF map'); colorbar;
    subplot(245); imagesc(sort_pfMap_sm(:,:,e)); colormap(map);
        xticks([1 Nbins/2 Nbins]); xticklabels([1 50 100]); xlabel('Position (cm)');
        if Ncells>1 
            yticks([1 Ncells]); yticklabels([1 Ncells ]);
        else 
            yticks(1); yticklabels(1); 
        end
        ylabel('Cell #'); 
        title('PF map'); colorbar;
    subplot(242); histogram(spkPeak,'Normalization','probability','FaceColor',[0.5 0.5 0.5],...
        'LineStyle','none','BinMethod','fd'); 
        hold on;
        plot(mean(spkPeak),0.2,'kv','markerfacecolor','k','markersize',6); hold off;
        xlabel('PSD peak'); 
        ylabel('Fraction of Cells');
    subplot(243); histogram(spkMean,'Normalization','probability','FaceColor',[0.5 0.5 0.5],...
        'LineStyle','none','BinMethod','fd'); 
        hold on;
        plot(mean(spkMean),0.3,'kv','markerfacecolor','k','markersize',6); hold off;
        xlabel('PSD mean'); 
        ylabel('Fraction of Cells');
    subplot(244); histogram(infoMap,'Normalization','probability','FaceColor',[0.5 0.5 0.5],...
        'LineStyle','none','BinMethod','fd'); 
        hold on;
        plot(mean(infoMap),0.3,'kv','markerfacecolor','k','markersize',6); hold off;
        xlabel('Info content');
        ylabel('Fraction of Cells');
    subplot(246); histogram(fieldSize(:,e),'Normalization','probability','FaceColor',[0.5 0.5 0.5],...
        'LineStyle','none','BinMethod','sqrt');
        hold on;
        mean_fieldsize = mean(fieldSize(:,e));
        plot(mean_fieldsize,0.2,'kv','markerfacecolor','k','markersize',6); hold off;
        xlabel('Field size(cm)');
        ylabel('Fraction of Cells');
    subplot(247); histogram(pfLoc(:,e),'FaceColor',[0.5 0.5 0.5],...
        'LineStyle','none','BinMethod','sqrt');
        hold on;
        mean_pfLoc = mean(pfLoc(:,e));
        plot(mean_pfLoc,10,'kv','markerfacecolor','k','markersize',6); hold off;
        xlabel('Position (cm)');
        ylabel('No. of PCs');
    subplot(248); histogram(bin_phi(:,e),Nbins,'Normalization','probability','FaceColor',[0.5 0.5 0.5],...
        'LineStyle','none');
        % hold on;
        % bin_phi_mean = mean(bin_phi(:,e));
        % plot(bin_phi_mean,0.03,'kv','markerfacecolor','k','markersize',6); hold off;
        xlabel('Position (cm)');
        ylabel('Fraction of Time');

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

end