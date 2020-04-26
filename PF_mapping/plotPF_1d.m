function plotPF_1d(hist, asd, PFdata, fclose, fsave, sdir, fname_pref)
    if nargin < 7, fsave = false; fname_pref = '_'; sdir = '_'; end
    if nargin < 6, fsave = false; end
    if nargin < 5, fsave = false; end
    if nargin < 4, fclose = true; fsave = false; end

    plotNon_norm = true;
    Nepochs = size(PFdata.occMap,3);
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
        plot_pfChar( Npcs, hist.SIsec.normpfMap_sm, hist.SIsec.pfLoc, hist.SIsec.pfSize, 'PC', fsave, fname, fclose )
        
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
        plot_pfChar( Nnonpcs, hist.normrMap_sm(hist.SIsec.nonpcIdx,:), [], [], 'Non-PC', fsave, fname, fclose )
    end
    

    % asd
    % plots for PLACE cells 
    if fsave
        if ~exist([sdir 'asd_SI_bitspersec/'],'dir'), mkdir([sdir 'asd_SI_bitspersec/']); end
    end
    if ~isempty(asd.SIsec.pcIdx)
        Npcs = numel(asd.SIsec.pcIdx);
        fname = [sdir 'asd_SI_bitspersec/' fname_pref '_normspkRaster_PCs'];
        plotRaster( Npcs, asd.SIsec.normspkRaster_pc, ytick_files, 'PC', fsave, fname, fclose );
        
        fname = [sdir 'asd_SI_bitspersec/' fname_pref '_rateplot_PCs'];
        plot_pfChar( Npcs, asd.SIsec.normpfMap, asd.SIsec.pfLoc, asd.SIsec.pfSize, 'PC', fsave, fname, fclose )
        
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
        plot_pfChar( Nnonpcs, asd.normrMap(asd.SIsec.nonpcIdx,:), [], [], 'Non-PC', fsave, fname, fclose )
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
        plot_pfChar( Npcs, hist.SIspk.normpfMap_sm, hist.SIspk.pfLoc, hist.SIspk.pfSize, 'PC', fsave, fname, fclose )
        
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
        plot_pfChar( Nnonpcs, hist.normrMap_sm(hist.SIspk.nonpcIdx,:), [], [], 'Non-PC', fsave, fname, fclose )
    end
    

    % asd
    % plots for PLACE cells 
    if fsave
        if ~exist([sdir 'asd_SI_bitsperspk/'],'dir'), mkdir([sdir 'asd_SI_bitsperspk/']); end
    end
    if ~isempty(asd.SIspk.pcIdx)
        Npcs = numel(asd.SIspk.pcIdx);
        fname = [sdir 'asd_SI_bitsperspk/' fname_pref '_normspkRaster_PCs'];
        plotRaster( Npcs, asd.SIspk.normspkRaster_pc, ytick_files, 'PC', fsave, fname, fclose );
        
        fname = [sdir 'asd_SI_bitsperspk/' fname_pref '_rateplot_PCs'];
        plot_pfChar( Npcs, asd.SIspk.normpfMap, asd.SIspk.pfLoc, asd.SIspk.pfSize, 'PC', fsave, fname, fclose )
        
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
        plot_pfChar( Nnonpcs, asd.normrMap(asd.SIspk.nonpcIdx,:), [], [], 'Non-PC', fsave, fname, fclose )
    end
    
function plotRaster( Ncells, spkRaster, ytick_files, title_str, fsave, fname, fclose )
    [nRow, nCol] = getnRownCol(Ncells);
    nPlot = nRow*nCol;
    Nfig = ceil((Ncells/nPlot))-1;
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

function plot_pfChar( Ncells, normpfMap_sm, pfLoc, pfSize, title_str, fsave, fname, fclose )
    [nRow, nCol] = getnRownCol(Ncells);
    nPlot = nRow*nCol;
    Nfig = ceil((Ncells/nPlot))-1;
    if Nfig<0, Nfig = 0; end

    for ii=0:Nfig
        fh = figure;
        ha = tight_subplot(nRow,nCol,[.02 .015],[.01 .05],[.01 .01]);
        for jj=0:nPlot-1
            if (ii*nPlot+jj+1) <= Ncells
                axes(ha(+jj+1));
                map = viridisMap; 
                plot(normpfMap_sm(ii*nPlot+jj+1,:)); colormap(map);
                yticks([]); xticks([]); 
                title([title_str ' ' num2str(ii*nPlot+jj+1)],'fontsize',12);
                if ~isempty(pfLoc)
                    hold on
                    plot(pfLoc(ii*nPlot+jj+1),1,'kv','markerfacecolor','k','markersize',8); 
                    plot([1, pfSize(ii*nPlot+jj+1)*Nbins/103],[0.5,0.5]); 
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

function plot_populSummary(sort_pfMap_sm, sort_normpfMap_sm, spkPeak, spkMean, infoMap, pfSize, pfLoc, bin_phi, fsave, fname, fclose )
    % summary of population data
    Ncells = size(sort_normpfMap_sm,1);
    for e = 1:Nepochs
        fh = figure;
        map = viridisMap; 
        subplot(241); imagesc(sort_normpfMap_sm(:,:,e)); colormap(map);
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
        subplot(246); histogram(pfSize(:,e),'Normalization','probability','FaceColor',[0.5 0.5 0.5],...
            'LineStyle','none','BinMethod','sqrt');
            hold on;
            mean_fieldsize = mean(pfSize(:,e));
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