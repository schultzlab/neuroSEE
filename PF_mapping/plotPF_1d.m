function plotPF_1d(hist, asd, PFdata, fsave, sdir, fname_pref, fclose)
    if nargin < 7, fclose = false; end
    if nargin < 6, fsave = false; end
    if nargin < 5, fsave = false; end
    if nargin < 4, fsave = false; end
    
    plotNon_norm = true;
    Nepochs = size(PFdata.occMap,3);
    Nposbins = size(PFdata.occMap,2);
    ytick_files = PFdata.ytick_files;
    
%% MI (bits)
    % histogram estimation
    % spike raster plots for PLACE cells 
%     if ~exist([sdir 'hist_MI/'],'dir'), mkdir([sdir 'hist_MI/']); end
%     if ~isempty(hist.pcIdx_MI)
%         Npcs = numel(hist.pcIdx_MI);
%         fname = [sdir 'hist_MI/' fname_pref '_normspkRaster_PCs_hist_MI'];
%         plotRaster( Npcs, hist.normspkRaster_MI_pc, ytick_files, 'PC', fsave, fname, fclose );
%         
%         if plotNon_norm
%             fname = [sdir 'hist_MI/' fname_pref '_spkRaster_PCs_hist_MI'];
%             plotRaster( Npcs, hist.spkRaster_MI_pc, ytick_files, 'PC', fsave, fname, fclose );
%         end
%         
%         fname = [sdir 'hist_MI/' fname_pref  '_populSummary'];
%         plot_populSummary( hist.sort_normpfMap_MI, hist.spkPeak_MI, hist.spkMean_MI, hist.infoMap_MI,...
%             hist.fieldsize_MI, hist.centroid_MI, PFdata.bin_phi, fsave, fname, fclose, hist.sort_normpfMap_MI_sm );
%             
%         if Nepochs > 1
%             fname = [sdir 'hist_MI/' fname_pref '_remapping_hist_MI'];
%             plotRemapping( normpfMap, sortIdx, fname );
%         end
%     end
%     
%     % spike raster plots for NON-place cells 
%     if ~isempty(hist.nonpcIdx_MI)
%         Nnonpcs = numel(hist.nonpcIdx_MI);
%         fname = [sdir 'hist_MI/' fname_pref '_normspkRaster_nonPCs_hist_MI'];
%         plotRaster( Nnonpcs, hist.normspkRaster_MI_nonpc, ytick_files, 'Non-PC', fsave, fname, fclose );
%         
%         if plotNon_norm
%             fname = [sdir 'hist_MI/' fname_pref '_spkRaster_nonPCs_hist_MI'];
%             plotRaster( Nnonpcs, hist.spkRaster_MI_nonpc, ytick_files, 'Non-PC', fsave, fname, fclose );
%         end
%     end
%     
% 
%     % asd, MI (bits)
%     % spike raster plots for PLACE cells 
%     if ~exist([sdir 'asd_MI/'],'dir'), mkdir([sdir 'asd_MI/']); end
%     if ~isempty(asd.pcIdx_MI)
%         Npcs = numel(asd.pcIdx_MI);
%         fname = [sdir 'asd_MI/' fname_pref '_normspkRaster_PCs_asd_MI'];
%         plotRaster( Npcs, asd.normspkRaster_MI_pc, ytick_files, 'PC', fsave, fname, fclose );
%         
%         if plotNon_norm
%             fname = [sdir 'asd_MI/' fname_pref '_spkRaster_PCs_asd_MI'];
%             plotRaster( Npcs, asd.spkRaster_MI_pc, ytick_files, 'PC', fsave, fname, fclose );
%         end
%         
%         fname = [sdir 'asd_MI/' fname_pref  '_populSummary'];
%         plot_populSummary( asd.sort_normpfMap_MI, asd.spkPeak_MI, asd.spkMean_MI, asd.infoMap_MI,...
%             asd.fieldsize_MI, asd.centroid_MI, PFdata.bin_phi, fsave, fname, fclose );
%             
%         if Nepochs > 1
%             fname = [sdir 'asd_MI/' fname_pref '_remapping_asd_MI'];
%             plotRemapping( normpfMap, sortIdx, fname );
%         end
%     end
%     
%     % spike raster plots for NON-place cells 
%     if ~isempty(asd.nonpcIdx_MI)
%         Nnonpcs = numel(asd.nonpcIdx_MI);
%         fname = [sdir 'asd_MI/' fname_pref '_normspkRaster_nonPCs_asd_MI'];
%         plotRaster( Nnonpcs, asd.normspkRaster_MI_nonpc, ytick_files, 'Non-PC', fsave, fname, fclose );
%         
%         if plotNon_norm
%             fname = [sdir 'asd_MI/' fname_pref '_spkRaster_nonPCs_asd_MI'];
%             plotRaster( Nnonpcs, asd.spkRaster_MI_nonpc, ytick_files, 'Non-PC', fsave, fname, fclose );
%         end
%     end
   
%% SI (bits/sec)
    % histogram estimation
    % spike raster plots for PLACE cells 
    if ~exist([sdir 'hist_SIsec/'],'dir'), mkdir([sdir 'hist_SIsec/']); end
    if ~isempty(hist.pcIdx_SIsec)
        Npcs = numel(hist.pcIdx_SIsec);
        fname = [sdir 'hist_SIsec/' fname_pref '_normspkRaster_PCs_hist_SIsec'];
        plotRaster( Npcs, hist.normspkRaster_SIsec_pc, ytick_files, 'PC', fsave, fname, fclose );
        
        if plotNon_norm
            fname = [sdir 'hist_SIsec/' fname_pref '_spkRaster_PCs_hist_SIsec'];
            plotRaster( Npcs, hist.spkRaster_SIsec_pc, ytick_files, 'PC', fsave, fname, fclose );
        end
        
        fname = [sdir 'hist_SIsec/' fname_pref  '_populSummary'];
        plot_populSummary( hist.sort_normpfMap_SIsec, hist.spkPeak_SIsec, hist.spkMean_SIsec, hist.infoMap_SIsec,...
            hist.fieldsize_SIsec, hist.centroid_SIsec, PFdata.bin_phi, fsave, fname, fclose, hist.sort_normpfMap_SIsec_sm );
            
        if Nepochs > 1
            fname = [sdir 'hist_SIsec/' fname_pref '_remapping_hist_SIsec'];
            plotRemapping( normpfMap, sortIdx, fname );
        end
    end
    
    % spike raster plots for NON-place cells 
    if ~isempty(hist.nonpcIdx_SIsec)
        Nnonpcs = numel(hist.nonpcIdx_SIsec);
        fname = [sdir 'hist_SIsec/' fname_pref '_normspkRaster_nonPCs_hist_SIsec'];
        plotRaster( Nnonpcs, hist.normspkRaster_SIsec_nonpc, ytick_files, 'Non-PC', fsave, fname, fclose );
        
        if plotNon_norm
            fname = [sdir 'hist_SIsec/' fname_pref '_spkRaster_nonPCs_hist_SIsec'];
            plotRaster( Nnonpcs, hist.spkRaster_SIsec_nonpc, ytick_files, 'Non-PC', fsave, fname, fclose );
        end
    end
    

    % asd
    % spike raster plots for PLACE cells 
%     if ~exist([sdir 'asd_SIsec/'],'dir'), mkdir([sdir 'asd_SIsec/']); end
%     if ~isempty(asd.pcIdx_SIsec)
%         Npcs = numel(asd.pcIdx_SIsec);
%         fname = [sdir 'asd_SIsec/' fname_pref '_normspkRaster_PCs_asd_SIsec'];
%         plotRaster( Npcs, asd.normspkRaster_SIsec_pc, ytick_files, 'PC', fsave, fname, fclose );
%         
%         if plotNon_norm
%             fname = [sdir 'asd_SIsec/' fname_pref '_spkRaster_PCs_asd_SIsec'];
%             plotRaster( Npcs, asd.spkRaster_SIsec_pc, ytick_files, 'PC', fsave, fname, fclose );
%         end
%         
%         fname = [sdir 'asd_SIsec/' fname_pref  '_populSummary'];
%         plot_populSummary( asd.sort_normpfMap_SIsec, asd.spkPeak_SIsec, asd.spkMean_SIsec, asd.infoMap_SIsec,...
%             asd.fieldsize_SIsec, asd.centroid_SIsec, PFdata.bin_phi, fsave, fname, fclose );
%             
%         if Nepochs > 1
%             fname = [sdir 'asd_SIsec/' fname_pref '_remapping_asd_SIsec'];
%             plotRemapping( normpfMap, sortIdx, fname );
%         end
%     end
    
    % spike raster plots for NON-place cells 
%     if ~isempty(asd.nonpcIdx_SIsec)
%         Nnonpcs = numel(asd.nonpcIdx_SIsec);
%         fname = [sdir 'asd_SIsec/' fname_pref '_normspkRaster_nonPCs_asd_SIsec'];
%         plotRaster( Nnonpcs, asd.normspkRaster_SIsec_nonpc, ytick_files, 'Non-PC', fsave, fname, fclose );
%         
%         if plotNon_norm
%             fname = [sdir 'asd_SIsec/' fname_pref '_spkRaster_nonPCs_asd_SIsec'];
%             plotRaster( Nnonpcs, asd.spkRaster_SIsec_nonpc, ytick_files, 'Non-PC', fsave, fname, fclose );
%         end
%     end

%% SI (bits/spike)
    % histogram estimation
    % spike raster plots for PLACE cells 
%     if ~exist([sdir 'hist_SIspk/'],'dir'), mkdir([sdir 'hist_SIspk/']); end
%     if ~isempty(hist.pcIdx_SIspk)
%         Npcs = numel(hist.pcIdx_SIspk);
%         fname = [sdir 'hist_SIspk/' fname_pref '_normspkRaster_PCs_hist_SIspk'];
%         plotRaster( Npcs, hist.normspkRaster_SIspk_pc, ytick_files, 'PC', fsave, fname, fclose );
%         
%         if plotNon_norm
%             fname = [sdir 'hist_SIspk/' fname_pref '_spkRaster_PCs_hist_SIspk'];
%             plotRaster( Npcs, hist.spkRaster_SIspk_pc, ytick_files, 'PC', fsave, fname, fclose );
%         end
%         
%         fname = [sdir 'hist_SIspk/' fname_pref  '_populSummary'];
%         plot_populSummary( hist.sort_normpfMap_SIspk, hist.spkPeak_SIspk, hist.spkMean_SIspk, hist.infoMap_SIspk,...
%             hist.fieldsize_SIspk, hist.centroid_SIspk, PFdata.bin_phi, fsave, fname, fclose, hist.sort_normpfMap_SIspk_sm );
%             
%         if Nepochs > 1
%             fname = [sdir 'hist_SIspk/' fname_pref '_remapping_hist_SIspk'];
%             plotRemapping( normpfMap, sortIdx, fname );
%         end
%     end
    
    % spike raster plots for NON-place cells 
%     if ~isempty(hist.nonpcIdx_SIspk)
%         Nnonpcs = numel(hist.nonpcIdx_SIspk);
%         fname = [sdir 'hist_SIspk/' fname_pref '_normspkRaster_nonPCs_hist_SIspk'];
%         plotRaster( Nnonpcs, hist.normspkRaster_SIspk_nonpc, ytick_files, 'Non-PC', fsave, fname, fclose );
%         
%         if plotNon_norm
%             fname = [sdir 'hist_SIspk/' fname_pref '_spkRaster_nonPCs_hist_SIspk'];
%             plotRaster( Nnonpcs, hist.spkRaster_SIspk_nonpc, ytick_files, 'Non-PC', fsave, fname, fclose );
%         end
%     end
%     

    % asd
    % spike raster plots for PLACE cells 
%     if ~exist([sdir 'asd_SIspk/'],'dir'), mkdir([sdir 'asd_SIspk/']); end
%     if ~isempty(asd.pcIdx_SIspk)
%         Npcs = numel(asd.pcIdx_SIspk);
%         fname = [sdir 'asd_SIspk/' fname_pref '_normspkRaster_PCs_asd_SIspk'];
%         plotRaster( Npcs, asd.normspkRaster_SIspk_pc, ytick_files, 'PC', fsave, fname, fclose );
%         
%         if plotNon_norm
%             fname = [sdir 'asd_SIspk/' fname_pref '_spkRaster_PCs_asd_SIspk'];
%             plotRaster( Npcs, asd.spkRaster_SIspk_pc, ytick_files, 'PC', fsave, fname, fclose );
%         end
%         
%         fname = [sdir 'asd_SIspk/' fname_pref  '_populSummary'];
%         plot_populSummary( asd.sort_normpfMap_SIspk, asd.spkPeak_SIspk, asd.spkMean_SIspk, asd.infoMap_SIspk,...
%             asd.fieldsize_SIspk, asd.centroid_SIspk, PFdata.bin_phi, fsave, fname, fclose );
%             
%         if Nepochs > 1
%             fname = [sdir 'asd_SIspk/' fname_pref '_remapping_asd_SIspk'];
%             plotRemapping( normpfMap, sortIdx, fname );
%         end
%     end
    
    % spike raster plots for NON-place cells 
%     if ~isempty(asd.nonpcIdx_SIspk)
%         Nnonpcs = numel(asd.nonpcIdx_SIspk);
%         fname = [sdir 'asd_SIspk/' fname_pref '_normspkRaster_nonPCs_asd_SIspk'];
%         plotRaster( Nnonpcs, asd.normspkRaster_SIspk_nonpc, ytick_files, 'Non-PC', fsave, fname, fclose );
%         
%         if plotNon_norm
%             fname = [sdir 'asd_SIspk/' fname_pref '_spkRaster_nonPCs_asd_SIspk'];
%             plotRaster( Nnonpcs, asd.spkRaster_SIspk_nonpc, ytick_files, 'Non-PC', fsave, fname, fclose );
%         end
%     end
    
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

function plot_populSummary(sort_normpfMap, spkPeak, spkMean, infoMap, fieldsize, centroids, bin_phi, fsave, fname, fclose, sort_normpfMap_sm )
    % summary of population data
    Ncells = size(sort_normpfMap,1);
    for e = 1:Nepochs
        fh = figure;
        map = viridisMap; 
        subplot(241); imagesc(sort_normpfMap(:,:,e)); colormap(map);
            xticks([1 Nposbins/2 Nposbins]); xticklabels([1 50 100]); xlabel('Position (cm)');
            if Ncells>1 
                yticks([1 Ncells]); yticklabels([1 Ncells]);
            else 
                yticks(1); yticklabels(1); 
            end
            ylabel('Cell #'); 
            title('PF map'); colorbar;
        if nargin > 10
            subplot(245); imagesc(sort_normpfMap_sm(:,:,e)); colormap(map);
            xticks([1 Nposbins/2 Nposbins]); xticklabels([1 50 100]); xlabel('Position (cm)');
            if Ncells>1 
                yticks([1 Ncells]); yticklabels([1 Ncells ]);
            else 
                yticks(1); yticklabels(1); 
            end
            ylabel('Cell #'); 
            title('Smoothened PF map'); colorbar;
        end
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
        subplot(246); histogram(fieldsize(:,e),'Normalization','probability','FaceColor',[0.5 0.5 0.5],...
            'LineStyle','none','BinMethod','sqrt');
            hold on;
            mean_fieldsize = mean(fieldsize(:,e));
            plot(mean_fieldsize,0.2,'kv','markerfacecolor','k','markersize',6); hold off;
            xlabel('Field size(cm)');
            ylabel('Fraction of Cells');
        subplot(247); histogram(centroids(:,e),Nposbins,'FaceColor',[0.5 0.5 0.5],...
            'LineStyle','none');
            hold on;
            mean_centroid = mean(centroids(:,e));
            plot(mean_centroid,10,'kv','markerfacecolor','k','markersize',6); hold off;
            xlabel('Position (cm)');
            ylabel('No. of PCs');
        subplot(248); histogram(bin_phi(:,e),Nposbins,'Normalization','probability','FaceColor',[0.5 0.5 0.5],...
            'LineStyle','none');
            hold on;
            bin_phi_mean = mean(bin_phi(:,e));
            plot(bin_phi_mean,0.03,'kv','markerfacecolor','k','markersize',6); hold off;
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