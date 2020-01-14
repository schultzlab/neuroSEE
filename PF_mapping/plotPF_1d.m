function plotPF_1d(hist, asd, PFdata, fsave, sdir, fname_pref, fclose)
    if nargin < 7, fclose = false; end
    if nargin < 6, fsave = false; end
    if nargin < 5, fsave = false; end
    if nargin < 4, fsave = false; end
    
    Nepochs = size(PFdata.occMap,3);
    ytick_files = PFdata.ytick_files;
    
    %% histogram estimation, MI (bits)
    Npcs = numel(hist.pcIdx_MI);

    if Npcs > 0
        % spike raster plots for place cells
        [nRow, nCol] = getnRownCol(Npcs);
        nPlot = nRow*nCol;
        Nfig = round((Npcs/nPlot))-1;
        if Nfig<0, Nfig = 0; end
        if ~exist([sdir 'hist_MI/'],'dir'), mkdir([sdir 'hist_MI/']); end

        for ii=0:Nfig
            fh = figure;
            ha = tight_subplot(nRow,nCol,[.02 .015],[.01 .05],[.01 .01]);
            for jj=0:nPlot-1
                if (ii*nPlot+jj+1) <= Npcs
                    axes(ha(+jj+1));
                    imagesc(hist.normspkRaster_MI_pc{ii*nPlot+jj+1}); 
                    yticks(ytick_files); yticklabels(ytick_files); % ylabel('Trial #');
                    xticks([]); % xticklabels([1 50 100]); xlabel('Position (cm)');
                    title(['PC ' num2str(ii*nPlot+jj+1)],'fontsize',15);
    %                 hold on
    %                 yyaxis right; plot(hist.meanspkRaster_MI(ii*nPlot+jj+1,:),'w-'); hold off
                end
            end
            if fsave
                if Npcs/nPlot <= 1
                    fname_fig = [sdir 'hist_MI/' fname_pref '_normspkRaster_PCs_hist_MI'];
                else
                    fname_fig = [sdir 'hist_MI/' fname_pref '_normspkRaster_PCs_hist_MI_' num2str(ii+1)];
                end
                savefig( fh, fname_fig );
                saveas( fh, fname_fig, 'png' );
                if fclose, close( fh ); end
            end
        end 

        for ii=0:Nfig
            fh = figure;
            ha = tight_subplot(nRow,nCol,[.02 .015],[.01 .05],[.01 .01]);
            for jj=0:nPlot-1
                if (ii*nPlot+jj+1) <= Npcs
                    axes(ha(+jj+1));
                    imagesc(hist.spkRaster_MI_pc{ii*nPlot+jj+1}); 
                    yticks(ytick_files); yticklabels(ytick_files); % ylabel('Trial #');
                    xticks([]); % xticklabels([1 50 100]); xlabel('Position (cm)');
                    title(['PC ' num2str(ii*nPlot+jj+1)],'fontsize',15);
                    colorbar;
    %                 hold on
    %                 yyaxis right; plot(hist.meanspkRaster_MI(ii*nPlot+jj+1,:),'w-'); hold off
                end
            end
            if fsave
                if Npcs/nPlot <= 1
                    fname_fig = [sdir 'hist_MI/' fname_pref '_spkRaster_PCs_hist_MI'];
                else
                    fname_fig = [sdir 'hist_MI/' fname_pref '_spkRaster_PCs_hist_MI_' num2str(ii+1)];
                end
                savefig( fh, fname_fig );
                saveas( fh, fname_fig, 'png' );
                if fclose, close( fh ); end
            end
        end 
        
        % summary of population data
        for e = 1:Nepochs
            fh = figure;
            subplot(241); imagesc(hist.sort_normpfMap_MI(:,:,e));
                xticks([1 15 30]); xticklabels([1 50 100]); xlabel('Position (cm)');
                if Npcs>1 
                    yticks([1 Npcs]); yticklabels([Npcs 1]);
                else 
                    yticks(1); yticklabels(1); 
                end
                ylabel('Cell #'); 
                title('PF map'); colorbar;
            subplot(245); imagesc(hist.sort_normpfMap_MI_sm(:,:,e));
                xticks([1 15 30]); xticklabels([1 50 100]); xlabel('Position (cm)');
                if Npcs>1 
                    yticks([1 Npcs]); yticklabels([Npcs 1]);
                else 
                    yticks(1); yticklabels(1); 
                end
                ylabel('Cell #'); 
                title('Smoothened PF map'); colorbar;
            subplot(242); histogram(hist.spkPeak_MI,100,'Normalization','probability'); 
                hold on;
                plot(mean(hist.spkPeak_MI),0.06,'kv','markerfacecolor','k','markersize',6); hold off;
                xlabel('PSD peak'); set(gca,'xscale','log');
                ylabel('Cells (%)');
            subplot(243); histogram(hist.spkMean_MI,40,'Normalization','probability');
                hold on;
                plot(mean(hist.spkMean_MI),0.08,'kv','markerfacecolor','k','markersize',6); hold off;
                xlabel('PSD mean');
                ylabel('Cells (%)');
            subplot(244); histogram(hist.infoMap_MI,40,'Normalization','probability');
                hold on;
                plot(mean(hist.infoMap_MI),0.12,'kv','markerfacecolor','k','markersize',6); hold off;
                xlabel('Mutual info (bits)');
                ylabel('Cells (%)');
            subplot(246); histogram(hist.fieldsize_MI(:,e),40,'Normalization','probability');
                hold on;
                mean_fieldsize_MI = mean(hist.fieldsize_MI(:,e));
                plot(mean_fieldsize_MI,0.1,'kv','markerfacecolor','k','markersize',6); hold off;
                xlabel('Field size(cm)');
                ylabel('Cells (%)');
            subplot(247); histogram(hist.centroid_MI(:,e),30);
                hold on;
                mean_centroid_MI = mean(hist.centroid_MI(:,e));
                plot(mean_centroid_MI,10,'kv','markerfacecolor','k','markersize',6); hold off;
                xlabel('Position (cm)');
                ylabel('No. of PCs');
            subplot(248); histogram(PFdata.bin_phi(:,e),30,'Normalization','probability');
                hold on;
                mean_bin_phi = mean(PFdata.bin_phi(:,e));
                plot(mean_bin_phi,0.04,'kv','markerfacecolor','k','markersize',6); hold off;
                xlabel('Position (cm)');
                ylabel('Time (%)');
            
            if fsave
                if ~exist(sdir,'dir'), mkdir(sdir); end

                if Nepochs == 1
                    fname_fig = [sdir 'hist_MI/' fname_pref  '_populSummary'];
                else
                    fname_fig = [sdir 'hist_MI/' fname_pref  '_populSummary_' num2str(e) 'of' num2str(Nepochs) 'ep'];
                end
                savefig( fh, fname_fig );
                saveas( fh, fname_fig, 'png' );
                if fclose, close( fh ); end
            end
        end
        
        % remapping within a session
        if Nepochs > 1
            fh = figure;
            for ei = 1:Nepochs % rows: sorting
                for ej = 1:Nepochs % cols: epochs 
                    subplot(Nepochs, Nepochs, (ei-1)*Nepochs + ej); imagesc(hist.normpfMap_MI(hist.sortIdx_MI(:,ei),:,ej)); 
                    title(['Epoch ' num2str(ej)]); ylabel(['Epoch' num2str(ei) ' sorting']);
                end
            end
            if fsave
                fname_fig = [sdir 'hist_MI/' fname_pref '_remapping_hist_MI'];
                savefig( fh, fname_fig );
                saveas( fh, fname_fig, 'png' );
                if fclose, close( fh ); end
            end
        end
    end
    
    % spike raster plots for NON-place cells 
    Nnonpcs = numel(hist.nonpcIdx_MI);
    if Nnonpcs > 0
        [nRow, nCol] = getnRownCol(Nnonpcs);
        nPlot = nRow*nCol;
        Nfig = round((Nnonpcs/nPlot))-1;
        if Nfig<0, Nfig = 0; end

        for ii=0:Nfig
            fh = figure;
            ha = tight_subplot(nRow,nCol,[.02 .015],[.01 .05],[.01 .01]);
            for jj=0:nPlot-1
                if (ii*nPlot+jj+1) <= Nnonpcs
                    axes(ha(+jj+1));
                    imagesc(hist.normspkRaster_MI_nonpc{ii*nPlot+jj+1}); 
                    yticks(ytick_files); yticklabels(ytick_files); % ylabel('Trial #');
                    xticks([]); % xticklabels([1 50 100]); xlabel('Position (cm)');
                    title(['Non-PC ' num2str(ii*nPlot+jj+1)],'fontsize',15);
                end
            end
            if fsave
                if Nnonpcs/nPlot <= 1
                    fname_fig = [sdir 'hist_MI/' fname_pref '_normspkRaster_nonPCs_hist_MI'];
                else
                    fname_fig = [sdir 'hist_MI/' fname_pref '_normspkRaster_nonPCs_hist_MI_' num2str(ii+1)];
                end
                savefig( fh, fname_fig );
                saveas( fh, fname_fig, 'png' );
                if fclose, close( fh ); end
            end
        end 
    
        for ii=0:Nfig
            fh = figure;
            ha = tight_subplot(nRow,nCol,[.02 .015],[.01 .05],[.01 .01]);
            for jj=0:nPlot-1
                if (ii*nPlot+jj+1) <= Nnonpcs
                    axes(ha(+jj+1));
                    imagesc(hist.spkRaster_MI_nonpc{ii*nPlot+jj+1}); 
                    yticks(ytick_files); yticklabels(ytick_files); % ylabel('Trial #');
                    xticks([]); % xticklabels([1 50 100]); xlabel('Position (cm)');
                    title(['Non-PC ' num2str(ii*nPlot+jj+1)],'fontsize',15);
                    colorbar;
                end
            end
            if fsave
                if Nnonpcs/nPlot <= 1
                    fname_fig = [sdir 'hist_MI/' fname_pref '_spkRaster_nonPCs_hist_MI'];
                else
                    fname_fig = [sdir 'hist_MI/' fname_pref '_spkRaster_nonPCs_hist_MI_' num2str(ii+1)];
                end
                savefig( fh, fname_fig );
                saveas( fh, fname_fig, 'png' );
                if fclose, close( fh ); end
            end
        end 
    end
    
    %% histogram estimation, SI (bits/sec)
    Npcs = numel(hist.pcIdx_SIsec);

    if Npcs > 0
        % spike raster plots for place cells
        [nRow, nCol] = getnRownCol(Npcs);
        nPlot = nRow*nCol;
        Nfig = round((Npcs/nPlot))-1;
        if Nfig<0, Nfig = 0; end
        if ~exist([sdir 'hist_SIsec/'],'dir'), mkdir([sdir 'hist_SIsec/']); end

        for ii=0:Nfig
            fh = figure;
            ha = tight_subplot(nRow,nCol,[.02 .015],[.01 .05],[.01 .01]);
            for jj=0:nPlot-1
                if (ii*nPlot+jj+1) <= Npcs
                    axes(ha(+jj+1));
                    imagesc(hist.normspkRaster_SIsec_pc{ii*nPlot+jj+1}); 
                    yticks(ytick_files); yticklabels(ytick_files); % ylabel('Trial #');
                    xticks([]); % xticklabels([1 50 100]); xlabel('Position (cm)');
                    title(['PC ' num2str(ii*nPlot+jj+1)],'fontsize',15);
    %                 hold on
    %                 yyaxis right; plot(hist.meanspkRaster_SIsec(ii*nPlot+jj+1,:),'w-'); hold off
                end
            end
            if fsave
                if Npcs/nPlot <= 1
                    fname_fig = [sdir 'hist_SIsec/' fname_pref '_normspkRaster_PCs_hist_SIsec'];
                else
                    fname_fig = [sdir 'hist_SIsec/' fname_pref '_normspkRaster_PCs_hist_SIsec_' num2str(ii+1)];
                end
                savefig( fh, fname_fig );
                saveas( fh, fname_fig, 'png' );
                if fclose, close( fh ); end
            end
        end
        
        for ii=0:Nfig
            fh = figure;
            ha = tight_subplot(nRow,nCol,[.02 .015],[.01 .05],[.01 .01]);
            for jj=0:nPlot-1
                if (ii*nPlot+jj+1) <= Npcs
                    axes(ha(+jj+1));
                    imagesc(hist.spkRaster_SIsec_pc{ii*nPlot+jj+1}); 
                    yticks(ytick_files); yticklabels(ytick_files); % ylabel('Trial #');
                    xticks([]); % xticklabels([1 50 100]); xlabel('Position (cm)');
                    title(['PC ' num2str(ii*nPlot+jj+1)],'fontsize',15);
                    colorbar;
    %                 hold on
    %                 yyaxis right; plot(hist.meanspkRaster_SIsec(ii*nPlot+jj+1,:),'w-'); hold off
                end
            end
            if fsave
                if Npcs/nPlot <= 1
                    fname_fig = [sdir 'hist_SIsec/' fname_pref '_spkRaster_PCs_hist_SIsec'];
                else
                    fname_fig = [sdir 'hist_SIsec/' fname_pref '_spkRaster_PCs_hist_SIsec_' num2str(ii+1)];
                end
                savefig( fh, fname_fig );
                saveas( fh, fname_fig, 'png' );
                if fclose, close( fh ); end
            end
        end

        % summary of population data
        for e = 1:Nepochs
            fh = figure;
            subplot(241); imagesc(hist.sort_normpfMap_SIsec(:,:,e));
                xticks([1 15 30]); xticklabels([1 50 100]); xlabel('Position (cm)');
                if Npcs>1 
                    yticks([1 Npcs]); yticklabels([Npcs 1]);
                else 
                    yticks(1); yticklabels(1); 
                end
                ylabel('Cell #'); 
                title('PF map'); colorbar;
            subplot(245); imagesc(hist.sort_normpfMap_SIsec_sm(:,:,e));
                xticks([1 15 30]); xticklabels([1 50 100]); xlabel('Position (cm)');
                if Npcs>1 
                    yticks([1 Npcs]); yticklabels([Npcs 1]);
                else 
                    yticks(1); yticklabels(1); 
                end
                ylabel('Cell #'); 
                title('Smoothened PF map'); colorbar;
            subplot(242); histogram(hist.spkPeak_SIsec,100,'Normalization','probability'); 
                hold on;
                plot(mean(hist.spkPeak_SIsec),0.06,'kv','markerfacecolor','k','markersize',6); hold off;
                xlabel('PSD peak'); set(gca,'xscale','log');
                ylabel('Cells (%)');
            subplot(243); histogram(hist.spkMean_SIsec,40,'Normalization','probability');
                hold on;
                plot(mean(hist.spkMean_SIsec),0.08,'kv','markerfacecolor','k','markersize',6); hold off;
                xlabel('PSD mean');
                ylabel('Cells (%)');
            subplot(244); histogram(hist.infoMap_SIsec,40,'Normalization','probability');
                hold on;
                plot(mean(hist.infoMap_SIsec),0.12,'kv','markerfacecolor','k','markersize',6); hold off;
                xlabel('Mutual info (bits)');
                ylabel('Cells (%)');
            subplot(246); histogram(hist.fieldsize_SIsec(:,e),40,'Normalization','probability');
                hold on;
                mean_fieldsize_SIsec = mean(hist.fieldsize_SIsec(:,e));
                plot(mean_fieldsize_SIsec,0.1,'kv','markerfacecolor','k','markersize',6); hold off;
                xlabel('Field size(cm)');
                ylabel('Cells (%)');
            subplot(247); histogram(hist.centroid_SIsec(:,e),30);
                hold on;
                mean_centroid_SIsec = mean(hist.centroid_SIsec(:,e));
                plot(mean_centroid_SIsec,10,'kv','markerfacecolor','k','markersize',6); hold off;
                xlabel('Position (cm)');
                ylabel('No. of PCs');
            subplot(248); histogram(PFdata.bin_phi(:,e),30,'Normalization','probability');
                hold on;
                mean_bin_phi = mean(PFdata.bin_phi(:,e));
                plot(mean_bin_phi,0.04,'kv','markerfacecolor','k','markersize',6); hold off;
                xlabel('Position (cm)');
                ylabel('Time (%)');
           
            if fsave
                if ~exist(sdir,'dir'), mkdir(sdir); end

                if Nepochs == 1
                    fname_fig = [sdir 'hist_SIsec/' fname_pref  '_populSummary'];
                else
                    fname_fig = [sdir 'hist_SIsec/' fname_pref  '_populSummary_' num2str(e) 'of' num2str(Nepochs) 'ep'];
                end
                savefig( fh, fname_fig );
                saveas( fh, fname_fig, 'png' );
                if fclose, close( fh ); end
            end
        end
        
        % remapping within a session
        if Nepochs > 1
            fh = figure;
            for ei = 1:Nepochs % rows: sorting
                for ej = 1:Nepochs % cols: epochs 
                    subplot(Nepochs, Nepochs, (ei-1)*Nepochs + ej); imagesc(hist.normpfMap_SIsec(hist.sortIdx_SIsec(:,ei),:,ej)); 
                    title(['Epoch ' num2str(ej)]); ylabel(['Epoch' num2str(ei) ' sorting']);
                end
            end
            if fsave
                fname_fig = [sdir 'hist_SIsec/' fname_pref '_remapping_hist_SIsec'];
                savefig( fh, fname_fig );
                saveas( fh, fname_fig, 'png' );
                if fclose, close( fh ); end
            end
        end
    end
    
    % spike raster plots for NON-place cells 
    Nnonpcs = numel(hist.nonpcIdx_SIsec);
    if Nnonpcs > 0
        [nRow, nCol] = getnRownCol(Nnonpcs);
        nPlot = nRow*nCol;
        Nfig = round((Nnonpcs/nPlot))-1;
        if Nfig<0, Nfig = 0; end

        for ii=0:Nfig
            fh = figure;
            ha = tight_subplot(nRow,nCol,[.02 .015],[.01 .05],[.01 .01]);
            for jj=0:nPlot-1
                if (ii*nPlot+jj+1) <= Nnonpcs
                    axes(ha(+jj+1));
                    imagesc(hist.normspkRaster_SIsec_nonpc{ii*nPlot+jj+1}); 
                    yticks(ytick_files); yticklabels(ytick_files); % ylabel('Trial #');
                    xticks([]); % xticklabels([1 50 100]); xlabel('Position (cm)');
                    title(['Non-PC ' num2str(ii*nPlot+jj+1)],'fontsize',15);
                end
            end
            if fsave
                if Nnonpcs/nPlot <= 1
                    fname_fig = [sdir 'hist_SIsec/' fname_pref '_normspkRaster_nonPCs_hist_SIsec'];
                else
                    fname_fig = [sdir 'hist_SIsec/' fname_pref '_normspkRaster_nonPCs_hist_SIsec_' num2str(ii+1)];
                end
                savefig( fh, fname_fig );
                saveas( fh, fname_fig, 'png' );
                if fclose, close( fh ); end
            end
        end 
        
        for ii=0:Nfig
            fh = figure;
            ha = tight_subplot(nRow,nCol,[.02 .015],[.01 .05],[.01 .01]);
            for jj=0:nPlot-1
                if (ii*nPlot+jj+1) <= Nnonpcs
                    axes(ha(+jj+1));
                    imagesc(hist.spkRaster_SIsec_nonpc{ii*nPlot+jj+1}); 
                    yticks(ytick_files); yticklabels(ytick_files); % ylabel('Trial #');
                    xticks([]); % xticklabels([1 50 100]); xlabel('Position (cm)');
                    title(['Non-PC ' num2str(ii*nPlot+jj+1)],'fontsize',15);
                    colorbar;
                end
            end
            if fsave
                if Nnonpcs/nPlot <= 1
                    fname_fig = [sdir 'hist_SIsec/' fname_pref '_spkRaster_nonPCs_hist_SIsec'];
                else
                    fname_fig = [sdir 'hist_SIsec/' fname_pref '_spkRaster_nonPCs_hist_SIsec_' num2str(ii+1)];
                end
                savefig( fh, fname_fig );
                saveas( fh, fname_fig, 'png' );
                if fclose, close( fh ); end
            end
        end
    end
        
    
    %% histogram estimation, SI (bits/spike)
    Npcs = numel(hist.pcIdx_SIspk);

    if Npcs > 0
        % spike raster plots for place cells
        [nRow, nCol] = getnRownCol(Npcs);
        nPlot = nRow*nCol;
        Nfig = round((Npcs/nPlot))-1;
        if Nfig<0, Nfig = 0; end
        if ~exist([sdir 'hist_SIspk/'],'dir'), mkdir([sdir 'hist_SIspk/']); end

        for ii=0:Nfig
            fh = figure;
            ha = tight_subplot(nRow,nCol,[.02 .015],[.01 .05],[.01 .01]);
            for jj=0:nPlot-1
                if (ii*nPlot+jj+1) <= Npcs
                    axes(ha(+jj+1));
                    imagesc(hist.normspkRaster_SIspk_pc{ii*nPlot+jj+1}); 
                    yticks(ytick_files); yticklabels(ytick_files); % ylabel('Trial #');
                    xticks([]); % xticklabels([1 50 100]); xlabel('Position (cm)');
                    title(['PC ' num2str(ii*nPlot+jj+1)],'fontsize',15);
    %                 hold on
    %                 yyaxis right; plot(hist.meanspkRaster_SIspk(ii*nPlot+jj+1,:),'w-'); hold off
                end
            end
            if fsave
                if Npcs/nPlot <= 1
                    fname_fig = [sdir 'hist_SIspk/' fname_pref '_normspkRaster_PCs_hist_SIspk'];
                else
                    fname_fig = [sdir 'hist_SIspk/' fname_pref '_normspkRaster_PCs_hist_SIspk_' num2str(ii+1)];
                end
                savefig( fh, fname_fig );
                saveas( fh, fname_fig, 'png' );
                if fclose, close( fh ); end
            end
        end
        
        for ii=0:Nfig
            fh = figure;
            ha = tight_subplot(nRow,nCol,[.02 .015],[.01 .05],[.01 .01]);
            for jj=0:nPlot-1
                if (ii*nPlot+jj+1) <= Npcs
                    axes(ha(+jj+1));
                    imagesc(hist.spkRaster_SIspk_pc{ii*nPlot+jj+1}); 
                    yticks(ytick_files); yticklabels(ytick_files); % ylabel('Trial #');
                    xticks([]); % xticklabels([1 50 100]); xlabel('Position (cm)');
                    title(['PC ' num2str(ii*nPlot+jj+1)],'fontsize',15);
                    colorbar;
    %                 hold on
    %                 yyaxis right; plot(hist.meanspkRaster_SIspk(ii*nPlot+jj+1,:),'w-'); hold off
                end
            end
            if fsave
                if Npcs/nPlot <= 1
                    fname_fig = [sdir 'hist_SIspk/' fname_pref '_spkRaster_PCs_hist_SIspk'];
                else
                    fname_fig = [sdir 'hist_SIspk/' fname_pref '_spkRaster_PCs_hist_SIspk_' num2str(ii+1)];
                end
                savefig( fh, fname_fig );
                saveas( fh, fname_fig, 'png' );
                if fclose, close( fh ); end
            end
        end 

        % summary of population data
        for e = 1:Nepochs
            fh = figure;
            subplot(241); imagesc(hist.sort_normpfMap_SIspk(:,:,e));
                xticks([1 15 30]); xticklabels([1 50 100]); xlabel('Position (cm)');
                if Npcs>1 
                    yticks([1 Npcs]); yticklabels([Npcs 1]);
                else 
                    yticks(1); yticklabels(1); 
                end
                ylabel('Cell #'); 
                title('PF map'); colorbar;
            subplot(245); imagesc(hist.sort_normpfMap_SIspk_sm(:,:,e));
                xticks([1 15 30]); xticklabels([1 50 100]); xlabel('Position (cm)');
                if Npcs>1 
                    yticks([1 Npcs]); yticklabels([Npcs 1]);
                else 
                    yticks(1); yticklabels(1); 
                end
                ylabel('Cell #'); 
                title('Smoothened PF map'); colorbar;
            subplot(242); histogram(hist.spkPeak_SIspk,100,'Normalization','probability'); 
                hold on;
                plot(mean(hist.spkPeak_SIspk),0.06,'kv','markerfacecolor','k','markersize',6); hold off;
                xlabel('PSD peak'); set(gca,'xscale','log');
                ylabel('Cells (%)');
            subplot(243); histogram(hist.spkMean_SIspk,40,'Normalization','probability');
                hold on;
                plot(mean(hist.spkMean_SIspk),0.08,'kv','markerfacecolor','k','markersize',6); hold off;
                xlabel('PSD mean');
                ylabel('Cells (%)');
            subplot(244); histogram(hist.infoMap_SIspk,40,'Normalization','probability');
                hold on;
                plot(mean(hist.infoMap_SIspk),0.12,'kv','markerfacecolor','k','markersize',6); hold off;
                xlabel('Mutual info (bits)');
                ylabel('Cells (%)');
            subplot(246); histogram(hist.fieldsize_SIspk(:,e),40,'Normalization','probability');
                hold on;
                mean_fieldsize_SIspk = mean(hist.fieldsize_SIspk(:,e));
                plot(mean_fieldsize_SIspk,0.1,'kv','markerfacecolor','k','markersize',6); hold off;
                xlabel('Field size(cm)');
                ylabel('Cells (%)');
            subplot(247); histogram(hist.centroid_SIspk(:,e),30);
                hold on;
                mean_centroid_SIspk = mean(hist.centroid_SIspk(:,e));
                plot(mean_centroid_SIspk,10,'kv','markerfacecolor','k','markersize',6); hold off;
                xlabel('Position (cm)');
                ylabel('No. of PCs');
            subplot(248); histogram(PFdata.bin_phi(:,e),30,'Normalization','probability');
                hold on;
                mean_bin_phi = mean(PFdata.bin_phi(:,e));
                plot(mean_bin_phi,0.04,'kv','markerfacecolor','k','markersize',6); hold off;
                xlabel('Position (cm)');
                ylabel('Time (%)');
            
            if fsave
                if ~exist(sdir,'dir'), mkdir(sdir); end

                if Nepochs == 1
                    fname_fig = [sdir 'hist_SIspk/' fname_pref  '_populSummary'];
                else
                    fname_fig = [sdir 'hist_SIspk/' fname_pref  '_populSummary_' num2str(e) 'of' num2str(Nepochs) 'ep'];
                end
                savefig( fh, fname_fig );
                saveas( fh, fname_fig, 'png' );
                if fclose, close( fh ); end
            end
        end
        
        % remapping within a session
        if Nepochs > 1
            fh = figure;
            for ei = 1:Nepochs % rows: sorting
                for ej = 1:Nepochs % cols: epochs 
                    subplot(Nepochs, Nepochs, (ei-1)*Nepochs + ej); imagesc(hist.normpfMap_SIspk(hist.sortIdx_SIspk(:,ei),:,ej)); 
                    title(['Epoch ' num2str(ej)]); ylabel(['Epoch' num2str(ei) ' sorting']);
                end
            end
            if fsave
                fname_fig = [sdir 'hist_SIspk/' fname_pref '_remapping_hist_SIspk'];
                savefig( fh, fname_fig );
                saveas( fh, fname_fig, 'png' );
                if fclose, close( fh ); end
            end
        end
    end
    
    % spike raster plots for NON-place cells 
    Nnonpcs = numel(hist.nonpcIdx_SIspk);
    if Nnonpcs > 0 
        [nRow, nCol] = getnRownCol(Nnonpcs);
        nPlot = nRow*nCol;
        Nfig = round((Nnonpcs/nPlot))-1;
        if Nfig<0, Nfig = 0; end

        for ii=0:Nfig
            fh = figure;
            ha = tight_subplot(nRow,nCol,[.02 .015],[.01 .05],[.01 .01]);
            for jj=0:nPlot-1
                if (ii*nPlot+jj+1) <= Nnonpcs
                    axes(ha(+jj+1));
                    imagesc(hist.normspkRaster_SIspk_nonpc{ii*nPlot+jj+1}); 
                    yticks(ytick_files); yticklabels(ytick_files); % ylabel('Trial #');
                    xticks([]); % xticklabels([1 50 100]); xlabel('Position (cm)');
                    title(['Non-PC ' num2str(ii*nPlot+jj+1)],'fontsize',15);
                end
            end
            if fsave
                if Nnonpcs/nPlot <= 1
                    fname_fig = [sdir 'hist_SIspk/' fname_pref '_normspkRaster_nonPCs_hist_SIspk'];
                else
                    fname_fig = [sdir 'hist_SIspk/' fname_pref '_normspkRaster_nonPCs_hist_SIspk_' num2str(ii+1)];
                end
                savefig( fh, fname_fig );
                saveas( fh, fname_fig, 'png' );
                if fclose, close( fh ); end
            end
        end 
        
        for ii=0:Nfig
            fh = figure;
            ha = tight_subplot(nRow,nCol,[.02 .015],[.01 .05],[.01 .01]);
            for jj=0:nPlot-1
                if (ii*nPlot+jj+1) <= Nnonpcs
                    axes(ha(+jj+1));
                    imagesc(hist.spkRaster_SIspk_nonpc{ii*nPlot+jj+1}); 
                    yticks(ytick_files); yticklabels(ytick_files); % ylabel('Trial #');
                    xticks([]); % xticklabels([1 50 100]); xlabel('Position (cm)');
                    title(['Non-PC ' num2str(ii*nPlot+jj+1)],'fontsize',15);
                    colorbar;
                end
            end
            if fsave
                if Nnonpcs/nPlot <= 1
                    fname_fig = [sdir 'hist_SIspk/' fname_pref '_spkRaster_nonPCs_hist_SIspk'];
                else
                    fname_fig = [sdir 'hist_SIspk/' fname_pref '_spkRaster_nonPCs_hist_SIspk_' num2str(ii+1)];
                end
                savefig( fh, fname_fig );
                saveas( fh, fname_fig, 'png' );
                if fclose, close( fh ); end
            end
        end
    end
        
    
    %% asd, MI (bits)
    Npcs = numel(asd.pcIdx_MI);

    if Npcs > 0
        % spike raster plots for place cells
        [nRow, nCol] = getnRownCol(Npcs);
        nPlot = nRow*nCol;
        Nfig = round((Npcs/nPlot))-1;
        if Nfig<0, Nfig = 0; end
        if ~exist([sdir 'asd_MI/'],'dir'), mkdir([sdir 'asd_MI/']); end

        for ii=0:Nfig
            fh = figure;
            ha = tight_subplot(nRow,nCol,[.02 .015],[.01 .05],[.01 .01]);
            for jj=0:nPlot-1
                if (ii*nPlot+jj+1) <= Npcs
                    axes(ha(+jj+1));
                    imagesc(asd.normspkRaster_MI_pc{ii*nPlot+jj+1}); 
                    yticks(ytick_files); yticklabels(ytick_files); % ylabel('Trial #');
                    xticks([]); % xticklabels([1 50 100]); xlabel('Position (cm)');
                    title(['PC ' num2str(ii*nPlot+jj+1)],'fontsize',15);
    %                 hold on
    %                 yyaxis right; plot(asd.meanspkRaster_MI(ii*nPlot+jj+1,:),'w-'); hold off
                end
            end
            if fsave
                if Npcs/nPlot <= 1
                    fname_fig = [sdir 'asd_MI/' fname_pref '_normspkRaster_PCs_asd_MI'];
                else
                    fname_fig = [sdir 'asd_MI/' fname_pref '_normspkRaster_PCs_asd_MI_' num2str(ii+1)];
                end
                savefig( fh, fname_fig );
                saveas( fh, fname_fig, 'png' );
                if fclose, close( fh ); end
            end
        end 

        for ii=0:Nfig
            fh = figure;
            ha = tight_subplot(nRow,nCol,[.02 .015],[.01 .05],[.01 .01]);
            for jj=0:nPlot-1
                if (ii*nPlot+jj+1) <= Npcs
                    axes(ha(+jj+1));
                    imagesc(asd.spkRaster_MI_pc{ii*nPlot+jj+1}); 
                    yticks(ytick_files); yticklabels(ytick_files); % ylabel('Trial #');
                    xticks([]); % xticklabels([1 50 100]); xlabel('Position (cm)');
                    title(['PC ' num2str(ii*nPlot+jj+1)],'fontsize',15);
                    colorbar;
    %                 hold on
    %                 yyaxis right; plot(asd.meanspkRaster_MI(ii*nPlot+jj+1,:),'w-'); hold off
                end
            end
            if fsave
                if Npcs/nPlot <= 1
                    fname_fig = [sdir 'asd_MI/' fname_pref '_spkRaster_PCs_asd_MI'];
                else
                    fname_fig = [sdir 'asd_MI/' fname_pref '_spkRaster_PCs_asd_MI_' num2str(ii+1)];
                end
                savefig( fh, fname_fig );
                saveas( fh, fname_fig, 'png' );
                if fclose, close( fh ); end
            end
        end
        
        % summary of population data
        for e = 1:Nepochs
            fh = figure;
            subplot(241); imagesc(asd.sort_normpfMap_MI(:,:,e));
                xticks([1 15 30]); xticklabels([1 50 100]); xlabel('Position (cm)');
                if Npcs>1 
                    yticks([1 Npcs]); yticklabels([Npcs 1]);
                else 
                    yticks(1); yticklabels(1); 
                end
                ylabel('Cell #'); 
                title('PF map'); colorbar;
            subplot(242); histogram(asd.spkPeak_MI,100,'Normalization','probability'); 
                hold on;
                plot(mean(asd.spkPeak_MI),0.06,'kv','markerfacecolor','k','markersize',6); hold off;
                xlabel('PSD peak'); set(gca,'xscale','log');
                ylabel('Cells (%)');
            subplot(243); histogram(asd.spkMean_MI,40,'Normalization','probability');
                hold on;
                plot(mean(asd.spkMean_MI),0.08,'kv','markerfacecolor','k','markersize',6); hold off;
                xlabel('PSD mean');
                ylabel('Cells (%)');
            subplot(244); histogram(asd.infoMap_MI,40,'Normalization','probability');
                hold on;
                plot(mean(asd.infoMap_MI),0.12,'kv','markerfacecolor','k','markersize',6); hold off;
                xlabel('Mutual info (bits)');
                ylabel('Cells (%)');
            subplot(246); histogram(asd.fieldsize_MI(:,e),40,'Normalization','probability');
                hold on;
                mean_fieldsize_MI = mean(asd.fieldsize_MI(:,e));
                plot(mean_fieldsize_MI,0.1,'kv','markerfacecolor','k','markersize',6); hold off;
                xlabel('Field size(cm)');
                ylabel('Cells (%)');
            subplot(247); histogram(asd.centroid_MI(:,e),30);
                hold on;
                mean_centroid_MI = mean(asd.centroid_MI(:,e));
                plot(mean_centroid_MI,10,'kv','markerfacecolor','k','markersize',6); hold off;
                xlabel('Position (cm)');
                ylabel('No. of PCs');
            subplot(248); histogram(PFdata.bin_phi(:,e),30,'Normalization','probability');
                hold on;
                mean_bin_phi = mean(PFdata.bin_phi(:,e));
                plot(mean_bin_phi,0.04,'kv','markerfacecolor','k','markersize',6); hold off;
                xlabel('Position (cm)');
                ylabel('Time (%)');
            
            if fsave
                if ~exist(sdir,'dir'), mkdir(sdir); end

                if Nepochs == 1
                    fname_fig = [sdir 'asd_MI/' fname_pref  '_populSummary'];
                else
                    fname_fig = [sdir 'asd_MI/' fname_pref  '_populSummary_' num2str(e) 'of' num2str(Nepochs) 'ep'];
                end
                savefig( fh, fname_fig );
                saveas( fh, fname_fig, 'png' );
                if fclose, close( fh ); end
            end
        end
        
        % remapping within a session
        if Nepochs > 1
            fh = figure;
            for ei = 1:Nepochs % rows: sorting
                for ej = 1:Nepochs % cols: epochs 
                    subplot(Nepochs, Nepochs, (ei-1)*Nepochs + ej); imagesc(asd.normpfMap_MI(asd.sortIdx_MI(:,ei),:,ej)); 
                    title(['Epoch ' num2str(ej)]); ylabel(['Epoch' num2str(ei) ' sorting']);
                end
            end
            if fsave
                fname_fig = [sdir 'asd_MI/' fname_pref '_remapping_asd_MI'];
                savefig( fh, fname_fig );
                saveas( fh, fname_fig, 'png' );
                if fclose, close( fh ); end
            end
        end
    end
    
    % spike raster plots for NON-place cells 
    Nnonpcs = numel(asd.nonpcIdx_MI);
    if Nnonpcs > 0
        [nRow, nCol] = getnRownCol(Nnonpcs);
        nPlot = nRow*nCol;
        Nfig = round((Nnonpcs/nPlot))-1;
        if Nfig<0, Nfig = 0; end

        for ii=0:Nfig
            fh = figure;
            ha = tight_subplot(nRow,nCol,[.02 .015],[.01 .05],[.01 .01]);
            for jj=0:nPlot-1
                if (ii*nPlot+jj+1) <= Nnonpcs
                    axes(ha(+jj+1));
                    imagesc(asd.normspkRaster_MI_nonpc{ii*nPlot+jj+1}); 
                    yticks(ytick_files); yticklabels(ytick_files); % ylabel('Trial #');
                    xticks([]); % xticklabels([1 50 100]); xlabel('Position (cm)');
                    title(['Non-PC ' num2str(ii*nPlot+jj+1)],'fontsize',15);
                end
            end
            if fsave
                if Nnonpcs/nPlot <= 1
                    fname_fig = [sdir 'asd_MI/' fname_pref '_normspkRaster_nonPCs_asd_MI'];
                else
                    fname_fig = [sdir 'asd_MI/' fname_pref '_normspkRaster_nonPCs_asd_MI_' num2str(ii+1)];
                end
                savefig( fh, fname_fig );
                saveas( fh, fname_fig, 'png' );
                if fclose, close( fh ); end
            end
        end 
    
        for ii=0:Nfig
            fh = figure;
            ha = tight_subplot(nRow,nCol,[.02 .015],[.01 .05],[.01 .01]);
            for jj=0:nPlot-1
                if (ii*nPlot+jj+1) <= Nnonpcs
                    axes(ha(+jj+1));
                    imagesc(asd.spkRaster_MI_nonpc{ii*nPlot+jj+1}); 
                    yticks(ytick_files); yticklabels(ytick_files); % ylabel('Trial #');
                    xticks([]); % xticklabels([1 50 100]); xlabel('Position (cm)');
                    title(['Non-PC ' num2str(ii*nPlot+jj+1)],'fontsize',15);
                    colorbar;
                end
            end
            if fsave
                if Nnonpcs/nPlot <= 1
                    fname_fig = [sdir 'asd_MI/' fname_pref '_spkRaster_nonPCs_asd_MI'];
                else
                    fname_fig = [sdir 'asd_MI/' fname_pref '_spkRaster_nonPCs_asd_MI_' num2str(ii+1)];
                end
                savefig( fh, fname_fig );
                saveas( fh, fname_fig, 'png' );
                if fclose, close( fh ); end
            end
        end
    end
    
    %% asd, SI (bits/sec)
    Npcs = numel(asd.pcIdx_SIsec);

    if Npcs > 0
        % spike raster plots for place cells
        [nRow, nCol] = getnRownCol(Npcs);
        nPlot = nRow*nCol;
        Nfig = round((Npcs/nPlot))-1;
        if Nfig<0, Nfig = 0; end
        if ~exist([sdir 'asd_SIsec/'],'dir'), mkdir([sdir 'asd_SIsec/']); end

        for ii=0:Nfig
            fh = figure;
            ha = tight_subplot(nRow,nCol,[.02 .015],[.01 .05],[.01 .01]);
            for jj=0:nPlot-1
                if (ii*nPlot+jj+1) <= Npcs
                    axes(ha(+jj+1));
                    imagesc(asd.normspkRaster_SIsec_pc{ii*nPlot+jj+1}); 
                    yticks(ytick_files); yticklabels(ytick_files); % ylabel('Trial #');
                    xticks([]); % xticklabels([1 50 100]); xlabel('Position (cm)');
                    title(['PC ' num2str(ii*nPlot+jj+1)],'fontsize',15);
    %                 hold on
    %                 yyaxis right; plot(asd.meanspkRaster_SIsec(ii*nPlot+jj+1,:),'w-'); hold off
                end
            end
            if fsave
                if Npcs/nPlot <= 1
                    fname_fig = [sdir 'asd_SIsec/' fname_pref '_normspkRaster_PCs_asd_SIsec'];
                else
                    fname_fig = [sdir 'asd_SIsec/' fname_pref '_normspkRaster_PCs_asd_SIsec_' num2str(ii+1)];
                end
                savefig( fh, fname_fig );
                saveas( fh, fname_fig, 'png' );
                if fclose, close( fh ); end
            end
        end 
        
        for ii=0:Nfig
            fh = figure;
            ha = tight_subplot(nRow,nCol,[.02 .015],[.01 .05],[.01 .01]);
            for jj=0:nPlot-1
                if (ii*nPlot+jj+1) <= Npcs
                    axes(ha(+jj+1));
                    imagesc(asd.spkRaster_SIsec_pc{ii*nPlot+jj+1}); 
                    yticks(ytick_files); yticklabels(ytick_files); % ylabel('Trial #');
                    xticks([]); % xticklabels([1 50 100]); xlabel('Position (cm)');
                    title(['PC ' num2str(ii*nPlot+jj+1)],'fontsize',15);
                    colorbar;
    %                 hold on
    %                 yyaxis right; plot(asd.meanspkRaster_SIsec(ii*nPlot+jj+1,:),'w-'); hold off
                end
            end
            if fsave
                if Npcs/nPlot <= 1
                    fname_fig = [sdir 'asd_SIsec/' fname_pref '_spkRaster_PCs_asd_SIsec'];
                else
                    fname_fig = [sdir 'asd_SIsec/' fname_pref '_spkRaster_PCs_asd_SIsec_' num2str(ii+1)];
                end
                savefig( fh, fname_fig );
                saveas( fh, fname_fig, 'png' );
                if fclose, close( fh ); end
            end
        end

        % summary of population data
        for e = 1:Nepochs
            fh = figure;
            subplot(241); imagesc(asd.sort_normpfMap_SIsec(:,:,e));
                xticks([1 15 30]); xticklabels([1 50 100]); xlabel('Position (cm)');
                if Npcs>1 
                    yticks([1 Npcs]); yticklabels([Npcs 1]);
                else 
                    yticks(1); yticklabels(1); 
                end
                ylabel('Cell #'); 
                title('PF map'); colorbar;
            subplot(242); histogram(asd.spkPeak_SIsec,100,'Normalization','probability'); 
                hold on;
                plot(mean(asd.spkPeak_SIsec),0.06,'kv','markerfacecolor','k','markersize',6); hold off;
                xlabel('PSD peak'); set(gca,'xscale','log');
                ylabel('Cells (%)');
            subplot(243); histogram(asd.spkMean_SIsec,40,'Normalization','probability');
                hold on;
                plot(mean(asd.spkMean_SIsec),0.08,'kv','markerfacecolor','k','markersize',6); hold off;
                xlabel('PSD mean');
                ylabel('Cells (%)');
            subplot(244); histogram(asd.infoMap_SIsec,40,'Normalization','probability');
                hold on;
                plot(mean(asd.infoMap_SIsec),0.12,'kv','markerfacecolor','k','markersize',6); hold off;
                xlabel('Mutual info (bits)');
                ylabel('Cells (%)');
            subplot(246); histogram(asd.fieldsize_SIsec(:,e),40,'Normalization','probability');
                hold on;
                mean_fieldsize_SIsec = mean(asd.fieldsize_SIsec(:,e));
                plot(mean_fieldsize_SIsec,0.1,'kv','markerfacecolor','k','markersize',6); hold off;
                xlabel('Field size(cm)');
                ylabel('Cells (%)');
            subplot(247); histogram(asd.centroid_SIsec(:,e),30);
                hold on;
                mean_centroid_SIsec = mean(asd.centroid_SIsec(:,e));
                plot(mean_centroid_SIsec,10,'kv','markerfacecolor','k','markersize',6); hold off;
                xlabel('Position (cm)');
                ylabel('No. of PCs');
            subplot(248); histogram(PFdata.bin_phi(:,e),30,'Normalization','probability');
                hold on;
                mean_bin_phi = mean(PFdata.bin_phi(:,e));
                plot(mean_bin_phi,0.04,'kv','markerfacecolor','k','markersize',6); hold off;
                xlabel('Position (cm)');
                ylabel('Time (%)');
            
            if fsave
                if ~exist(sdir,'dir'), mkdir(sdir); end

                if Nepochs == 1
                    fname_fig = [sdir 'asd_SIsec/' fname_pref  '_populSummary'];
                else
                    fname_fig = [sdir 'asd_SIsec/' fname_pref  '_populSummary_' num2str(e) 'of' num2str(Nepochs) 'ep'];
                end
                savefig( fh, fname_fig );
                saveas( fh, fname_fig, 'png' );
                if fclose, close( fh ); end
            end
        end
        
        % remapping within a session
        if Nepochs > 1
            fh = figure;
            for ei = 1:Nepochs % rows: sorting
                for ej = 1:Nepochs % cols: epochs 
                    subplot(Nepochs, Nepochs, (ei-1)*Nepochs + ej); imagesc(asd.normpfMap_SIsec(asd.sortIdx_SIsec(:,ei),:,ej)); 
                    title(['Epoch ' num2str(ej)]); ylabel(['Epoch' num2str(ei) ' sorting']);
                end
            end
            if fsave
                fname_fig = [sdir 'asd_SIsec/' fname_pref '_remapping_asd_SIsec'];
                savefig( fh, fname_fig );
                saveas( fh, fname_fig, 'png' );
                if fclose, close( fh ); end
            end
        end
    end
    
    % spike raster plots for NON-place cells 
    Nnonpcs = numel(asd.nonpcIdx_SIsec);
    if Nnonpcs > 0 
        [nRow, nCol] = getnRownCol(Nnonpcs);
        nPlot = nRow*nCol;
        Nfig = round((Nnonpcs/nPlot))-1;
        if Nfig<0, Nfig = 0; end

        for ii=0:Nfig
            fh = figure;
            ha = tight_subplot(nRow,nCol,[.02 .015],[.01 .05],[.01 .01]);
            for jj=0:nPlot-1
                if (ii*nPlot+jj+1) <= Nnonpcs
                    axes(ha(+jj+1));
                    imagesc(asd.normspkRaster_SIsec_nonpc{ii*nPlot+jj+1}); 
                    yticks(ytick_files); yticklabels(ytick_files); % ylabel('Trial #');
                    xticks([]); % xticklabels([1 50 100]); xlabel('Position (cm)');
                    title(['Non-PC ' num2str(ii*nPlot+jj+1)],'fontsize',15);
                end
            end
            if fsave
                if Nnonpcs/nPlot <= 1
                    fname_fig = [sdir 'asd_SIsec/' fname_pref '_normspkRaster_nonPCs_asd_SIsec'];
                else
                    fname_fig = [sdir 'asd_SIsec/' fname_pref '_normspkRaster_nonPCs_asd_SIsec_' num2str(ii+1)];
                end
                savefig( fh, fname_fig );
                saveas( fh, fname_fig, 'png' );
                if fclose, close( fh ); end
            end
        end 
        
        for ii=0:Nfig
            fh = figure;
            ha = tight_subplot(nRow,nCol,[.02 .015],[.01 .05],[.01 .01]);
            for jj=0:nPlot-1
                if (ii*nPlot+jj+1) <= Nnonpcs
                    axes(ha(+jj+1));
                    imagesc(asd.spkRaster_SIsec_nonpc{ii*nPlot+jj+1}); 
                    yticks(ytick_files); yticklabels(ytick_files); % ylabel('Trial #');
                    xticks([]); % xticklabels([1 50 100]); xlabel('Position (cm)');
                    title(['Non-PC ' num2str(ii*nPlot+jj+1)],'fontsize',15);
                    colorbar;
                end
            end
            if fsave
                if Nnonpcs/nPlot <= 1
                    fname_fig = [sdir 'asd_SIsec/' fname_pref '_spkRaster_nonPCs_asd_SIsec'];
                else
                    fname_fig = [sdir 'asd_SIsec/' fname_pref '_spkRaster_nonPCs_asd_SIsec_' num2str(ii+1)];
                end
                savefig( fh, fname_fig );
                saveas( fh, fname_fig, 'png' );
                if fclose, close( fh ); end
            end
        end 
    end
        
    %% asd, SI (bits/spike)
    Npcs = numel(asd.pcIdx_SIspk);

    if Npcs > 0
        % spike raster plots for place cells
        [nRow, nCol] = getnRownCol(Npcs);
        nPlot = nRow*nCol;
        Nfig = round((Npcs/nPlot))-1;
        if Nfig<0, Nfig = 0; end
        if ~exist([sdir 'asd_SIspk/'],'dir'), mkdir([sdir 'asd_SIspk/']); end

        for ii=0:Nfig
            fh = figure;
            ha = tight_subplot(nRow,nCol,[.02 .015],[.01 .05],[.01 .01]);
            for jj=0:nPlot-1
                if (ii*nPlot+jj+1) <= Npcs
                    axes(ha(+jj+1));
                    imagesc(asd.normspkRaster_SIspk_pc{ii*nPlot+jj+1}); 
                    yticks(ytick_files); yticklabels(ytick_files); % ylabel('Trial #');
                    xticks([]); % xticklabels([1 50 100]); xlabel('Position (cm)');
                    title(['PC ' num2str(ii*nPlot+jj+1)],'fontsize',15);
    %                 hold on
    %                 yyaxis right; plot(asd.meanspkRaster_SIspk(ii*nPlot+jj+1,:),'w-'); hold off
                end
            end
            if fsave
                if Npcs/nPlot <= 1
                    fname_fig = [sdir 'asd_SIspk/' fname_pref '_normspkRaster_PCs_asd_SIspk'];
                else
                    fname_fig = [sdir 'asd_SIspk/' fname_pref '_normspkRaster_PCs_asd_SIspk_' num2str(ii+1)];
                end
                savefig( fh, fname_fig );
                saveas( fh, fname_fig, 'png' );
                if fclose, close( fh ); end
            end
        end 

        for ii=0:Nfig
            fh = figure;
            ha = tight_subplot(nRow,nCol,[.02 .015],[.01 .05],[.01 .01]);
            for jj=0:nPlot-1
                if (ii*nPlot+jj+1) <= Npcs
                    axes(ha(+jj+1));
                    imagesc(asd.spkRaster_SIspk_pc{ii*nPlot+jj+1}); 
                    yticks(ytick_files); yticklabels(ytick_files); % ylabel('Trial #');
                    xticks([]); % xticklabels([1 50 100]); xlabel('Position (cm)');
                    title(['PC ' num2str(ii*nPlot+jj+1)],'fontsize',15);
                    colorbar;
    %                 hold on
    %                 yyaxis right; plot(asd.meanspkRaster_SIspk(ii*nPlot+jj+1,:),'w-'); hold off
                end
            end
            if fsave
                if Npcs/nPlot <= 1
                    fname_fig = [sdir 'asd_SIspk/' fname_pref '_spkRaster_PCs_asd_SIspk'];
                else
                    fname_fig = [sdir 'asd_SIspk/' fname_pref '_spkRaster_PCs_asd_SIspk_' num2str(ii+1)];
                end
                savefig( fh, fname_fig );
                saveas( fh, fname_fig, 'png' );
                if fclose, close( fh ); end
            end
        end
        
        % summary of population data
        for e = 1:Nepochs
            fh = figure;
            subplot(241); imagesc(asd.sort_normpfMap_SIspk(:,:,e));
                xticks([1 15 30]); xticklabels([1 50 100]); xlabel('Position (cm)');
                if Npcs>1 
                    yticks([1 Npcs]); yticklabels([Npcs 1]);
                else 
                    yticks(1); yticklabels(1); 
                end
                ylabel('Cell #'); 
                title('PF map'); colorbar;
            subplot(242); histogram(asd.spkPeak_SIspk,100,'Normalization','probability'); 
                hold on;
                plot(mean(asd.spkPeak_SIspk),0.06,'kv','markerfacecolor','k','markersize',6); hold off;
                xlabel('PSD peak'); set(gca,'xscale','log');
                ylabel('Cells (%)');
            subplot(243); histogram(asd.spkMean_SIspk,40,'Normalization','probability');
                hold on;
                plot(mean(asd.spkMean_SIspk),0.08,'kv','markerfacecolor','k','markersize',6); hold off;
                xlabel('PSD mean');
                ylabel('Cells (%)');
            subplot(244); histogram(asd.infoMap_SIspk,40,'Normalization','probability');
                hold on;
                plot(mean(asd.infoMap_SIspk),0.12,'kv','markerfacecolor','k','markersize',6); hold off;
                xlabel('Mutual info (bits)');
                ylabel('Cells (%)');
            subplot(246); histogram(asd.fieldsize_SIspk(:,e),40,'Normalization','probability');
                hold on;
                mean_fieldsize_SIspk = mean(asd.fieldsize_SIspk(:,e));
                plot(mean_fieldsize_SIspk,0.1,'kv','markerfacecolor','k','markersize',6); hold off;
                xlabel('Field size(cm)');
                ylabel('Cells (%)');
            subplot(247); histogram(asd.centroid_SIspk(:,e),30);
                hold on;
                mean_centroid_SIspk = mean(asd.centroid_SIspk(:,e));
                plot(mean_centroid_SIspk,10,'kv','markerfacecolor','k','markersize',6); hold off;
                xlabel('Position (cm)');
                ylabel('No. of PCs');
            subplot(248); histogram(PFdata.bin_phi(:,e),30,'Normalization','probability');
                hold on;
                mean_bin_phi = mean(PFdata.bin_phi(:,e));
                plot(mean_bin_phi,0.04,'kv','markerfacecolor','k','markersize',6); hold off;
                xlabel('Position (cm)');
                ylabel('Time (%)');
            
            if fsave
                if ~exist(sdir,'dir'), mkdir(sdir); end

                if Nepochs == 1
                    fname_fig = [sdir 'asd_SIspk/' fname_pref  '_populSummary'];
                else
                    fname_fig = [sdir 'asd_SIspk/' fname_pref  '_populSummary_' num2str(e) 'of' num2str(Nepochs) 'ep'];
                end
                savefig( fh, fname_fig );
                saveas( fh, fname_fig, 'png' );
                if fclose, close( fh ); end
            end
        end
        
        % remapping within a session
        if Nepochs > 1
            fh = figure;
            for ei = 1:Nepochs % rows: sorting
                for ej = 1:Nepochs % cols: epochs 
                    subplot(Nepochs, Nepochs, (ei-1)*Nepochs + ej); imagesc(asd.normpfMap_SIspk(asd.sortIdx_SIspk(:,ei),:,ej)); 
                    title(['Epoch ' num2str(ej)]); ylabel(['Epoch' num2str(ei) ' sorting']);
                end
            end
            if fsave
                fname_fig = [sdir 'asd_SIspk/' fname_pref '_remapping_asd_SIspk'];
                savefig( fh, fname_fig );
                saveas( fh, fname_fig, 'png' );
                if fclose, close( fh ); end
            end
        end
    end
    
    % spike raster plots for NON-place cells 
    Nnonpcs = numel(asd.nonpcIdx_SIspk);
    if Nnonpcs > 0
        [nRow, nCol] = getnRownCol(Nnonpcs);
        nPlot = nRow*nCol;
        Nfig = round((Nnonpcs/nPlot))-1;
        if Nfig<0, Nfig = 0; end

        for ii=0:Nfig
            fh = figure;
            ha = tight_subplot(nRow,nCol,[.02 .015],[.01 .05],[.01 .01]);
            for jj=0:nPlot-1
                if (ii*nPlot+jj+1) <= Nnonpcs
                    axes(ha(+jj+1));
                    imagesc(asd.normspkRaster_SIspk_nonpc{ii*nPlot+jj+1}); 
                    yticks(ytick_files); yticklabels(ytick_files); % ylabel('Trial #');
                    xticks([]); % xticklabels([1 50 100]); xlabel('Position (cm)');
                    title(['Non-PC ' num2str(ii*nPlot+jj+1)],'fontsize',15);
                end
            end
            if fsave
                if Nnonpcs/nPlot <= 1
                    fname_fig = [sdir 'asd_SIspk/' fname_pref '_normspkRaster_nonPCs_asd_SIspk'];
                else
                    fname_fig = [sdir 'asd_SIspk/' fname_pref '_normspkRaster_nonPCs_asd_SIspk_' num2str(ii+1)];
                end
                savefig( fh, fname_fig );
                saveas( fh, fname_fig, 'png' );
                if fclose, close( fh ); end
            end
        end 
    
        for ii=0:Nfig
            fh = figure;
            ha = tight_subplot(nRow,nCol,[.02 .015],[.01 .05],[.01 .01]);
            for jj=0:nPlot-1
                if (ii*nPlot+jj+1) <= Nnonpcs
                    axes(ha(+jj+1));
                    imagesc(asd.spkRaster_SIspk_nonpc{ii*nPlot+jj+1}); 
                    yticks(ytick_files); yticklabels(ytick_files); % ylabel('Trial #');
                    xticks([]); % xticklabels([1 50 100]); xlabel('Position (cm)');
                    title(['Non-PC ' num2str(ii*nPlot+jj+1)],'fontsize',15);
                    colorbar;
                end
            end
            if fsave
                if Nnonpcs/nPlot <= 1
                    fname_fig = [sdir 'asd_SIspk/' fname_pref '_spkRaster_nonPCs_asd_SIspk'];
                else
                    fname_fig = [sdir 'asd_SIspk/' fname_pref '_spkRaster_nonPCs_asd_SIspk_' num2str(ii+1)];
                end
                savefig( fh, fname_fig );
                saveas( fh, fname_fig, 'png' );
                if fclose, close( fh ); end
            end
        end 
    end
        
end