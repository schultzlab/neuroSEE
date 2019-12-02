function plotPF_1d(list, occMap, hist, asd, normspkMap_pertrial, ytick_files, fsave, sdir, files, fclose)
    if nargin < 10, fclose = false; end
    if nargin < 9, fsave = false; end
    if nargin < 8, fsave = false; end
    if nargin < 7, fsave = false; end
    
    Ncells = numel(normspkMap_pertrial);
    Npcs = length(hist.pcIdx);
    Npcs_asd = length(asd.pcIdx);
    Nepochs = size(occMap,3);

    % MouseID and experiment name
    [ mouseid, expname ] = find_mouseIDexpname(list);

    % summary of occMap, spkMaps, pfMaps
    for e = 1:Nepochs
        fh = figure('Position',[1087 648 800 800]);
        subplot(10,8,2:4); imagesc(occMap(e,:));
            xticks([]); yticks([]); ylabel('Occ');
            title('Histogram estimation'); colorbar;
        
        if Npcs > 0     
            subplot(10,8,[10:12,18:20,26:28]);
                imagesc(hist.spkMap(hist.sortIdx(:,e),:,e));
                xticks([]);
                if Npcs >  1
                    yticks([1 Npcs]); 
                else
                    yticks(1);
                end
                ylabel('Cell #'); 
                title('Spike map'); colorbar;
            subplot(10,8,[33,41,49]); 
                imagesc(hist.infoMap(hist.sortIdx,1,e));
                xticks([]);
                if Npcs >  1
                    yticks([1 Npcs]); 
                else
                    yticks(1);
                end
                ylabel('Cell #');
                title('Max MI'); colorbar;
            subplot(10,8,[34:36,42:44,50:52]);    
                imagesc(hist.sort_pfMap(:,:,e)); 
                xticks([]); 
                if Npcs >  1
                    yticks([1 Npcs]); 
                else
                    yticks(1);
                end
                title('Place field map'); colorbar;
            subplot(10,8,[58:60,66:68,74:76]);    
                imagesc(hist.sort_pfMap_sm(:,:,e)); 
                if Npcs >  1
                    yticks([1 Npcs]); 
                else
                    yticks(1);
                end
                ylabel('Cell #');
                xticks([1 15 30]); xticklabels([1 50 100]);
                xlabel('Position (cm)');
                title('Smoothened pf map'); colorbar;
        end
            
        if Npcs_asd > 0    
            subplot(10,8,6:8); imagesc(occMap(e,:));
                xticks([]); yticks([]); ylabel('Occ');
                title('ASD'); colorbar;
            subplot(10,8,[14:16,22:24,30:32]);
                imagesc(asd.spkMap(asd.sortIdx(:,e),:,e));
                xticks([]);
                if Npcs_asd >  1
                    yticks([1 Npcs_asd]); 
                else
                    yticks(1);
                end
                ylabel('Cell #');
                title('Spike map'); colorbar;
            subplot(10,8,[37,45,53]); imagesc(asd.infoMap(asd.sortIdx,1,e));
                xticks([]);
                if Npcs_asd >  1
                    yticks([1 Npcs_asd]); 
                else
                    yticks(1);
                end 
                title('Max MI'); colorbar;
            subplot(10,8,[38:40,46:48,54:56]);    
                imagesc(asd.sort_pfMap(:,:,e)); 
                if Npcs_asd >  1
                    yticks([1 Npcs_asd]); 
                else
                    yticks(1);
                end
                xticks([1 15 30]); xticklabels([1 50 100]);
                xlabel('Position (cm)');
                title('Place field map'); colorbar;
        else
            subplot(10,8,6:8); title('No place cells for ASD method');
        end
        if fsave
            if ~exist([sdir 'PFmaps/'],'dir'), mkdir([sdir 'PFmaps/']); end
        
            if Nepochs == 1
                fname_fig = [sdir 'PFmaps/' mouseid '_' expname '_ref_' files(1,:)  '_PFmaps'];
            else
                fname_fig = [sdir 'PFmaps/' mouseid '_' expname '_ref_' files(1,:)  '_PFmaps_' num2str(e) 'of' num2str(Nepochs) 'ep'];
            end
            savefig( fh, fname_fig );
            saveas( fh, fname_fig, 'png' );
            if fclose, close( fh ); end
        end
    end

    % summary of occMap, normspkMaps, normpfMaps
    for e = 1:Nepochs
        fh = figure('Position',[1087 648 800 800]);
        subplot(10,8,2:4); imagesc(occMap(e,:));
            xticks([]); yticks([]); ylabel('Occ');
            title('Histogram estimation'); colorbar;
            
        if Npcs > 0
            subplot(10,8,[10:12,18:20,26:28]);
                imagesc(hist.normspkMap(hist.sortIdx(:,e),:,e));
                xticks([]);
                if Npcs >  1
                    yticks([1 Npcs]); 
                else
                    yticks(1);
                end
                ylabel('Cell #');
                title('Normalised spk map'); colorbar;
            subplot(10,8,[33,41,49]); 
                imagesc(hist.infoMap(hist.sortIdx,1,e));
                xticks([]);
                if Npcs >  1
                    yticks([1 Npcs]); 
                else
                    yticks(1);
                end
                ylabel('Cell #');
                title('Max MI'); colorbar;
            subplot(10,8,[34:36,42:44,50:52]);    
                imagesc(hist.sort_normpfMap(:,:,e)); 
                xticks([]);
                if Npcs >  1
                    yticks([1 Npcs]); 
                else
                    yticks(1);
                end
                title('Normalised pf map'); colorbar;
            subplot(10,8,[58:60,66:68,74:76]);    
                imagesc(hist.sort_normpfMap_sm(:,:,e)); 
                if Npcs >  1
                    yticks([1 Npcs]); 
                else
                    yticks(1);
                end
                ylabel('Cell #');
                xticks([1 15 30]); xticklabels([1 50 100]);
                xlabel('Position (cm)');
                title('Norm smooth pf map'); colorbar;
        end
        
        if Npcs_asd > 0
            subplot(10,8,6:8); imagesc(occMap(e,:));
                xticks([]); yticks([]); ylabel('Occ');
                title('ASD'); colorbar;
            subplot(10,8,[14:16,22:24,30:32]);
                imagesc(asd.normspkMap(asd.sortIdx(:,e),:,e));
                xticks([]);
                if Npcs_asd >  1
                    yticks([1 Npcs_asd]); 
                else
                    yticks(1);
                end
                ylabel('Cell #');
                title('Normalised spk map'); colorbar;
            subplot(10,8,[37,45,53]); 
                imagesc(asd.infoMap(asd.sortIdx,1,e));
                xticks([]);
                if Npcs_asd >  1
                    yticks([1 Npcs_asd]); 
                else
                    yticks(1);
                end 
                title('Max MI'); colorbar;
            subplot(10,8,[38:40,46:48,54:56]);    
                imagesc(asd.sort_normpfMap(:,:,e)); 
                if Npcs_asd >  1
                    yticks([1 Npcs_asd]); 
                else
                    yticks(1);
                end
                xticks([1 15 30]); xticklabels([1 50 100]);
                xlabel('Position (cm)');
                title('Normalised pf map'); colorbar;
        end

        if fsave
            if Nepochs == 1
                fname_fig = [sdir 'PFmaps/' mouseid '_' expname '_ref_' files(1,:) '_PFmaps_norm'];
            else
                fname_fig = [sdir 'PFmaps/' mouseid '_' expname '_ref_' files(1,:) '_PFmaps_norm_' num2str(e) 'of' num2str(Nepochs) 'ep'];
            end
            savefig( fh, fname_fig );
            saveas( fh, fname_fig, 'png' );
            if fclose, close( fh ); end
        end
    end

    % per trial spike maps for all cells
    [nRow, nCol] = getnRownCol(Ncells);
    nPlot = nRow*nCol;
    
    for ii=0:(Ncells/nPlot)-1
        fh = figure;
        ha = tight_subplot(nRow,nCol,[.02 .015],[.01 .05],[.01 .01]);
        for jj=0:nPlot-1
            if (ii*nPlot+jj+1) <= Ncells
                axes(ha(+jj+1));
                imagesc(normspkMap_pertrial{ii*nPlot+jj+1}); 
                yticks(ytick_files{ii*nPlot+jj+1}); yticklabels(ytick_files{ii*nPlot+jj+1}); % ylabel('Trial #');
                xticks([]); % xticklabels([1 50 100]); xlabel('Position (cm)');
                title(['Cell ' num2str(ii*nPlot+jj+1)],'fontsize',15);
            end
        end
        if fsave
            if Ncells/nPlot <= 1
                fname_fig = [sdir 'PFmaps/' mouseid '_' expname '_ref_' files(1,:) '_spk_pertrial_hist_norm'];
            else
                fname_fig = [sdir 'PFmaps/' mouseid '_' expname '_ref_' files(1,:) '_spk_pertrial_hist_norm_' num2str(ii+1)];
            end
            savefig( fh, fname_fig );
            saveas( fh, fname_fig, 'png' );
            if fclose, close( fh ); end
        end
    end 

    % per trial spike maps for place cells
    [nRow, nCol] = getnRownCol(Npcs);
    nPlot = nRow*nCol;
    
    % histogram
    if Npcs > 0 
        for ii=0:(Npcs/nPlot)-1
            fh = figure;
            ha = tight_subplot(nRow,nCol,[.01 .01],[.01 .05],[.01 .01]);
            for jj=0:nPlot-1
                if (ii*nPlot+jj+1) <= Npcs
                    axes(ha(+jj+1));
                    imagesc(hist.spkMap_pertrial{ii*nPlot+jj+1}); 
                    yticks(hist.ytick_files{ii*nPlot+jj+1}); yticklabels(hist.ytick_files{ii*nPlot+jj+1}); % ylabel('Trial #');
                    xticks([]); % xticklabels([1 50 100]); xlabel('Position (cm)');
                    title(['PC_h_i_s_t ' num2str(ii*nPlot+jj+1)],'fontsize',15);
                end
            end
            if fsave
                if Npcs/nPlot <= 1
                    fname_fig = [sdir 'PFmaps/' mouseid '_' expname '_ref_' files(1,:) '_PCs_hist_spk_pertrial'];
                else
                    fname_fig = [sdir 'PFmaps/' mouseid '_' expname '_ref_' files(1,:) '_PCs_hist_spk_pertrial' num2str(ii+1)];
                end
                savefig( fh, fname_fig );
                saveas( fh, fname_fig, 'png' );
                if fclose, close( fh ); end
            end
        end 

        for ii=0:(Npcs/nPlot)-1
            fh = figure;
            ha = tight_subplot(nRow,nCol,[.01 .01],[.01 .05],[.01 .01]);
            for jj=0:nPlot-1
                if (ii*nPlot+jj+1) <= Npcs
                    axes(ha(+jj+1));
                    imagesc(hist.normspkMap_pertrial{ii*nPlot+jj+1}); 
                    yticks(hist.ytick_files{ii*nPlot+jj+1}); yticklabels(hist.ytick_files{ii*nPlot+jj+1}); %ylabel('Trial #');
                    xticks([]); % xticklabels([1 50 100]); xlabel('Position (cm)');
                    title(['PC_h_i_s_t ' num2str(ii*nPlot+jj+1)],'fontsize',15);
                end
            end
            if fsave
                if Npcs/nPlot <= 1
                    fname_fig = [sdir 'PFmaps/' mouseid '_' expname '_ref_' files(1,:) '_PCs_hist_spk_pertrial_norm'];
                else
                    fname_fig = [sdir 'PFmaps/' mouseid '_' expname '_ref_' files(1,:) '_PCs_hist_spk_pertrial_norm' num2str(ii+1)];
                end
                savefig( fh, fname_fig );
                saveas( fh, fname_fig, 'png' );
                if fclose, close( fh ); end
            end
        end 
    end

    % asd
    [nRow, nCol] = getnRownCol(Npcs_asd);
    nPlot = nRow*nCol;
    
    if Npcs_asd > 0 
        for ii=0:(Npcs_asd/nPlot)-1
            fh = figure;
            ha = tight_subplot(nRow,nCol,[.01 .01],[.01 .05],[.01 .01]);
            for jj=0:nPlot-1
                if (ii*nPlot+jj+1) <= Npcs_asd
                    axes(ha(+jj+1));
                    imagesc(asd.spkMap_pertrial{ii*nPlot+jj+1}); 
                    yticks(asd.ytick_files{ii*nPlot+jj+1}); yticklabels(asd.ytick_files{ii*nPlot+jj+1}); % ylabel('Trial #');
                    xticks([]); % xticklabels([1 50 100]); xlabel('Position (cm)');
                    title(['PC_a_s_d ' num2str(ii*nPlot+jj+1)],'fontsize',15);
                end
            end
            if fsave
                if Npcs_asd/nPlot <= 1
                    fname_fig = [sdir 'PFmaps/' mouseid '_' expname '_ref_' files(1,:) '_PCs_asd_spk_pertrial'];
                else
                    fname_fig = [sdir 'PFmaps/' mouseid '_' expname '_ref_' files(1,:) '_PCs_asd_spk_pertrial' num2str(ii+1)];
                end
                savefig( fh, fname_fig );
                saveas( fh, fname_fig, 'png' );
                if fclose, close( fh ); end
            end
        end 

        for ii=0:(Npcs_asd/nPlot)-1
            fh = figure;
            ha = tight_subplot(nRow,nCol,[.01 .01],[.01 .05],[.01 .01]);
            for jj=0:nPlot-1
                if (ii*nPlot+jj+1) <= Npcs_asd
                    axes(ha(+jj+1));
                    imagesc(asd.normspkMap_pertrial{ii*nPlot+jj+1}); 
                    yticks(asd.ytick_files{ii*nPlot+jj+1}); yticklabels(asd.ytick_files{ii*nPlot+jj+1}); % ylabel('Trial #');
                    xticks([]); % xticklabels([1 50 100]); xlabel('Position (cm)');
                    title(['PC_a_s_d ' num2str(ii*nPlot+jj+1)],'fontsize',15);
                end
            end
            if fsave
                if Npcs_asd/nPlot <= 1
                    fname_fig = [sdir 'PFmaps/' mouseid '_' expname '_ref_' files(1,:) '_PCs_asd_spk_pertrial_norm'];
                else
                    fname_fig = [sdir 'PFmaps/' mouseid '_' expname '_ref_' files(1,:) '_PCs_asd_spk_pertrial_norm' num2str(ii+1)];
                end
                savefig( fh, fname_fig );
                saveas( fh, fname_fig, 'png' );
                if fclose, close( fh ); end
            end
        end 
    end

    % remapping within a session
    if Nepochs > 1
        if Npcs > 0
            fh = figure;
            for ei = 1:Nepochs % rows: sorting
                for ej = 1:Nepochs % cols: epochs 
                    subplot(Nepochs, Nepochs, (ei-1)*Nepochs + ej); imagesc(hist.normpfMap(hist.sortIdx(:,ei),:,ej)); 
                    title(['Epoch ' num2str(ej)]); ylabel(['Epoch' num2str(ei) ' sorting']);
                end
            end
            if fsave
                fname_fig = [sdir 'PFmaps/' mouseid '_' expname '_ref_' files(1,:) '_remapping_hist'];
                savefig( fh, fname_fig );
                saveas( fh, fname_fig, 'png' );
                if fclose, close( fh ); end
            end
        end

        if Npcs_asd > 0
            fh = figure;
            for ei = 1:Nepochs % rows: sorting
                for ej = 1:Nepochs % cols: epochs 
                    subplot(Nepochs, Nepochs, (ei-1)*Nepochs + ej); imagesc(asd.normpfMap(asd.sortIdx(:,ei),:,ej)); 
                    title(['Epoch ' num2str(ej)]); ylabel(['Epoch' num2str(ei) ' sorting']);
                end
            end
            if fsave
                fname_fig = [sdir 'PFmaps/' mouseid '_' expname '_ref_' files(1,:) '_remapping_asd'];
                savefig( fh, fname_fig );
                saveas( fh, fname_fig, 'png' );
                if fclose, close( fh ); end
            end
        end
    end
end