function plotSpikeRasterTrials( spkRaster, ytick_files, title_str, fsave, fname, fclose )
    Ncells = numel(spkRaster);
    Nbins = size(spkRaster{1},2);
    [nRow, nCol] = getnRownCol(Ncells);
    nPlot = nRow*nCol;
    Nfig = ceil(Ncells/nPlot)-1;
    if Nfig<0, Nfig = 0; end
    cmap0 = [0.9 0.9 0.9];
    cmap1 = [0 0 1];
    cmap = zeros(50,3);
    for j=1:3
        cmap(:,j) = linspace(cmap0(j),cmap1(j),50);
    end
    colormap(cmap);

    for ii=0:Nfig
        fh = figure;
        ha = tight_subplot(nRow,nCol,[.05 0.008],[.08 .05],[.06 .02]);
        for jj=0:nPlot-1
            if (ii*nPlot+jj+1) <= Ncells
                axes(ha(jj+1));
                imagesc(spkRaster{ii*nPlot+jj+1}); colormap(cmap);
                % only put ylabels for 1st column plots
                if ii*nPlot+jj+1 < 8
                    if ii*nPlot+jj+1 == 1
                        try
                            yticks(ytick_files); yticklabels(ytick_files); 
                            ylabel('Lap #');
                        catch
                        end
                    else
                        yticks([]);
                    end
                else
                    if mod(ii*nPlot+jj+1,nCol) == 1
                        try
                            yticks(ytick_files); yticklabels(ytick_files); 
                            ylabel('Lap #');
                        catch
                        end
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