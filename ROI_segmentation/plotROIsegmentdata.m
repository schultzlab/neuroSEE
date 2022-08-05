function plotROIsegmentdata(corr_image, masks, elim_masks, tsG, df_f, figname_pref, fclose)
    if nargin<7, fclose = true; end
    if nargin<6, figname_pref = []; end
    
    % ROIs overlayed on correlation image
    plotopts.plot_ids = 1; % set to 1 to view the ID number of the ROIs on the plot
    fig = plotContoursOnSummaryImage(corr_image, masks, plotopts);
    if ~isempty(figname_pref)
        savefig(fig, [figname_pref '_ROIs']);
        saveas(fig, [figname_pref '_ROIs'], 'png');
    end
    if fclose, close(fig); end

    % eliminated ROIs overlayed on correlation image
    if ~isempty(elim_masks)
        plotopts.plot_ids = 1; % set to 1 to view the ID number of the ROIs on the plot
        fig = plotContoursOnSummaryImage(corr_image, elim_masks, plotopts);
        if ~isempty(figname_pref)
            savefig(fig, [figname_pref '_elimROIs']);
            saveas(fig, [figname_pref '_elimROIs'], 'png');
        end
        if fclose, close(fig); end
    end

    % raw timeseries
    fig = figure;
    y = max(size(tsG,1),100);
    iosr.figures.multiwaveplot(1:size(tsG,2),1:y,tsG(1:y,:),'gain',5); yticks([]); xticks([]); 
    title('Raw timeseries','Fontweight','normal','Fontsize',12); 
    if ~isempty(figname_pref)
        savefig(fig, [figname_pref '_raw_timeseries']);
        saveas(fig, [figname_pref '_raw_timeseries'], 'png');
    end
    if fclose, close(fig); end

    % dF/F
    fig = figure;
    iosr.figures.multiwaveplot(1:size(df_f,2),1:y,df_f(1:y,:),'gain',5); yticks([]); xticks([]); 
    title('dF/F','Fontweight','normal','Fontsize',12); 
    if ~isempty(figname_pref)
        savefig(fig, [figname_pref '_df_f']);
        saveas(fig, [figname_pref '_df_f'], 'png');
    end
    if fclose, close(fig); end
end