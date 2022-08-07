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

    % raw timeseries (of up to 50 cells)
    multiplot_ts(tsG, [figname_pref '_raw_timeseries'], 'Raw timeseries', fclose, 50)

    % dF/F (of up to 50 cells)
    multiplot_ts(df_f, [figname_pref '_df_f'], 'dF/F', fclose, 50)
end