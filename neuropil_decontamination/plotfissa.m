function plotfissa(dtsG, ddf_f, fname_fig)
    % raw timeseries
    fig = figure;
    iosr.figures.multiwaveplot(1:size(dtsG,2),1:size(dtsG,1),dtsG,'gain',5); yticks([]); xticks([]); 
    title('Fissa-corrected raw timeseries','Fontweight','normal','Fontsize',12); 
    savefig(fig,[fname_fig '_fissa_result']);
    saveas(fig,[fname_fig '_fissa_result'],'png');
    close(fig);

    % dF/F
    fig = figure;
    iosr.figures.multiwaveplot(1:size(ddf_f,2),1:size(ddf_f,1),ddf_f,'gain',5); yticks([]); xticks([]); 
    title('Fissa-corrected dF/F','Fontweight','normal','Fontsize',12); 
    savefig(fig,[fname_fig '_fissa_df_f']);
    saveas(fig,[fname_fig '_fissa_df_f'],'png');
    close(fig);
end
