function plotspikes(spikes, fname_fig)
    fig = figure;
    iosr.figures.multiwaveplot(1:size(spikes,2),1:size(spikes,1),spikes,'gain',5); yticks([]); xticks([]);
    title('dF/F','Fontweight','normal','Fontsize',12);
    savefig(fig,[fname_fig '_spikes']);
    saveas(fig,[fname_fig '_spikes'],'png');
    close(fig);
end

