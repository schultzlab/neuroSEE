% fname_fig does not include file type (e.g. .fig, .jpg)

function multiplot_ts(ts, fname_fig, title_str, fclose)
    if nargin<4, fclose = true; end
    if nargin<3, title_str = []; end
    if nargin<2, fname_fig = []; end
    
    fig = figure;
    iosr.figures.multiwaveplot(1:size(ts,2),1:size(ts,1),ts,'gain',5); yticks([]); xticks([]);
    if ~isempty(title_str)
        title(title_str,'Fontweight','normal','Fontsize',12);
    end
    if ~isempty(fname_fig)
        savefig(fig, fname_fig);
        saveas(fig, fname_fig,'png');
    end
    if fclose
        close(fig);
    end
end

