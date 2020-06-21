% function for plotting neural events
% spikes: matrix of # of cells x time points
% fname_fig does not include file type (e.g. .fig, .jpg)

function plotSpikes(spikes, fname_fig, fclose)
    if nargin<3, fclose = true; end
    if nargin<2, fname_fig = []; end

    fig = figure;
    spklist = convert_mat2eventlist(spikes);
    plot_amplraster(spklist,'b',10);
    title('Neural events','Fontweight','normal','Fontsize',12);
    ylabel('Cell #');
    if ~isempty(fname_fig)
        savefig(fig, fname_fig);
        saveas(fig, fname_fig,'png');
    end
    if fclose
        close(fig);
    end
end