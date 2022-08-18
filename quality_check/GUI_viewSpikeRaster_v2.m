% Written by Ann Go
% GUI for viewing spike rasters

function GUI_viewSpikeRaster_v2(PFdata)

%% GUI structures
% GUI figure
hfig = figure('Name','Spikes','NumberTitle','off','Resize','off',...
    'Position',[1000 150 200 240]);

% image 
ax_image = axes('Units','pixels','Position',[30 20 150 150],'Box','off','Visible','on','clipping','off'); 

% data
raster = PFdata.normspkRaster;
ytick_files = PFdata.ytick_files;
Nrois = numel(raster);

% Cell number
str = sprintf('ROI number (1 to %g):', Nrois);
text_roiNum = uicontrol('style','text','Units','pixels','string',str,'Fontsize',12,...
    'FontWeight','bold','Position',[20 210 171 23]); 
edit_roiNum = uicontrol('style','edit','Units','pixels','string','1','Fontsize',14,...
    'Position',[20 185 45 26],'Callback',@edit_roiNum_callback); 
slider_roiNum = uicontrol('style','slider','Units','pixels','Position',[85 185 101 22],...
    'Min',1,'Max',Nrois,'Value',1,'SliderStep',[1/Nrois 1/Nrois],'Callback',@slider_roiNum_callback); 

% Initialise GUI
plotRaster(1);

%% Subfunctions and Callback functions
function plotRaster(id)
    axes(ax_image);
    imagesc(raster{id}); colormap(flipud(gray));
    yticks(ytick_files);
    xticks([]);
end

function edit_roiNum_callback(varargin)
    id = str2double(get(edit_roiNum,'string'));
    if isnan(id)
        set(edit_roiNum,'string','1');
        set(slider_roiNum,'Value',1);
        plotRaster(1);
    else
        if ~isinteger(id)
            id = round(id);
            set(edit_roiNum,'String',num2str(id));
        end
        if id < 1 
            set(edit_roiNum,'string','1');
            set(slider_roiNum,'Value',1);
            plotRaster(1);
        elseif id > Nrois
            set(edit_roiNum,'string',num2str(Nrois));
            set(slider_roiNum,'Value',Nrois);
            plotRaster(Nrois);
        else
            set(edit_roiNum,'String',num2str(id));
            set(slider_roiNum,'Value',id);
            plotRaster(id);
        end
    end
end

function slider_roiNum_callback(varargin)
    id = round(get(slider_roiNum,'Value'));
    plotRaster(id);
    set(edit_roiNum,'string',num2str(id));
end


end