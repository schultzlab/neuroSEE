% Written by Ann Go
% 
% This GUI display mean 2P image and the ROI mask, raw Ca transient, ratiometric Ca
% transient (delta R / R) and the extracted spikes for each cell.
% User can choose cell and mode of spike extraction (NND vs OASIS).

function GUI_viewROIsSpikes( mean_imratio, masks, cell_tsG, R, spikes )
    
    % Create figure
    hfig_h = 500;
    hfig_w = 1200;
    hdl_gui = figure('MenuBar','none','Name','ROIs and spikes','NumberTitle','off',...%'Resize','off',...
    'Position',[1680-hfig_w,1050-hfig_h,hfig_w,hfig_h]);
    
    % Display cell number and allow user to change cell number either by
    % using a slider or an edit box
    Numcells = size(masks,3);
    str = sprintf('Cell number: 1 to %g', Numcells);
    text_cellNum = uicontrol('Parent',hdl_gui,'style','text','Units','normalized','string',str,'Fontsize',12); 
        set(text_cellNum,'Position',[0.45 0.16 0.1 0.05]);
    edit_cellNum = uicontrol('Parent',hdl_gui,'style','edit','Units','normalized','string','1','Callback',@edit_cellNum_callback); 
        set(edit_cellNum,'Position',[0.45 0.11 0.04 0.05]);
    slider_cellNum = uicontrol('Parent',hdl_gui,'style','slider','Units','normalized',...
        'Min',1,'Max',Numcells,'Value',1,'SliderStep',[1/Numcells 1/Numcells],'Callback',@slider_cellNum_callback); 
        set(slider_cellNum,'Position',[0.5 0.10 0.1 0.05]);

    % Allow option to show all ROIs
    check_allROI = uicontrol('Parent',hdl_gui,'style','check','Units','normalized',...
        'string','Show all ROIs','Fontsize',11,'Position',[0.68 0.11 0.08 0.05],...
        'Value',0,'Callback',@check_allROI_callback);
    
    % Show 2P image and overlay mask of specific cell
    ax_masks = axes('Position',[0.01 0.03 0.4 0.93],'Xlim',[0 1],'Ylim',[0  1],'Box','off',...
                'Visible','off','Units','normalized', 'clipping' , 'off'); 
    imagesc(mean_imratio); colormap(gray); axis off;
    hold on
    outline = bwboundaries(masks(:,:,1));
    trace = outline{1};
    plot(trace(:,2),trace(:,1),'g','Linewidth',3);
    hold off
    
    area1 = bwarea(masks(:,:,1));
    str = sprintf('ROI area: %g pixels', area1);
    text_ROIarea = uicontrol('Parent',hdl_gui,'style','text','Units','normalized','string',str,'Fontsize',12); 
        set(text_ROIarea,'Position',[0.45 0.04 0.15 0.05]);

    % Plot Ca time series for specific cell
    ax_ts = axes('Position',[0.45 0.78 0.52 0.15],'Xlim',[0 1],'Ylim',[0  1],'Box','off',...
                'Visible','off','Units','normalized', 'clipping' , 'off'); 
    plot(cell_tsG(1,:)); 
    axes('Position',[0.45 0.96 0.52 0.05],'Xlim',[0 1],'Ylim',[0  1],'Box','off',...
                'Visible','off','Units','normalized', 'clipping' , 'off'); 
        text(0.4, 0, 'Raw Ca^{2+} time series','Fontsize',12,'Fontweight','bold');

    % Plot ratiometric Ca time series (delta R/R) for specific cell
    ax_R = axes('Position',[0.45 0.54 0.52 0.15],'Xlim',[0 1],'Ylim',[0  1],'Box','off',...
                'Visible','off','Units','normalized', 'clipping' , 'off'); 
    plot(R(1,:)); 
    axes('Position',[0.45 0.71 0.52 0.05],'Xlim',[0 1],'Ylim',[0  1],'Box','off',...
                'Visible','off','Units','normalized', 'clipping' , 'off'); 
        text(0.46, 0, 'deltaR / R','Fontsize',12,'Fontweight','bold');

    % Plot spikes for specific cell. Allow user to choose between 
    ax_spikes = axes('Position',[0.45 0.30 0.52 0.15],'Xlim',[0 1],'Ylim',[0  1],'Box','off',...
                'Visible','off','Units','normalized', 'clipping' , 'off'); 
    plot(spikes(1,:)); 
    axes('Position',[0.45 0.47 0.52 0.05],'Xlim',[0 1],'Ylim',[0  1],'Box','off',...
                'Visible','off','Units','normalized', 'clipping' , 'off'); 
        text(0.47, 0, 'Spikes','Fontsize',12,'Fontweight','bold');

    
    function plotROIsSpikes(id)
        % Display neuron image overlaid with specified ROI mask
        axes(ax_masks);
        imagesc(mean_imratio); colormap(gray); axis off;
        hold on
        if get(check_allROI,'Value') == 1
            for j = 1:Numcells
                outline = bwboundaries(masks(:,:,j));
                trace = outline{1};
                plot(trace(:,2),trace(:,1),'b','Linewidth',2);
            end
        end
        outline = bwboundaries(masks(:,:,id));
        trace = outline{1};
        plot(trace(:,2),trace(:,1),'g','Linewidth',3);
        hold off
    
        % Plot Ca time series for specific cell
        plot(ax_ts, cell_tsG(id,:)); 

        % Plot ratiometric Ca time series (delta R/R) for specific cell
        plot(ax_R, R(id,:)); 
    
        % Plot spikes for specific cell. Allow user to choose between 
        plot(ax_spikes, spikes(id,:)); 

    end

    function edit_cellNum_callback(varargin)
        id = str2double(get(edit_cellNum,'string'));
        
        if isnan(id)
            set(edit_cellNum,'string','1');
            set(slider_cellNum,'Value',1);
        else
            if id < 1 
                set(edit_cellNum,'string','1');
                set(slider_cellNum,'Value',1);
                str = sprintf('ROI area: %g pixels', area1);
                    set(text_ROIarea,'string',str);
            elseif id > Numcells
                set(edit_cellNum,'string',num2str(Numcells));
                set(slider_cellNum,'Value',Numcells);
                area = bwarea(masks(:,:,Numcells));
                    str = sprintf('ROI area: %g pixels', area);
                    set(text_ROIarea,'string',str);
            else
                plotROIsSpikes(id);
                set(slider_cellNum,'Value',id);
                area = bwarea(masks(:,:,id));
                    str = sprintf('ROI area: %g pixels', area);
                    set(text_ROIarea,'string',str);
            end
        end
    end

    function slider_cellNum_callback(varargin)
        id = round(get(slider_cellNum,'Value'));
        plotROIsSpikes(id);
        set(edit_cellNum,'string',num2str(id));
        area = bwarea(masks(:,:,id));
            str = sprintf('ROI area: %g pixels', area);
            set(text_ROIarea,'string',str);
    end

    function check_allROI_callback(varargin)
        id = str2double(get(edit_cellNum,'string'));
        plotROIsSpikes(id);  
    end
        
end