function GUI_viewFISSAout(raw, result)

hdl_gui = figure('Name','FISSA result','NumberTitle','off','Resize','off',...
    'Position',[1087 1022 907 666]);

Numcells = numel(fieldnames(result));

str = sprintf('Cell number (1 to %g)', Numcells);
text_cellNum = uicontrol('Parent',hdl_gui,'style','text','Units','pixels','string',str,'Fontsize',16,...
    'FontWeight','bold'); 
    set(text_cellNum,'Position',[30 557 171 27]);
edit_cellNum = uicontrol('Parent',hdl_gui,'style','edit','Units','pixels','string','1','Fontsize',14,...
    'Callback',@edit_cellNum_callback); 
    set(edit_cellNum,'Position',[72 447 60 31]);
slider_cellNum = uicontrol('Parent',hdl_gui,'style','slider','Units','pixels',...
    'Min',1,'Max',Numcells,'Value',1,'SliderStep',[1/Numcells 1/Numcells],'Callback',@slider_cellNum_callback); 
    set(slider_cellNum,'Position',[59 499 92 28]);

ax1 = axes('Units','pixels','Position',[233 366 617 272],'Xlim',[0 1],'Ylim',[0  1],'Box','off',...
    'Visible','off','clipping','off'); 
ax2 = axes('Units','pixels','Position',[233 44 617 272],'Xlim',[0 1],'Ylim',[0  1],'Box','off',...
    'Visible','off','clipping','off'); 
plotResults(1);
    
% Callback functions
    function edit_cellNum_callback(varargin)
        id = str2double(get(edit_cellNum,'string'));

        if isnan(id)
            set(edit_cellNum,'string','1');
            set(slider_cellNum,'Value',1);
            id = 1;
        else
            if ~isinteger(id)
                id = round(id);
                set(edit_cellNum,'String',num2str(id));
            end
            if id < 1 
                set(edit_cellNum,'string','1');
                set(slider_cellNum,'Value',1);
                id = 1;
            elseif id > Numcells
                set(edit_cellNum,'string',num2str(Numcells));
                set(slider_cellNum,'Value',Numcells);
                id = Numcells;
            else
                set(edit_cellNum,'String',num2str(id));
            end
        end
        
        plotResults(id);
    end

    function slider_cellNum_callback(varargin)
        id = round(get(slider_cellNum,'Value'));
        set(edit_cellNum,'String',num2str(id));
        plotResults(id);
    end

    function plotResults(id)
        axes(ax1); 
            plot(raw.(['cell' num2str(id)]).trial0(1,:),'c'); 
            % axis([0 7500 0 14000]); 
            hold on
            plot(result.(['cell' num2str(id)]).trial0(1,:),'k'); 
            legend('raw','FISSA'); hold off
        axes(ax2); 
            plot(result.(['cell' num2str(id)]).trial0(1,:),'k'); 
            % axis([0 7500 0 14000]); 
            hold on
            plot(result.(['cell' num2str(id)]).trial0(2,:),'r'); hold on
            plot(result.(['cell' num2str(id)]).trial0(3,:),'g'); hold on
            plot(result.(['cell' num2str(id)]).trial0(4,:),'b'); hold on
            plot(result.(['cell' num2str(id)]).trial0(5,:),'m'); hold on
            legend('FISSA','c1','c2','c3','c4');
            hold off
    end
    

end