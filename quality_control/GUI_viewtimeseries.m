function GUI_viewtimeseries(tsG, df_f, spikes, dtsG, ddf_f)

hdl_gui = figure('Name','FISSA result','NumberTitle','off','Resize','off',...
    'Position',[1087 1022 1600 550]);

Numcells = size(tsG,1);

str = sprintf('Cell number (1 to %g)', Numcells);
text_cellNum = uicontrol('Parent',hdl_gui,'style','text','Units','pixels','string',str,'Fontsize',14,...
    'FontWeight','bold'); 
    set(text_cellNum,'Position',[10 487 180 27]);
    edit_cellNum = uicontrol('Parent',hdl_gui,'style','edit','Units','pixels','string','1','Fontsize',14,...
        'Callback',@edit_cellNum_callback); 
        set(edit_cellNum,'Position',[22 457 60 31]);
    slider_cellNum = uicontrol('Parent',hdl_gui,'style','slider','Units','pixels',...
        'Min',1,'Max',Numcells,'Value',1,'SliderStep',[1/Numcells 1/Numcells],'Callback',@slider_cellNum_callback); 
        set(slider_cellNum,'Position',[93 449 92 28]);
% text_mouseID = uicontrol('Parent',hdl_gui,'style','text','Units','pixels','string','m62','Fontsize',12,...
%     'FontWeight','bold','HorizontalAlignment','left'); 
%     set(text_mouseID,'Position',[25 400 50 27]);
% text_fov = uicontrol('Parent',hdl_gui,'style','text','Units','pixels','string','fov 1of2','Fontsize',12,...
%     'FontWeight','bold','HorizontalAlignment','left'); 
%     set(text_fov,'Position',[25 380 50 27]);
% text_genotype = uicontrol('Parent',hdl_gui,'style','text','Units','pixels','string','WT','Fontsize',12,...
%     'FontWeight','bold','HorizontalAlignment','left'); 
%     set(text_genotype,'Position',[25 360 50 27]);
% text_ageMOS = uicontrol('Parent',hdl_gui,'style','text','Units','pixels','string','5.6 mos','Fontsize',12,...
%     'FontWeight','bold','HorizontalAlignment','left'); 
%     set(text_ageMOS,'Position',[25 340 50 27]);

ax1 = axes('Units','pixels','Position',[250 450 1300 72],'Xlim',[0 1],'Ylim',[0  1],'Box','off',...
    'Visible','on','clipping','off'); 
ax2 = axes('Units','pixels','Position',[250 350 1300 72],'Xlim',[0 1],'Ylim',[0  1],'Box','off',...
    'Visible','on','clipping','off'); 
ax3 = axes('Units','pixels','Position',[250 250 1300 72],'Xlim',[0 1],'Ylim',[0  1],'Box','off',...
    'Visible','on','clipping','off'); 
ax4 = axes('Units','pixels','Position',[250 150 1300 72],'Xlim',[0 1],'Ylim',[0  1],'Box','off',...
    'Visible','on','clipping','off'); 
ax5 = axes('Units','pixels','Position',[250 50 1300 72],'Xlim',[0 1],'Ylim',[0  1],'Box','off',...
    'Visible','on','clipping','off'); 

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
                set(slider_cellNum,'Value',id);
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
        x = (1:length(tsG(id,:)))/30.91;
        axes(ax1); 
            plot(x,tsG(id,:),'b');
            ylabel('Fluor');
        axes(ax2); 
            plot(x,df_f(id,:),'b');
            ylabel('df/f');
        axes(ax3); 
            plot(x,spikes(id,:),'b');
            ylabel('Spikes');
            xlabel('Time (s)');
%         axes(ax4); 
%             plot(x,dtsG(id,:),'r');
%             ylabel('df/f detrended');
%         axes(ax5); 
%             plot(x,ddf_f(id,:),'r');
%             ylabel('Spikes');
%         set(text_mouseID,'string',mouseID{id});
%         set(text_fov,'string',['fov ' fov{id}]);
%         set(text_genotype,'string',genotype{id});
%         set(text_ageMOS,'string',num2str(ageMOS{id}));
    end
    

end