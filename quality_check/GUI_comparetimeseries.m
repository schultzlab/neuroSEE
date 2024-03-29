function GUI_comparetimeseries(ts1, ts2, ts3, ts4, ts5, label1, label2, label3, label4, label5, colour1, colour2, colour3, colour4, colour5)
if nargin<15, colour{5} = 'b'; else, colour{5} = colour5; end % b = blue
if nargin<14, colour{4} = 'b'; else, colour{4} = colour4; end
if nargin<13, colour{3} = 'b'; else, colour{3} = colour3; end
if nargin<12, colour{2} = 'b'; else, colour{2} = colour2; end
if nargin<11, colour{1} = 'b'; else, colour{1} = colour1; end
if nargin<10, label{5} = []; else, label{5} = label5; end
if nargin<9, label{4} = []; else, label{4} = label4; end
if nargin<8, label{3} = []; else, label{3} = label3; end
if nargin<7, label{2} = []; else, label{2} = label2; end
if nargin<6, label{1} = []; else, label{1} = label1; end
if nargin<5, ts{5} = []; else, ts{5} = ts5; end
if nargin<4, ts{4} = []; else, ts{4} = ts4; end
if nargin<3, ts{3} = []; else, ts{3} = ts3; end
if nargin<2, ts{2} = []; else, ts{2} = ts2; end
ts{1} = ts1;

hdl_gui = figure('Name','','NumberTitle','off','Resize','off',...
    'Position',[1087 1022 1600 550]);

Numcells = size(ts{1},1);

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

ax{1} = axes('Units','pixels','Position',[250 450 1300 72],'Xlim',[0 1],'Ylim',[0  1],'Box','off',...
    'Visible','on','clipping','off'); 
ax{2} = axes('Units','pixels','Position',[250 350 1300 72],'Xlim',[0 1],'Ylim',[0  1],'Box','off',...
    'Visible','on','clipping','off'); 
ax{3} = axes('Units','pixels','Position',[250 250 1300 72],'Xlim',[0 1],'Ylim',[0  1],'Box','off',...
    'Visible','on','clipping','off'); 
ax{4} = axes('Units','pixels','Position',[250 150 1300 72],'Xlim',[0 1],'Ylim',[0  1],'Box','off',...
    'Visible','on','clipping','off'); 
ax{5} = axes('Units','pixels','Position',[250 50 1300 72],'Xlim',[0 1],'Ylim',[0  1],'Box','off',...
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
        for i = 1:numel(ax)
            if ~isempty(ts{i})
                axes(ax{i}); 
                x = (1:length(ts{i}(id,:)))/30.91;
                plot(x,ts{i}(id,:),colour{i});
                if ~isempty(label{i}), ylabel(label{i}); end
            end
        end
    end
end