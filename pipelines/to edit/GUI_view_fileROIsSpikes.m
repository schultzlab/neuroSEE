% Written by Ann Go
% 
% This GUI allows user to view correlation image, masks, raw Ca timeseries,
% dF/F and spikes
%
% INPUTS:
%   spikes      : spike estimates (optional)
%   tsG         : raw timeseries 
%   dtsG        : fissa-corrected timeseries (optional) 
%   df_f        : dF/F (optional)
%   ddf_f       : fissa-corrected dF/F (optional)
%   R           : ratiometric timeseries R = Fgreen/Fred (optional)
%   corr_image  : correlation image from green channel (optional)
%   masks       : ROI masks (optional)
%
% Note: one of these {df_f, ddf_f, R} must not be empty
% Note: for optional inputs, use empty matrix
%   e.g. [spikes, params] = GUI_manually_refine_spikes( spikes, tsG, [], df_f, [], [], [], [])


function GUI_view_ROIsSpikes( spikes, tsG, dtsG, df_f, ddf_f, R, corr_image, masks )

N = size(tsG,1); 

%% GUI structures

% GUI figure
hfig_h = 581;
hfig_w = 1346;
fig_gui = figure('Name','Manually refine spike estimates','NumberTitle','off',...%'Resize','off',...
'Position',[1680-hfig_w,1050-hfig_h,hfig_w,hfig_h]);

% correlation image 
ax_masks = axes('Units','pixels','Position',[27 122 448 435],'Xlim',[0 1],'Ylim',[0  1],'Box','off',...
            'Visible','on','clipping' , 'off'); 

% cell number 
str = sprintf('Cell number: 1 to %g', N);
text_cellNum = uicontrol('Parent',fig_gui,'style','text','Units','pixels','string',str,'Fontsize',14); 
    set(text_cellNum,'Position',[50 69 145 26]);
edit_cellNum = uicontrol('Parent',fig_gui,'style','edit','Units','pixels','string','1','Callback',@edit_cellNum_callback); 
    set(edit_cellNum,'Position',[56 31 40 28]);
slider_cellNum = uicontrol('Parent',fig_gui,'style','slider','Units','pixels',...
    'Min',1,'Max',N,'Value',1,'SliderStep',[1/N 1/N],'Callback',@slider_cellNum_callback); 
    set(slider_cellNum,'Position',[114 28 136 24]);

% show all ROIs
check_allROI = uicontrol('Parent',fig_gui,'style','check','Units','pixels',...
    'string','Show all ROIs','Fontsize',11,'Position',[367 34 102 23],...
    'Value',0,'Callback',@check_allROI_callback);
      
% Ca time series plot
ax_F = axes('Parent',fig_gui,'Units','pixels','Position',[534 445 785 95],'Xlim',[0 1],'Ylim',[0  1],'Box','off',...
            'Visible','on','clipping','off'); 

% df/f 
ax_dF = axes('Parent',fig_gui,'Units','pixels','Position',[534 313 785 95],'Xlim',[0 1],'Ylim',[0  1],'Box','off',...
            'Visible','on','clipping','off');

% spikes plot
ax_spikes = axes('Parent',fig_gui,'Units','pixels','Position',[534 180 785 95],'Xlim',[0 1],'Ylim',[0  1],'Box','off',...
            'Visible','on','clipping','off'); 

    
%% initial GUI data
if ~isempty(dtsG)
    C1 = dtsG; str_C1 = 'FISSA-corrected timeseries';
else
    C1 = tsG; str_C1 = 'Raw timeseries';
end
if ~isempty(ddf_f)
    C2 = ddf_f; str_C2 = 'FISSA-corrected dF/F';
elseif ~isempty(df_f)
    C2 = df_f; str_C2 = 'dF/F';
elseif ~isempty(R)
    C2 = R; str_C2 = 'dR/R';
else
    C2 = zeros(size(tsG)); str_C2 = 'dF/F No: data';
end
if ~isempty(spikes)
    str_spikes = 'Spikes';
else
    spikes = zeros(size(tsG)); str_spikes = 'Spikes: No data';
end

if nargin < 8, masks = []; end 
if nargin < 7, corr_image = zeros(512,512); end
displayROI(1); showTS(1);


%% GUI Subfunctions
function displayROI(id)
    % Display correlation image overlaid with specified ROI mask
    axes(ax_masks);
    imagesc(corr_image); axis off; colormap(gray);
    if ~isempty(masks)
        hold on
        if get(check_allROI,'Value') == 1
            for j = 1:N
                outline = bwboundaries(masks(:,:,j));
                if size(outline,1) > 0
                    trace = outline{1};
                    plot(trace(:,2),trace(:,1),'w','Linewidth',2);
                end
            end
        end
        outline = bwboundaries(masks(:,:,id));
        if size(outline,1) > 0
            trace = outline{1};
            plot(trace(:,2),trace(:,1),'g','Linewidth',3);
        end
        hold off
    end
end

function edit_cellNum_callback(varargin)
    id = str2double(get(edit_cellNum,'string'));

    if isnan(id)
        set(edit_cellNum,'string','1');
        set(slider_cellNum,'Value',1);
        displayROI(1); showTS(1);
    else
        if id < 1 
            set(edit_cellNum,'string','1');
            set(slider_cellNum,'Value',1);
            displayROI(1); showTS(1);
        elseif id > N
            set(edit_cellNum,'string',num2str(N));
            set(slider_cellNum,'Value',N);
            displayROI(N); showTS(N);
        else
            set(slider_cellNum,'Value',id);
            displayROI(id); showTS(id);
        end
    end
end

function slider_cellNum_callback(varargin)
    id = round(get(slider_cellNum,'Value'));
    set(edit_cellNum,'string',num2str(id));
    displayROI(id); showTS(id);
end

function check_allROI_callback(varargin)
    id = str2double(get(edit_cellNum,'string'));
    displayROI(id);    
end

function showTS(id)
    % Ca time series
    axes(ax_F);
    plot(C1(id,:)); title(str_C1);

    % df/f or dR/R
    axes(ax_dF);
    plot(C2(id,:)); title(str_C2);

    % spikes
    axes(ax_spikes);
    plot(spikes(id,:)); title(str_spikes);
end

end