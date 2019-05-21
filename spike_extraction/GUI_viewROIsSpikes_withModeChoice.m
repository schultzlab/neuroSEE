% Written by Ann Go
% 
% This GUI displays the mean 2P image and the ROI mask, raw Ca transient, ratiometric Ca
% transient (delta R / R) and the extracted spikes for each cell.
% User can choose cell and mode of spike extraction (NND vs OASIS) and can
% tweak parameters for ROIs and spikes.

function GUI_viewROIsSpikes_withModeChoice(mean_R,masks,tsG,R,spikes_nnd,spikes_oasis)
    
    %% GUI variables
    Fmin = min(min(tsG));
    Fmax= 12000; %max(max(tsG));
    Rmin = min(min(R));
    Rmax = 2; % max(max(R));
    spikemax = 0.8;
    satThr = 0.95;
    satTime = 0.3;
    
    %% GUI structures

    % Create figure
    hfig_h = 500;
    hfig_w = 1200;
    hdl_gui = figure('MenuBar','none','Name','ROIs and spikes','NumberTitle','off',...%'Resize','off',...
    'Position',[1680-hfig_w,1050-hfig_h,hfig_w,hfig_h]);

    % Initialise GUI data
    setappdata(hdl_gui,'curr_masks',masks);
    setappdata(hdl_gui,'curr_tsG',tsG);
    setappdata(hdl_gui,'curr_R',R);
    setappdata(hdl_gui,'ALLspikes_nnd',spikes_nnd);
    setappdata(hdl_gui,'ALLspikes_oasis',spikes_oasis);
    setappdata(hdl_gui,'curr_spikes_nnd',spikes_nnd);
    setappdata(hdl_gui,'curr_spikes_oasis',spikes_oasis);
    setappdata(hdl_gui,'curr_spikes',spikes_nnd);
    Numcells = size(masks,3); 
        setappdata(hdl_gui,'curr_Numcells',Numcells);
        setappdata(hdl_gui,'curr_ind',1:Numcells);
    
    
    % Display cell number and allow user to change cell number either by
    % using a slider or an edit box
    str = sprintf('Cell number: 1 to %g', Numcells);
    text_cellNum = uicontrol('Parent',hdl_gui,'style','text','Units','normalized','string',str,'Fontsize',12); 
        set(text_cellNum,'Position',[0.42 0.245 0.1 0.05]);
    edit_cellNum = uicontrol('Parent',hdl_gui,'style','edit','Units','normalized','string','1','Callback',@edit_cellNum_callback); 
        set(edit_cellNum,'Position',[0.42 0.20 0.04 0.04]);
    slider_cellNum = uicontrol('Parent',hdl_gui,'style','slider','Units','normalized',...
        'Min',1,'Max',Numcells,'Value',1,'SliderStep',[1/Numcells 1/Numcells],'Callback',@slider_cellNum_callback); 
        set(slider_cellNum,'Position',[0.47 0.185 0.08 0.05]);
    
    % Allow option to show all ROIs
    check_allROI = uicontrol('Parent',hdl_gui,'style','check','Units','normalized',...
        'string','Show all ROIs','Fontsize',11,'Position',[0.42 0.13 0.08 0.05],...
        'Value',0,'Callback',@check_allROI_callback);
    
    % Show 2P image and overlay mask of specific cell
    ax_masks = axes('Position',[0.01 0.04 0.4 0.93],'Xlim',[0 1],'Ylim',[0  1],'Box','off',...
                'Visible','off','Units','normalized', 'clipping' , 'off'); 
    displayROI(1);
    
    % Calculate the area and SNR for each ROI 
    maskArea = zeros(1,Numcells);
    for i = 1:Numcells
        maskArea(i) = bwarea(masks(:,:,i));
    end
    maskSNR = zeros(1,Numcells);
    for i = 1:Numcells
        y = tsG(i,:);
        maskSNR(i) = GetSn(y,[0.25,0.5],'logmexp');
    end
    setappdata(hdl_gui,'curr_maskArea',maskArea);
    setappdata(hdl_gui,'curr_maskSNR',maskSNR);
    
    % Display ROI area and SNR
    area1 = maskArea(1);
    str1 = sprintf('Area: %g pixels', round(area1,1));
    text_ROIarea = uicontrol('Parent',hdl_gui,'style','text','Units','normalized',...
    'string',str1,'Fontsize',12,'HorizontalAlignment','left'); 
        set(text_ROIarea,'Position',[0.58 0.245 0.15 0.05]);
    
    snr1 = maskSNR(1);
    str2 = sprintf('Noise: %g', round(snr1,2));
    text_ROIsnr = uicontrol('Parent',hdl_gui,'style','text','Units','normalized',...
    'string',str2,'Fontsize',12,'HorizontalAlignment','left'); 
        set(text_ROIsnr,'Position',[0.58 0.20 0.15 0.05]);
     
     uicontrol('Parent',hdl_gui,'style','text','Units','normalized','Position',[0.58 0.155 0.15 0.05],...
        'string','Saturation:','Fontsize',12,'HorizontalAlignment','left'); 
        %set(text_ROIsnr,);

        
    % Allow user to set threshold for area and/or SNR
    check_areaThr = uicontrol('Parent',hdl_gui,'style','check','Units','normalized',...
        'string','Thr','Fontsize',11,'Position',[0.66 0.25 0.05 0.05],...
        'Value',0,'Callback',@AreaSNRSatThr_callback);
    edit_areaThr = uicontrol('Parent',hdl_gui,'style','edit','Units','normalized',...
        'string','70','Enable','off','Callback',@AreaSNRSatThr_callback); 
        set(edit_areaThr,'Position',[0.7 0.255 0.04 0.04]);
        
    check_snrThr = uicontrol('Parent',hdl_gui,'style','check','Units','normalized',...
        'string','Thr','Fontsize',11,'Position',[0.66 0.205 0.12 0.05],...
        'Value',0,'Callback',@AreaSNRSatThr_callback);
    edit_snrThr = uicontrol('Parent',hdl_gui,'style','edit','Units','normalized',...
        'string','490','Enable','off','Callback',@AreaSNRSatThr_callback); 
        set(edit_snrThr,'Position',[0.7 0.21 0.04 0.04]);

    check_satThr = uicontrol('Parent',hdl_gui,'style','check','Units','normalized',...
        'string','Thr','Fontsize',11,'Position',[0.66 0.16 0.12 0.05],...
        'Value',0,'Callback',@AreaSNRSatThr_callback);
    edit_satThr = uicontrol('Parent',hdl_gui,'style','edit','Units','normalized',...
        'string','3000','Enable','off','Callback',@AreaSNRSatThr_callback); 
        set(edit_satThr,'Position',[0.7 0.165 0.04 0.04]);


    % Plot Ca time series for specific cell
    ax_ts = axes('Position',[0.45 0.82 0.52 0.13],'Xlim',[0 1],'Ylim',[0  1],'Box','off',...
                'Visible','off','Units','normalized', 'clipping' , 'off'); 
    plot(tsG(1,:)); axis([0 7500 Fmin Fmax]);
    axes('Position',[0.45 0.98 0.52 0.05],'Xlim',[0 1],'Ylim',[0  1],'Box','off',...
                'Visible','off','Units','normalized', 'clipping' , 'off'); 
        text(0.4, 0, 'Raw Ca^{2+} time series','Fontsize',12,'Fontweight','bold');

    % Plot ratiometric Ca time series (delta R/R) for specific cell
    ax_R = axes('Position',[0.45 0.6 0.52 0.13],'Xlim',[0 1],'Ylim',[0  1],'Box','off',...
                'Visible','off','Units','normalized', 'clipping' , 'off'); 
    plot(R(1,:)); axis([0 7500 Rmin Rmax]);
    axes('Position',[0.45 0.75 0.52 0.05],'Xlim',[0 1],'Ylim',[0  1],'Box','off',...
                'Visible','off','Units','normalized', 'clipping' , 'off'); 
        text(0.46, 0, 'deltaR / R','Fontsize',12,'Fontweight','bold');

    % Plot spikes for specific cell. 
    ax_spikes = axes('Position',[0.45 0.38 0.52 0.13],'Xlim',[0 1],'Ylim',[0  1],'Box','off',...
                'Visible','off','Units','normalized', 'clipping' , 'off'); 
    plot(spikes_nnd(1,:)); axis([0 7500 0 spikemax]);
    axes('Position',[0.45 0.53 0.52 0.05],'Xlim',[0 1],'Ylim',[0  1],'Box','off',...
                'Visible','off','Units','normalized', 'clipping' , 'off'); 
        text(0.47, 0, 'Spikes','Fontsize',12,'Fontweight','bold');
        
    % Allow user to choose between NND and OASIS. NND is the default
    axes('Position',[0.82 0.28 0.1 0.05],'Xlim',[0 1],'Ylim',[0  1],'Box','off',...
                'Visible','off','Units','normalized', 'clipping' , 'off'); 
        text(0, 0, 'Choose spike extraction mode','Fontsize',12);
    
    % Allow user to choose NND and to tweak del F/ F threshold for Ca transient
    panel_nnd = uipanel('Parent',hdl_gui,'Units','normalized','Position',[0.83 0.16 0.14 0.09]);
    check_nnd = uicontrol('Parent',panel_nnd,'style','check','Units','normalized',...
        'string','NND','Fontsize',11,'Position',[0.02 0.6 0.8 0.34],...
        'Value',0,'Callback',@check_nnd_callback);
    text_Rthr = uicontrol('Parent',panel_nnd,'style','text','Units','normalized','string','R threshold',...
        'Enable','off','Fontsize',11); 
        set(text_Rthr,'Position',[0.07 0.1 0.5 0.42]);
    edit_Rthr = uicontrol('Parent',panel_nnd,'style','edit','Units','normalized',...
        'Enable','off','string','0.022','Callback',@edit_Rthr_callback); 
        set(edit_Rthr,'Position',[0.68 0.1 0.28 0.42]);

    % Allow user to choose OASIS and to tweak OASIS parameters
    panel_oasis = uipanel('Parent',hdl_gui,'Units','normalized','Position',[0.83 0.02 0.14 0.13]);
    check_oasis = uicontrol('Parent',panel_oasis,'style','check','Units','normalized',...
        'string','OASIS','Fontsize',11,'Position',[0.02 0.73 0.8 0.2],...
        'Value',1,'Callback',@check_oasis_callback);     
    text_fluorImpulse = uicontrol('Parent',panel_oasis,'style','text','Units','normalized',...
        'string','Fluor impulse','Enable','on','Fontsize',11); 
        set(text_fluorImpulse,'Position',[0.14 0.41 0.5 0.27]);
    edit_fluorImpulse = uicontrol('Parent',panel_oasis,'style','edit','Units','normalized',...
        'string','0.94','Enable','on','Callback',@edit_oasis_callback); 
        set(edit_fluorImpulse,'Position',[0.68 0.42 0.28 0.28]);
    text_sparsPenalty = uicontrol('Parent',panel_oasis,'style','text','Units','normalized',...
        'string','Sparsity penalty','Enable','on','Fontsize',11); 
        set(text_sparsPenalty,'Position',[0.15 0.06 0.5 0.27]);
    edit_sparsPenalty = uicontrol('Parent',panel_oasis,'style','edit','Units','normalized',...
        'string','2.4','Enable','on','Callback',@edit_oasis_callback); 
        set(edit_sparsPenalty,'Position',[0.68 0.06 0.28 0.28]);
    
        
    %% GUI Subfunctions
    function displayNumcells(Numcells)
        str = sprintf('Cell number: 1 to %g', Numcells);
        set(text_cellNum,'string',str);
        if get(slider_cellNum,'Value') > Numcells
            set(slider_cellNum,'Value',Numcells);
            set(edit_cellNum,'String',num2str(Numcells));
        end

        set(slider_cellNum,'Max',Numcells);
        set(slider_cellNum,'SliderStep',[1/Numcells 1/Numcells]);
    end

    function displayROI(id)
        % Display neuron image overlaid with specified ROI mask
        curr_masks = getappdata(hdl_gui,'curr_masks');
        curr_Numcells = getappdata(hdl_gui,'curr_Numcells');
        axes(ax_masks);
        imagesc(mean_R); axis off; colormap(jet);
        hold on
        if get(check_allROI,'Value') == 1
            for j = 1:curr_Numcells
                outline = bwboundaries(curr_masks(:,:,j));
                if size(outline,1) > 0
                    trace = outline{1};
                    plot(trace(:,2),trace(:,1),'w','Linewidth',2);
                end
            end
        end
        outline = bwboundaries(curr_masks(:,:,id));
        if size(outline,1) > 0
            trace = outline{1};
            plot(trace(:,2),trace(:,1),'g','Linewidth',3);
        end
        hold off
    end

    function displayROIAreaSNR(id)
        curr_maskArea = getappdata(hdl_gui,'curr_maskArea');
        area = curr_maskArea(id);
        str1 = sprintf('Area: %g pixels', round(area,1));
            set(text_ROIarea,'string',str1);
        curr_maskSNR = getappdata(hdl_gui,'curr_maskSNR');
        snr = curr_maskSNR(id);
        str2 = sprintf('Noise: %g', round(snr,2));
            set(text_ROIsnr,'string',str2);
    end

    function showTS(id)
        % Ca time series
        curr_tsG = getappdata(hdl_gui,'curr_tsG');
        axes(ax_ts);
        plot(curr_tsG(id,:)); axis([0 7500 0 Fmax]);

        % delta R/R
        curr_R = getappdata(hdl_gui,'curr_R');
        axes(ax_R);
        plot(curr_R(id,:)); axis([0 7500 Rmin Rmax]);
    
        % Spikes
        curr_spikes = getappdata(hdl_gui,'curr_spikes');
        axes(ax_spikes);
        plot(curr_spikes(id,:)); axis([0 7500 0 spikemax]);
    end

    function edit_cellNum_callback(varargin)
        id = str2double(get(edit_cellNum,'string'));
        curr_Numcells = getappdata(hdl_gui,'curr_Numcells');
        
        if isnan(id)
            set(edit_cellNum,'string','1');
            set(slider_cellNum,'Value',1);
            displayROIAreaSNR(area1);
            showTS(1);
        else
            if id < 1 
                set(edit_cellNum,'string','1');
                set(slider_cellNum,'Value',1);
                displayROIAreaSNR(1);
                showTS(1);
            elseif id > curr_Numcells
                set(edit_cellNum,'string',num2str(curr_Numcells));
                set(slider_cellNum,'Value',curr_Numcells);
                displayROIAreaSNR(curr_Numcells);
                showTS(curr_Numcells);
            else
                displayROI(id);
                displayROIAreaSNR(id);
                showTS(id);
                set(slider_cellNum,'Value',id);
            end
        end
    end

    function slider_cellNum_callback(varargin)
        id = round(get(slider_cellNum,'Value'));
        displayROI(id);
        showTS(id);
        displayROIAreaSNR(id);
        set(edit_cellNum,'string',num2str(id));
    end

    function check_allROI_callback(varargin)
        id = str2double(get(edit_cellNum,'string'));
        displayROI(id);    
    end

    function AreaSNRSatThr_callback(varargin)
        yn_areaThr = get(check_areaThr,'Value');
        if yn_areaThr
            set(edit_areaThr,'Enable','on');
        else
            set(edit_areaThr,'Enable','off');
        end
        yn_snrThr = get(check_snrThr,'Value');
        if yn_snrThr
            set(edit_snrThr,'Enable','on');
        else
            set(edit_snrThr,'Enable','off');
        end
        yn_satThr = get(check_satThr,'Value');
        if yn_satThr
            set(edit_satThr,'Enable','on');
        else
            set(edit_satThr,'Enable','off');
        end

        areaThr  = str2double(get(edit_areaThr,'string'));
        snrThr   = str2double(get(edit_snrThr,'string'));
        satValue = str2double(get(edit_satThr,'string'));
        
        if yn_areaThr 
            i_area = find(maskArea>=areaThr);
        else
            i_area = 1:Numcells;
        end
        if yn_snrThr
            i_snr = find(maskSNR<=snrThr);
        else
            i_snr = 1:Numcells;
        end
        if yn_satThr
            Fsat = (tsG >= satThr*satValue);
            i_sat = find(mean(Fsat,ndims(tsG))<satTime);
        else
            i_sat = 1:Numcells;
        end
        
        ind = intersect(intersect(i_area,i_snr),i_sat);
        displayNumcells(length(ind));
        
        spikes_nnd = getappdata(hdl_gui,'ALLspikes_nnd');
        spikes_oasis = getappdata(hdl_gui,'ALLspikes_oasis');
        
        setappdata(hdl_gui,'curr_ind',ind);
        setappdata(hdl_gui,'curr_Numcells',length(ind));
        setappdata(hdl_gui,'curr_masks',masks(:,:,ind));
        setappdata(hdl_gui,'curr_tsG',tsG(ind,:));
        setappdata(hdl_gui,'curr_R',R(ind,:));
        setappdata(hdl_gui,'curr_spikes_nnd',spikes_nnd(ind,:));
        setappdata(hdl_gui,'curr_spikes_oasis',spikes_oasis(ind,:));
        setappdata(hdl_gui,'curr_maskArea',maskArea(ind));
        setappdata(hdl_gui,'curr_maskSNR',maskSNR(ind));        
        
        edit_cellNum_callback;
        check_nnd_callback;
    end


    function check_nnd_callback(varargin)
        % NND is checked
        if get(check_nnd,'Value') == 1
            set(check_oasis,'Value',0);
            set(edit_fluorImpulse,'Enable','off');
            set(edit_sparsPenalty,'Enable','off');
            % set(edit_Rthr,'Enable','on');
            set(text_fluorImpulse,'Enable','off');
            set(text_sparsPenalty,'Enable','off');
            set(text_Rthr,'Enable','on');
            
            curr_spikes_nnd = getappdata(hdl_gui,'curr_spikes_nnd');
            setappdata(hdl_gui,'curr_spikes',curr_spikes_nnd);
            id = str2double(get(edit_cellNum,'string'));
            plot(ax_spikes, curr_spikes_nnd(id,:)); 
        else
            % NND is unchecked
            set(check_oasis,'Value',1);
            set(edit_fluorImpulse,'Enable','on');
            set(edit_sparsPenalty,'Enable','on');
            set(edit_Rthr,'Enable','off');
            set(text_fluorImpulse,'Enable','on');
            set(text_sparsPenalty,'Enable','on');
            set(text_Rthr,'Enable','off');
            
            curr_spikes_oasis = getappdata(hdl_gui,'curr_spikes_oasis');
            setappdata(hdl_gui,'curr_spikes',curr_spikes_oasis);
            id = str2double(get(edit_cellNum,'string'));
            plot(ax_spikes, curr_spikes_oasis(id,:)); 
        end
        
    end

    function edit_Rthr_callback(varargin)
        Rthr = str2double(get(edit_Rthr,'string'));
        % Recalculate spikes_nnd
        spikes_nnd = nndORoasis(R, 1, 0.91, 2.4);
        ind = getappdata(hdl_gui,'curr_ind');
        curr_spikes_nnd = spikes_nnd(ind,:);
        id = str2double(get(edit_cellNum,'string'));
        plot(ax_spikes, curr_spikes_nnd(id,:)); 
        
        setappdata(hdl_gui,'ALLspikes_nnd',spikes_nnd);
        setappdata(hdl_gui,'curr_spikes_nnd',curr_spikes_nnd);
        setappdata(hdl_gui,'curr_spikes',curr_spikes_nnd);
    end

    function check_oasis_callback(varargin)
        % OASIS is checked
        if get(check_oasis,'Value') == 1
            set(check_nnd,'Value',0);
            set(edit_Rthr,'Enable','off');
            set(edit_fluorImpulse,'Enable','on');
            set(edit_sparsPenalty,'Enable','on');
            set(text_Rthr,'Enable','off');
            set(text_fluorImpulse,'Enable','on');
            set(text_sparsPenalty,'Enable','on');
            
            curr_spikes_oasis = getappdata(hdl_gui,'curr_spikes_oasis');
            setappdata(hdl_gui,'curr_spikes',curr_spikes_oasis);
            id = str2double(get(edit_cellNum,'string'));
            plot(ax_spikes, curr_spikes_oasis(id,:)); 
        else
            % OASIS is unchecked
            set(check_nnd,'Value',1);
            set(edit_Rthr,'Enable','on');
            set(edit_fluorImpulse,'Enable','off');
            set(edit_sparsPenalty,'Enable','off');
            set(text_Rthr,'Enable','on');
            set(text_fluorImpulse,'Enable','off');
            set(text_sparsPenalty,'Enable','off');
            
            curr_spikes_nnd = getappdata(hdl_gui,'curr_spikes_nnd');
            setappdata(hdl_gui,'curr_spikes',curr_spikes_nnd);
            id = str2double(get(edit_cellNum,'string'));
            plot(ax_spikes, curr_spikes_nnd(id,:)); 
        end
    end

    function edit_oasis_callback(varargin)
        g1 = str2double(get(edit_fluorImpulse,'string'));
        lam1 = str2double(get(edit_sparsPenalty,'string'));
        
        % Recalculate spikes_oasis
        spikes_oasis = nndORoasis(R, 2, g1, lam1);
        ind = getappdata(hdl_gui,'curr_ind');
        curr_spikes_oasis = spikes_oasis(ind,:);
        id = str2double(get(edit_cellNum,'string'));
        plot(ax_spikes, curr_spikes_oasis(id,:)); 
        
        setappdata(hdl_gui,'ALLspikes_oasis',spikes_oasis);
        setappdata(hdl_gui,'curr_spikes_oasis',curr_spikes_oasis);
        setappdata(hdl_gui,'curr_spikes',curr_spikes_oasis);
    end
        
end