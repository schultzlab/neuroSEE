% Written by Ann Go
% 
% This GUI display mean 2P image and the ROI mask, raw Ca transient, ratiometric Ca
% transient (delta R / R) and the extracted spikes for each cell.
% User can choose cell and mode of spike extraction (NND vs OASIS).

function GUI_viewROIsSpikes_withSave( file, mean_imratio, masks, cell_tsG, cell_tsR, R, spikes )
    
    setAxisLimits = 0;
    Ncells = length(masks);
    tmasks = zeros(512,512,Ncells);
    for i = 1:Ncells
        mask = zeros(512,512);
        mask(masks{i}) = 1;
        tmasks(:,:,i) = mask;
    end

    % Create figure
    hfig_h = 620;
    hfig_w = 1250;
    hdl_gui = figure('MenuBar','none','Name','ROIs and spikes','NumberTitle','off',...%'Resize','off',...
    'Position',[1680-hfig_w,1050-hfig_h,hfig_w,hfig_h]);

    % Display cell number and allow user to change cell number either by
    % using a slider or an edit box
    str = sprintf('Cell number: 1 to %g', Ncells);
    text_cellNum = uicontrol('Parent',hdl_gui,'style','text','Units','pixels','string',str,'Fontsize',12); 
        set(text_cellNum,'Position',[24 32 96 51]);
    edit_cellNum = uicontrol('Parent',hdl_gui,'style','edit','Units','pixels','string','1','Callback',@edit_cellNum_callback); 
        set(edit_cellNum,'Position',[128 55 40 28],'Fontsize',12);
    slider_cellNum = uicontrol('Parent',hdl_gui,'style','slider','Units','pixels',...
        'Min',1,'Max',Ncells,'Value',1,'SliderStep',[1/Ncells 1/Ncells],'Callback',@slider_cellNum_callback); 
        set(slider_cellNum,'Position',[186 52 136 24]);

    % Allow option to show all ROIs
    check_allROI = uicontrol('Parent',hdl_gui,'style','check','Units','pixels',...
        'string','Show all ROIs','Fontsize',11,'Position',[25 566 102 23],...
        'Value',0,'Callback',@check_allROI_callback);
    
    % Show 2P image and overlay mask of specific cell
    ax_masks = axes('Units','pixels','Position',[25 107 450 450],'Xlim',[0 1],'Ylim',[0  1],'Box','off',...
                'Visible','off','clipping' , 'off'); 
    imagesc(mean_imratio); colormap(jet); axis off;
    hold on
    outline = bwboundaries(tmasks(:,:,1));
    trace = outline{1};
    plot(trace(:,2),trace(:,1),'g','Linewidth',3);
    hold off
    
    % Plot Ca time series for specific cell
    ax_tsG = axes('Units','pixels', 'Position',[534 498 690 81],'Xlim',[0 1],'Ylim',[0  1],'Box','off',...
                'Visible','off','clipping' , 'off'); 
    plot(cell_tsG(1,:)); title('Raw Ca^{2+} time series','Fontsize',12,'Fontweight','bold');
    if setAxisLimits, axis([0 7400 10000 60000]); end
    
    % Plot red channel time series for specific cell
    ax_tsR = axes('Units','pixels', 'Position',[534 380 690 81],'Xlim',[0 1],'Ylim',[0  1],'Box','off',...
                'Visible','off','clipping' , 'off'); 
    plot(cell_tsR(1,:)); title('mRuby signal','Fontsize',12,'Fontweight','bold');
    %axis([0 7400 4000 9000]);
    
    % Plot ratiometric Ca time series (delta R/R) for specific cell
    ax_R = axes('Units','pixels', 'Position',[536 260 690 81],'Xlim',[0 1],'Ylim',[0  1],'Box','off',...
                'Visible','off','clipping' , 'off'); 
    plot(R(1,:)); title('deltaR / R','Fontsize',12,'Fontweight','bold');
    if setAxisLimits, axis([0 7400 -0.2 1]); end
    
    % Plot spikes for specific cell. Allow user to choose between 
    ax_spikes = axes('Units','pixels', 'Position',[535 139 690 81],'Xlim',[0 1],'Ylim',[0  1],'Box','off',...
                'Visible','off','clipping' , 'off'); 
    plot(spikes(1,:)); title('Spikes','Fontsize',12,'Fontweight','bold');
    if setAxisLimits, axis([0 7400 -0.02 0.1]); end
    
    % Allow changing of parameters for R and spikes
    text_Rsmooth = uicontrol('Parent',hdl_gui,'style','text','Units','pixels','string','Rsmoothing','Fontsize',12); 
        set(text_Rsmooth,'Position',[554 72 69 20]);
    edit_Rsmooth = uicontrol('Parent',hdl_gui,'style','edit','Units','pixels','string','23','Callback',@edit_Rsmooth_callback); 
        set(edit_Rsmooth,'Position',[632 64 63 30],'Enable','off');

    text_g = uicontrol('Parent',hdl_gui,'style','text','Units','pixels','string','g','Fontsize',12); 
        set(text_g,'Position',[762 70 24 20]);
    edit_g = uicontrol('Parent',hdl_gui,'style','edit','Units','pixels','string','0.91','Callback',@edit_glambda_callback); 
        set(edit_g,'Position',[786 63 63 30]);

    text_lambda = uicontrol('Parent',hdl_gui,'style','text','Units','pixels','string','lambda','Fontsize',12); 
        set(text_lambda,'Position',[733 29 51 20]);
    edit_lambda = uicontrol('Parent',hdl_gui,'style','edit','Units','pixels','string','2.4','Callback',@edit_glambda_callback); 
        set(edit_lambda,'Position',[786 23 63 30]);
    
    % Allow saving of all variables
    uicontrol('Parent',hdl_gui,'style','pushbutton','Units','pixels',...
        'string','Reset','Fontsize',12,'Position',[939 49 85 33],...
        'Enable','off','Value',0,'Callback',@pbutton_reset_callback);
    uicontrol('Parent',hdl_gui,'style','pushbutton','Units','pixels',...
        'string','SAVE','Fontsize',12,'Position',[1075 47 125 44],...
        'Value',0,'Callback',@pbutton_save_callback);

    % Initialise GUI data
    setappdata(hdl_gui,'curr_R',R);
    setappdata(hdl_gui,'curr_spikes',spikes);
  
    Rsmooth = ones(1,Ncells)*23;
    g = ones(1,Ncells)*0.91;
    lambda = ones(1,Ncells)*2.4;
    setappdata(hdl_gui,'curr_Rsmooth',Rsmooth);
    setappdata(hdl_gui,'curr_g',g);
    setappdata(hdl_gui,'curr_lambda',lambda);
    
    % functions
    function displayROI(id)
        % Display neuron image overlaid with specified ROI mask
        axes(ax_masks);
        imagesc(mean_imratio); axis off; colormap(jet);
        hold on
        if get(check_allROI,'Value') == 1
            for j = 1:Ncells
                outline = bwboundaries(tmasks(:,:,id));
                if size(outline,1) > 0
                    trace = outline{1};
                    plot(trace(:,2),trace(:,1),'w','Linewidth',2);
                end
            end
        end
        outline = bwboundaries(tmasks(:,:,id));
        if size(outline,1) > 0
            trace = outline{1};
            plot(trace(:,2),trace(:,1),'g','Linewidth',3);
        end
        hold off
    end

    function plotTS(id)
        % Plot Ca time series for specific cell
        axes(ax_tsG); plot(cell_tsG(id,:)); title('Raw Ca^{2+} time series','Fontsize',12,'Fontweight','bold');
        if setAxisLimits, axis([0 7400 10000 60000]); end
        
        % Plot red channel time series for specific cell
        axes(ax_tsR); plot(cell_tsR(id,:)); title('mRuby signal','Fontsize',12,'Fontweight','bold');
        % axis([0 7400 4000 9000]);
        
        % Plot ratiometric Ca time series (delta R/R) for specific cell
        curr_R = getappdata(hdl_gui,'curr_R');
        axes(ax_R);
        plot(curr_R(id,:)); title('deltaR / R','Fontsize',12,'Fontweight','bold');
        if setAxisLimits, axis([0 7400 -0.2 1]); end
        
        % Plot spikes for specific cell. Allow user to choose between 
        curr_spikes = getappdata(hdl_gui,'curr_spikes');
        axes(ax_spikes);
        plot(curr_spikes(id,:)); title('Spikes','Fontsize',12,'Fontweight','bold');
        if setAxisLimits, axis([0 7400 -0.02 0.1]); end
    end

    function displayParams(id)
        curr_Rsmooth = getappdata(hdl_gui,'curr_Rsmooth');
            set(edit_Rsmooth,'string',num2str(curr_Rsmooth(id)));
        curr_g = getappdata(hdl_gui,'curr_g');
            set(edit_g,'string',num2str(curr_g(id)));
        curr_lambda = getappdata(hdl_gui,'curr_lambda');
            set(edit_lambda,'string',num2str(curr_lambda(id)));
    end

    function edit_cellNum_callback(varargin)
        id = str2double(get(edit_cellNum,'string'));
        
        if isnan(id)
            set(edit_cellNum,'string','1');
            set(slider_cellNum,'Value',1);
            id = 1;
        else
            if id < 1 
                set(edit_cellNum,'string','1');
                set(slider_cellNum,'Value',1);
                id = 1;
            elseif id > Ncells
                set(edit_cellNum,'string',num2str(Ncells));
                set(slider_cellNum,'Value',Ncells);
                id = Ncells;
            else
                set(slider_cellNum,'Value',id);
            end
        end
        
        displayROI(id);
        plotTS(id);
        displayParams(id)
    end

    function slider_cellNum_callback(varargin)
        id = round(get(slider_cellNum,'Value'));
        displayROI(id);
        plotTS(id);
        displayParams(id)
        set(edit_cellNum,'string',num2str(id));
    end

    function check_allROI_callback(varargin)
        id = str2double(get(edit_cellNum,'string'));
        displayROI(id);  
    end

%     function edit_Rsmooth_callback(varargin)
%         Rsmoothi = str2double(get(edit_Rsmooth,'string'));
%         id = round(get(slider_cellNum,'Value'));
%         
%         % Recalculate R
%         RR = ratiometric_Ca( cell_tsG, cell_tsR, Rsmoothi );
%         curr_R = getappdata(hdl_gui,'curr_R');
%         curr_R(id,:) = RR(id,:);
%         setappdata(hdl_gui,'curr_R',curr_R);
%         axes(ax_R); plot(curr_R(id,:)); 
%         edit_glambda_callback;
%     end
%     
    function edit_Rsmooth_callback(varargin)
        Rsmoothi = str2double(get(edit_Rsmooth,'string'));
        id = round(get(slider_cellNum,'Value'));
        
        % Recalculate R
        Ri = ratiometric_Ca( cell_tsG(id,:), cell_tsR(id,:), Rsmoothi );
        curr_R = getappdata(hdl_gui,'curr_R');
        curr_R(id,:) = Ri;
        curr_Rsmooth = getappdata(hdl_gui,'curr_Rsmooth');
        curr_Rsmooth(id) = Rsmoothi;
        setappdata(hdl_gui,'curr_R',curr_R);
        setappdata(hdl_gui,'curr_Rsmooth',curr_Rsmooth);
        cprintf('R: Success!\n')
        axes(ax_R); plot(curr_R(id,:)); title('deltaR / R','Fontsize',12,'Fontweight','bold');
        if setAxisLimits, axis([0 7400 -0.2 1]); end
        edit_glambda_callback;
    end

%     function edit_glambda_callback(varargin)
%         gi = str2double(get(edit_g,'string'));
%         lambdai = str2double(get(edit_lambda,'string'));
%         
%         id = round(get(slider_cellNum,'Value'));
%         curr_R = getappdata(hdl_gui,'curr_R');
%         
%         % Recalculate spike
%         spikeSS = nndORoasis(curr_R, 2, gi, lambdai);
%         
%         curr_g = getappdata(hdl_gui,'curr_g');
%         curr_lambda = getappdata(hdl_gui,'curr_lambda');
%         curr_spikes = getappdata(hdl_gui,'curr_spikes');
%         curr_g(id) = gi;
%         curr_lambda(id) = lambdai;
%         curr_spikes(id,:) = spikeSS(id,:);
%         
%         setappdata(hdl_gui,'curr_g',curr_g);
%         setappdata(hdl_gui,'curr_lambda',curr_lambda);
%         setappdata(hdl_gui,'curr_spikes',curr_spikes);
% 
%         axes(ax_spikes); plot(curr_spikes(id,:));  
%     end

    function edit_glambda_callback(varargin)
        gi = str2double(get(edit_g,'string'));
        lambdai = str2double(get(edit_lambda,'string'));
        
        id = round(get(slider_cellNum,'Value'));
%         curr_R = getappdata(hdl_gui,'curr_R');
        
        % Recalculate spike
        spikesi = nndORoasis(R(id,:), 2, gi, lambdai);
        
        curr_g = getappdata(hdl_gui,'curr_g');
        curr_lambda = getappdata(hdl_gui,'curr_lambda');
        curr_spikes = getappdata(hdl_gui,'curr_spikes');
        curr_g(id) = gi;
        curr_lambda(id) = lambdai;
        curr_spikes(id,:) = spikesi;
        
        setappdata(hdl_gui,'curr_g',curr_g);
        setappdata(hdl_gui,'curr_lambda',curr_lambda);
        setappdata(hdl_gui,'curr_spikes',curr_spikes);

        cprintf('Spikes: Success!\n')
        axes(ax_spikes); plot(curr_spikes(id,:)); title('Spikes','Fontsize',12,'Fontweight','bold'); 
        if setAxisLimits, axis([0 7400 -0.02 0.1]); end 
    end

    function pbutton_reset_callback(varargin)
        
    end

    function pbutton_save_callback(varargin)
%         R = getappdata(hdl_gui,'curr_R');
        spikes = getappdata(hdl_gui,'curr_spikes');
        Rsmooth = getappdata(hdl_gui,'curr_Rsmooth');
        g = getappdata(hdl_gui,'curr_g');
        lambda = getappdata(hdl_gui,'curr_lambda');
        save([file '_tsRspikes.mat'],'cell_tsG','cell_tsR','R','spikes','masks','mean_imratio','Rsmooth','g','lambda')
    end
end