function GUI_view_mouseData( file )

% Written by Ann Go
%
% Script for viewing data for an individual video. 

%% Load file
[data_locn,~,err] = load_neuroSEEmodules(false);

if ~isempty(err)
   beep; errordlg(err,'ERROR');
    return
end

filedir = [data_locn 'Data/' file(1:8) '/Processed/' file '/'];
matfiles = dir(fullfile(filedir,['*.','mat']));
if numel(matfiles) > 0 
    [fileN,OK] = listdlg('PromptString','Select a file:','SelectionMode','single',...
                        'ListSize',[330 160],'ListString',{matfiles.name});
else
    beep; errordlg('No .mat file for processed data exists.','ERROR');
    return
end

if OK
    datafile = matfiles(fileN).name;
    data = load(fullfile(filedir, datafile));
    params = data.params;
    masks = data.masks;
        Ncells = size(masks,3);
    spikes = data.spikes;
    downData = data.downData;
    activeData = data.activeData;
    occMap = data.occMap;
    
    if isfield(data,'hist'); hist = data.hist; end
    if isfield(data,'asd'); asd = data.asd; end
    if isfield(data,'sorted_placeMap'); sort_pfMap = data.sorted_placeMap; end
    if isfield(data,'pfMap_sm'); pfMap_sm = data.pfMap_sm; end
    if isfield(data,'corr_image'); corr_image = data.corr_image; end
    if isfield(data,'mean_imratio'); mean_imratio = data.mean_imratio; end
    if isfield(data,'mcorr_output'); mcorr_output = data.mcorr_output; end
    if isfield(data,'ddf_f'); ddf_f = data.ddf_f; end
    if isfield(data,'df_f'); df_f = data.df_f; end
    if isfield(data,'R'); R = data.R; end
    if isfield(data,'dtsG'); dtsG = data.dtsG; end
    if isfield(data,'tsG'); tsG = data.tsG; end
    
    if isfield(params,'methods')
        mcorr_method = params.methods.mcorr_method;
        segment_method = params.methods.segment_method;
        dofissa = params.methods.dofissa;
        Nepochs = params.PFmap.Nepochs;
        Nbins = params.PFmap.Nbins;
        Vthr = params.PFmap.Vthr;
    else
        mcorr_method = 'fftRigid';
        segment_method = 'ABLE';
        dofissa = false;
        Nepochs = params.Nepochs;
        Nbins = params.Nbins;
        Vthr = params.Vthr;
    end
    if dofissa
        str_fissa = 'FISSA';
    else
        str_fissa = 'no FISSA';
    end
        
    %% GUI variables    
    mode_dim = params.mode_dim;
    if isnumeric(mode_dim)
        if mode_dim == 1
            mode_dim = '1D';
        else
            mode_dim = '2D';
        end
    end
    FOV = params.FOV;

    if exist('corr_image','var')
        im = corr_image; im_str = 'Correlation image';
    elseif exist('mean_imratio','var')
        im = mean_imratio; im_str = 'Green/red ratio image';
    elseif exist('mcorr_output','var')
        im = mcorr_output.green.meanregframe; im_str = 'Green channel image';
    end
    
    if exist('dtsG','var')
        C1 = dtsG;
    else
        C1 = tsG;
    end
    
    if exist('ddf_f','var')
        C2 = ddf_f; C2_str = 'dF/F';
    elseif exist('df_f','var')
        C2 = df_f; C2_str = 'dF/F';
    else
        C2 = R; C2_str = 'dR/R';
    end
    
    % all ROIs
    for j = 1:Ncells
        outline{:,:,j} = bwboundaries(masks(:,:,j));    % boundary of each ROI
        c = regionprops(masks(:,:,j),'centroid');
        centroids(j,:) = c.Centroid;                    % centroid of each ROI
    end
    
    % pcs from hist   
    if exist('hist','var') && strcmpi(mode_dim,'1D')
        masksH = masks(:,:,hist.pcIdx);
        Nhist = numel(hist.pcIdx);
            for j = 1:Nhist
                outlineH{:,:,j} = bwboundaries(masksH(:,:,j));    % boundary of each ROI
                c = regionprops(masksH(:,:,j),'centroid');
                centroidsH(j,:) = c.Centroid;                    % centroid of each ROI
            end
        activespkH = activeData.spikes_hist;
        
        masksA = masks(:,:,asd.pcIdx);
        Nasd = numel(asd.pcIdx);
            for j = 1:Nasd
                outlineA{:,:,j} = bwboundaries(masksA(:,:,j));    % boundary of each ROI
                c = regionprops(masksA(:,:,j),'centroid');
                centroidsA(j,:) = c.Centroid;                    % centroid of each ROI
            end
        activespkA = activeData.spikes_asd;
    else
        activespk = activeData.spikes;
    end
    
    downx = downData.x;
    downy = downData.y;
    downphi = downData.phi;
        dphi = [0; diff(downphi)];
        iLoop = [0; find( abs(dphi) > 190 ); length(downphi)];
    downt = downData.t;
        T = downt(end);
    activex = activeData.x;
    activey = activeData.y;
    pix2um = FOV/size(im,1);


    %% GUI structures

    % GUI figure
    str = [file ': ' mcorr_method ', ' segment_method ', ' str_fissa];
    hdl_gui = figure('MenuBar','none','Name',str,'NumberTitle','off','Resize','off',...
        'Position',[1087 648 1290 1000]); 

    % 2P image
    check_allROI = uicontrol('Parent',hdl_gui,'style','check','Units','pixels',...
        'string','Show all ROIs','Fontsize',12,'Position',[110 480 123 23],...
        'Value',0,'Callback',@check_allROI_callback);
    check_showLabels = uicontrol('Parent',hdl_gui,'style','check','Units','pixels',...
        'string','Show labels','Fontsize',12,'Position',[280 480 103 23],...
        'Enable','off','Value',0,'Callback',@check_showLabels_callback);
    ax_2P = axes('Units','pixels','Position',[50 520 450 450],'Xlim',[0 1],'Ylim',[0  1],'Box','off',...
        'Visible','on','clipping','off'); 
        
    % PF maps 
    ax_pfH = axes('Units','pixels','Position',[580 660 310 310],'Xlim',[0 1],'Ylim',[0  1],'Box','on',...
        'Visible','on','clipping','off'); 
        if exist('hist','var')
            if strcmpi(mode_dim,'1D')
                imagesc(hist.sort_normpfMap_sm(:,:,1)); 
                xticks([1 Nbins/2 Nbins]); xticklabels([1 50 100]); xlabel('Position (cm)');
                yticks([1 Nhist]); ylabel('Cell #');
            else
                imagesc(hist.pfMap_sm(:,:,1));
                axis off;
            end
        else
            if strcmpi(mode_dim,'1D')
                imagesc(sort_pfMap(:,:,1));
                xticks([1 Nbins/2 Nbins]); xticklabels([1 50 100]); xlabel('Position (cm)');
                yticks([1 size(sort_pfMap,1)]); ylabel('Cell #');
            else
                imagesc(pfMap_sm(:,:,1));
                axis off;
            end
        end
        title('HIST','Fontweight','normal','Fontsize',12);
    
    if exist('asd','var')    
        ax_pfA = axes('Units','pixels','Position',[930 660 310 310],'Xlim',[0 1],'Ylim',[0  1],'Box','on',...
            'Visible','on','clipping','off'); 
            if strcmpi(mode_dim,'1D')
                imagesc(asd.sort_normpfMap(:,:,1)); 
                xticks([1 Nbins/2 Nbins]); xticklabels([1 50 100]); xlabel('Position (cm)');
                yticks([1 Nasd]); 
            else
                imagesc(asd.pfMap(:,:,1));
                axis off;
            end
            title('ASD','Fontweight','normal','Fontsize',12);
    end
    
    % PF epoch control 
    if Nepochs > 1
        %label_Nepochs = 
        uicontrol('Parent',hdl_gui,'style','text','Units','pixels','string','Nepoch:',...
            'FontWeight','bold','Fontsize',12,'Position',[600 590 50 25]); 
        slider_Nepoch = uicontrol('Parent',hdl_gui,'style','slider','Units','pixels',...
            'Min',1,'Max',Nepochs,'Value',1,'SliderStep',[1/Nepochs 1/Nepochs],...
            'Position',[660 600 85 15],'Callback',@slider_Nepoch_callback); 
    end

    % Time series for specified cell
    % Raw calcium
    ax_rawCa = axes('Units','pixels','Position',[75 375 890 70],'Xlim',[0 1],'Ylim',[0  1],'Box','off',...
        'Visible','on','clipping','off'); 

    % Smoothened dR/R
    ax_df_f = axes('Units','pixels','Position',[75 285 890 70],'Xlim',[0 1],'Ylim',[0  1],'Box','off',...
        'Visible','on','clipping','off'); 

    % Spikes
    ax_spikes = axes('Units','pixels','Position',[75 195 890 70],'Xlim',[0 1],'Ylim',[0  1],'Box','off',...
        'Visible','on','clipping','off'); 
   
    % Trajectory and spike locations
    ax_TrajSpikeLoc = axes('Units','pixels','Position',[1000 290 220 220],'Xlim',[0 1],'Ylim',[0  1],'Box','off',...
        'Visible','on','clipping','off');

    % Settings panel
    panel_settings = uipanel('Parent',hdl_gui,'Units','pixels','Position',[50 23 230 135]);
    
    % Cell number
    if ~exist('masksH','var')
        str = sprintf('Cell number (1 to %g):', Ncells);
    else
        str = sprintf('Cell number (1 to %g):', Nhist);
    end
    text_cellNum = uicontrol('Parent',panel_settings,'style','text','Units','pixels','string',str,'Fontsize',16,...
        'FontWeight','bold','Position',[10 98 171 23]); 
    edit_cellNum = uicontrol('Parent',panel_settings,'style','edit','Units','pixels','string','1','Fontsize',14,...
        'Position',[14 63 45 26],'Callback',@edit_cellNum_callback); 
    if strcmpi(mode_dim,'1D') && exist('hist','var')
        slider_cellNum = uicontrol('Parent',panel_settings,'style','slider','Units','pixels','Position',[72 60 101 22],...
        'Min',1,'Max',Nhist,'Value',1,'SliderStep',[1/Nhist 1/Nhist],'Callback',@slider_cellNum_callback); 
    else
        slider_cellNum = uicontrol('Parent',panel_settings,'style','slider','Units','pixels','Position',[72 60 101 22],...
        'Min',1,'Max',Ncells,'Value',1,'SliderStep',[1/Ncells 1/Ncells],'Callback',@slider_cellNum_callback); 
    end
    
    
    % Other settings: Nepochs, Nbins, PFsmoothFac, Vthr
    %text_Nepochs = 
    uicontrol('Parent',panel_settings,'style','text','Units','pixels','string',['No. of epochs: ' num2str(Nepochs)],...
        'Enable','on','HorizontalAlignment','left','Fontsize',12,'Position',[10 25 110 23]);
    %text_Vthr = 
    uicontrol('Parent',panel_settings,'style','text','Units','pixels','string',['Speed threshold (mm/s): ' num2str(Vthr)],...
        'Enable','on','HorizontalAlignment','left','Fontsize',12,'Position',[10 5 160 23]);

    % Cell image
    ax_cellIm = axes('Units','pixels','Position',[825 35 120 120],'Xlim',[0 1],'Ylim',[0  1],'Box','on',...
        'Visible','on','clipping','off');         
        
    if strcmpi(mode_dim,'1D')
        % Occupancy map
        ax_occ = axes('Units','pixels','Position',[580 510 310 40],'Xlim',[0 1],'Ylim',[0  1],'Box','on',...
            'Visible','on','clipping','off'); 
        imagesc(occMap); % colormap(ax_occ,hsv);
        xticks([1 Nbins/2 Nbins]); xticklabels([1 50 100]); xlabel('Position (cm)');
        yticks([]); 
        title('Occupancy','Fontweight','normal','Fontsize',12);
        
        if exist('hist','var')
            % switch: hist, asd
            set(panel_settings,'Position',[50 23 340 135]);
            uicontrol('Parent',panel_settings,'style','text','Units','pixels','string','Show data for:','Fontsize',16,...
                'FontWeight','bold','Position',[208 98 113 23]); 
            check_hist = uicontrol('Parent',panel_settings,'style','checkbox','Units','pixels','string','HIST',...
                'Enable','on','Value',1,'Fontsize',12,'Position',[217 65 83 23],'Callback',@check_hist_callback);
            check_asd = uicontrol('Parent',panel_settings,'style','checkbox','Units','pixels','string','ASD',...
                'Enable','on','Value',0,'Fontsize',12,'Position',[217 37 83 23],'Callback',@check_asd_callback);

            % Per trial spike map
            ax_TrialSpk = axes('Units','pixels','Position',[1020 45 200 200],'Xlim',[0 1],'Ylim',[0  1],'Box','on',...
            'Visible','on','clipping','off'); 
            displayTrialSpk(1);
        end
    end
    
    displayROI(1);
    showTS(1);
    displayTrajSpikeLoc(1);
    displayCellIm(1);

end

%% GUI subfunctions
function displayROI(id)
    % Display neuron image overlaid with specified ROI mask
    axes(ax_2P);
    imagesc(im); 
    yticks(round(size(im,2)/4:size(im,2)/4:size(im,2))); 
    yticklabels(round(pix2um*(size(im,2)/4:size(im,2)/4:size(im,2))));
    xticks([]); ylabel('Size (um)','Fontsize',11);
    title(im_str,'Fontweight','normal','Fontsize',12);
    % colormap(ax_2P,jet);
    colormap(ax_2P,gray);
    hold on
    if get(check_allROI,'Value') == 1
        % Outline all ROIs 
        for j = 1:Ncells
            plot(outline{1,1,j}{1}(:,2),outline{1,1,j}{1}(:,1),'w','Linewidth',1);
        end
    end

    if strcmpi(mode_dim,'2D') || ~exist('hist','var')
        % Outline ROI of current cell in green
        plot(outline{1,1,id}{1}(:,2),outline{1,1,id}{1}(:,1),'g','Linewidth',3);
        hold off
    else
        if get(check_hist,'Value') == 1
            % Outline ROIs of all pcs
            for j = 1:Nhist
                plot(outlineH{1,1,j}{1}(:,2),outlineH{1,1,j}{1}(:,1),'b','Linewidth',1);
            end
            % Outline ROI of current cell in green
            plot(outlineH{1,1,id}{1}(:,2),outlineH{1,1,id}{1}(:,1),'g','Linewidth',3);
            hold off
        else
            % Outline ROIs of all pcs
            for j = 1:Nasd
                plot(outlineA{1,1,j}{1}(:,2),outlineA{1,1,j}{1}(:,1),'b','Linewidth',1);
            end
            % Outline ROI of current cell in green
            plot(outlineA{1,1,id}{1}(:,2),outlineA{1,1,id}{1}(:,1),'g','Linewidth',3);
            hold off
        end
    end
end

function showTS(id)
    if ~exist('masksH','var')
        idx = id;
    else
        if get(check_hist,'Value') == 1
            idx = hist.pcIdx(id);
        else
            idx = asd.pcIdx(id);
        end
    end
    
    % Ca time series
    plot(ax_rawCa, downt, C1(idx,:), 'Color', 'b'); axis(ax_rawCa,[0 T 0 inf]); 
    ylabel(ax_rawCa,'F','Fontweight','bold','Fontsize',12); xticks(ax_rawCa,[]);

    % df/f
    cell_df = C2(idx,:);
    axes(ax_df_f);
    for i = 1:length(iLoop)-1
        if mod(i,2) == 0
            color = 'r';
        else
            color = 'b';
        end
        plot(downt(iLoop(i)+1:iLoop(i+1)), cell_df(iLoop(i)+1:iLoop(i+1)), 'Color',color);
        axis(ax_df_f,[0 T 0 inf]); hold on
    end
    hold off
    ylabel(ax_df_f,C2_str,'Fontweight','bold','Fontsize',12); xticks(ax_df_f,[]);

    % Spikes
    plot(ax_spikes, downt, spikes(idx,:), 'Color', 'b'); axis(ax_spikes,[0 T 0 inf]);
    ylabel(ax_spikes,'Spikes','Fontweight','bold','Fontsize',12); xlabel(ax_spikes,'Time (s)','Fontweight','bold');
end

function displayTrajSpikeLoc(id)
    axes(ax_TrajSpikeLoc);
    plot(downx,downy,'Color',[0.5 0.5 0.5],'linewidth',3); axis off; axis square;
    str = sprintf('Cell %g', id);
    title(str,'Fontweight','normal','Fontsize',12);
    hold on
    if strcmpi(mode_dim,'2D') || ~exist('masksH','var')
        ind = find(activespk(id,:)>0);
        scatter(activex(ind),activey(ind),40,activespk(id,ind),'filled'); 
    else
        if get(check_hist,'Value') == 1
            ind = find(activespkH(id,:)>0);
            scatter(activex(ind),activey(ind),40,activespkH(id,ind),'filled'); 
        else
            ind = find(activespkA(id,:)>0);
            scatter(activex(ind),activey(ind),40,activespkA(id,ind),'filled'); 
        end
    end
    c = colorbar('Location','east','Units','pixels','Position',[1230 290 16 205]);
    title(c,'Spks/s');
    hold off
end

function displayCellIm(id)
    if strcmpi(mode_dim,'2D') || ~exist('hist','var')
        bounds = outline{1,1,id}{1};
    else
        if get(check_hist,'Value') == 1
            bounds = outlineH{1,1,id}{1};
        else
            bounds = outlineA{1,1,id}{1};
        end
    end
    pad = 10;
    xmin = min(bounds(:,2))-pad;
        if xmin < 1
            xmin = 1;
        end
    ymin = min(bounds(:,1))-pad; 
        if ymin < 1
            ymin = 1;
        end
    xmax = max(bounds(:,2))+pad;
        if xmax > size(im,2)
            xmax = size(im,2);
        end
    ymax = max(bounds(:,1))+pad;
        if ymax > size(im,1)
            ymax = size(im,1);
        end
    xrange = xmax - xmin;
    yrange = ymax - ymin;
    if xrange < yrange
        xyrange = xrange;
    else 
        xyrange = yrange;
    end
    
    im_cell = imcrop(im,[xmin ymin xyrange xyrange]);
    axes(ax_cellIm);
    imagesc(im_cell); % colormap(ax_cellIm,hsv); 
    yticks(round(xyrange/2:xyrange/2:xyrange)); yticklabels(round(pix2um*(xyrange/2:xyrange/2:xyrange)));
    ylabel('Size (um)');
    xticks([]); 
    title(['Cell ' num2str(id)],'Fontweight','normal','Fontsize',12);
end

function displayTrialSpk(id)
    Ntrials = size(hist.normspkMap_pertrial,1);
    axes(ax_TrialSpk);
    if get(check_hist,'Value') == 1
        imagesc(hist.normspkMap_pertrial(:,:,id));
    else
        imagesc(asd.normspkMap_pertrial(:,:,id));
    end
    % colormap(ax_TrialSpk,hsv);
    yticks(1:Ntrials:Ntrials); yticklabels([1,Ntrials]); ylabel('Trial #');
    xticks([1 15 30]); xticklabels([1 50 100]); xlabel('Position (cm)');
    title(['Cell ' num2str(id)],'fontsize',12);
end


%% GUI callback functions
function check_allROI_callback(varargin)
    % Disable showing labels when all ROIs are not displayed
    id = str2double(get(edit_cellNum,'string'));
    displayROI(id);
    if get(check_allROI,'Value') == 0
        set(check_showLabels,'Enable','off');
    else
        set(check_showLabels,'Enable','on');
        check_showLabels_callback;
    end
end

function check_showLabels_callback(varargin)
    if get(check_showLabels,'Value') == 1
        if strcmpi(mode_dim,'2D') || ~exist('masksH','var')
            for j = 1:Ncells
                text(centroids(j,1),centroids(j,2),num2str(j),'Color','w');
            end
        else
            if get(check_hist,'Value') == 1
                for j = 1:Nhist
                    text(centroidsH(j,1),centroidsH(j,2),num2str(j),'Color','w');
                end
            else
                for j = 1:Nasd
                    text(centroidsA(j,1),centroidsA(j,2),num2str(j),'Color','w');
                end
            end
        end
    else
        id = str2double(get(edit_cellNum,'string'));
        displayROI(id);    
    end
end

function edit_cellNum_callback(varargin)
    id = str2double(get(edit_cellNum,'string'));
    if strcmpi(mode_dim,'2D') || ~exist('masksH','var')
        n = size(masks,3);
    else
        if get(check_hist,'Value') == 1
            n = Nhist;
        else
            n = Nasd;
        end
    end

    if isnan(id)
        set(edit_cellNum,'string','1');
        set(slider_cellNum,'Value',1);
        displayROI(1);
        showTS(1);
        displayTrajSpikeLoc(1);
        displayCellIm(1);
        if strcmpi(mode_dim,'1D') && exist('hist','var') 
            displayTrialSpk(1); 
        end
    else
        if ~isinteger(id)
            id = round(id);
            set(edit_cellNum,'String',num2str(id));
        end
        if id < 1 
            set(edit_cellNum,'string','1');
            set(slider_cellNum,'Value',1);
            displayROI(1);
            showTS(1);
            displayTrajSpikeLoc(1);
            displayCellIm(1);
            if strcmpi(mode_dim,'1D') && exist('hist','var') 
                displayTrialSpk(1); 
            end
        elseif id > n
            set(edit_cellNum,'string',num2str(n));
            set(slider_cellNum,'Value',n);
            displayROI(n);
            showTS(n);
            displayTrajSpikeLoc(n);
            displayCellIm(n);
            if strcmpi(mode_dim,'1D') && exist('hist','var') 
                displayTrialSpk(n);     
            end
        else
            set(edit_cellNum,'String',num2str(id));
            set(slider_cellNum,'Value',id);
            displayROI(id);
            showTS(id);
            displayTrajSpikeLoc(id);
            displayCellIm(id);
            if strcmpi(mode_dim,'1D') && exist('hist','var') 
                displayTrialSpk(id);     
            end
        end
    end
end

function slider_cellNum_callback(varargin)
    id = round(get(slider_cellNum,'Value'));
    displayROI(id);
    showTS(id);
    displayTrajSpikeLoc(id);
    displayCellIm(id);
    if strcmpi(mode_dim,'1D') && exist('hist','var')
        displayTrialSpk(id); 
    end
    set(edit_cellNum,'string',num2str(id));
end

function check_hist_callback(varargin)
    id = str2double(get(edit_cellNum,'String'));
    if get(check_hist,'Value') == 1
        set(check_asd,'Value',0);
        set(text_cellNum,'String',['Cell number (1 to ' num2str(Nhist) ')']);
        displayROI(id);
    else
        set(check_hist,'Value',1);
        set(text_cellNum,'String',['Cell number (1 to ' num2str(Nasd) ')']);
        displayROI(id);
    end
end

function check_asd_callback(varargin)
    id = str2double(get(edit_cellNum,'String'));
    if get(check_asd,'Value') == 1
        set(check_hist,'Value',0);
        set(text_cellNum,'String',['Cell number (1 to ' num2str(Nasd) ')']);
        displayROI(id);
    else
        set(check_hist,'Value',1);
        set(text_cellNum,'String',['Cell number (1 to ' num2str(Nhist) ')']);
        displayROI(id);
    end
end

function slider_Nepoch_callback(varargin)
    e = round(get(slider_Nepoch,'Value'));
    
    axes(ax_pfH);
        imagesc(hist.sort_normpfMap_sm(:,:,e)); % colormap(ax_pfH,hsv);
        if strcmpi(mode_dim,'1D')
            xticks([1 Nbins/2 Nbins]); xticklabels([1 50 100]); xlabel('Position (cm)');
            yticks([1 Nhist]); ylabel('Cell #');
        else
            axis off;
        end
        title(['HIST: Epoch ' num2str(e)],'Fontweight','normal','Fontsize',12);
        
    axes(ax_pfA);
        imagesc(asd.sort_normpfMap(:,:,e)); % colormap(ax_pfA,hsv);
        if strcmpi(mode_dim,'1D')
            xticks([1 Nbins/2 Nbins]); xticklabels([1 50 100]); xlabel('Position (cm)');
            yticks([1 Nasd]); 
        else
            axis off;
        end
        title(['ASD: Epoch ' num2str(e)],'Fontweight','normal','Fontsize',12);
    
end

end