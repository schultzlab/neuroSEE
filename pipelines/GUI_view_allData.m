% Written by Ann Go
%
% Script for viewing data for an individual video. Script asks user to open
% a mat file containing all data for the processed video and opens a gui

[filename, pathname, proceed] = uigetfile('*.mat', 'Choose data to view');

if proceed
    load(fullfile(pathname, filename))
else
    return
end

try
    mean_imG = mcorr_output.green.meanregframe;
    mean_imR = mcorr_output.red.meanregframe;
    FOV = params.FOV;
    Nbins = params.Nbins;

    displayDataGUI(file,mean_imG,mean_imR,mean_imratio,masks,tsG,R,spikes,downData,activeData,sorted_placeMap,FOV,Nbins);
catch
    beep;
    errordlg('Not a valid data file!');
    return
end

function displayDataGUI(file,mean_imG,mean_imR,mean_imratio,masks,tsG,R,spikes,downData,activeData,sorted_placeMap,FOV,Nbins)

addpath(genpath('../PF_mapping'));

%% GUI variables
Numcells = size(masks,3);                           % no. of cells
for j = 1:Numcells
    outline{:,:,j} = bwboundaries(masks(:,:,j));    % boundary of each ROI
    c = regionprops(masks(:,:,j),'centroid');
    centroids(j,:) = c.Centroid;                    % centroid of each ROI
end
downx = downData.x;
downy = downData.y;
downphi = downData.phi;
    dphi = [0; diff(downphi)];
    iLoop = [0; find(dphi<-150); length(downphi)];
downt = downData.t;
    T = downt(end);
activex = activeData.x;
activey = activeData.y;
activespikes = activeData.spikes;
pix2um = FOV/size(mean_imG,1);
degperbin = 360/Nbins;

% axis min and max
cmin = min(min(activespikes));                      % spike map colorbar
cmax = max(max(activespikes));
tsGmin = min(min(tsG));                             % raw Ca time series
tsGmax = 8000; %max(max(tsG));
Rmin = min(min(R));                                 % smoothened dR/R
Rmax = 2.0; %max(max(R));
spikesmin = min(min(spikes));                       % spikes
spikesmax = 0.5; %max(max(spikes));

%% GUI structures

% Create figure
str = sprintf('%s data',file);
hdl_gui = figure('MenuBar','none','Name',str,'NumberTitle','off','Resize','off',...
    'Position',[1087 648 1190 1040]); % for epochs, width = 1770

% 2P image
check_allROI = uicontrol('Parent',hdl_gui,'style','check','Units','pixels',...
    'string','Show all ROIs','Fontsize',12,'Position',[110 490 123 23],...
    'Value',0,'Callback',@check_allROI_callback);
check_showLabels = uicontrol('Parent',hdl_gui,'style','check','Units','pixels',...
    'string','Show labels','Fontsize',12,'Position',[280 490 103 23],...
    'Enable','off','Value',0,'Callback',@check_showLabels_callback);
ax_whole2P = axes('Units','pixels','Position',[50 516 500 500],'Xlim',[0 1],'Ylim',[0  1],'Box','off',...
    'Visible','off','clipping','off'); 
    
    displayROI(1);


% dR/R for all cells
check_showSubset = uicontrol('Parent',hdl_gui,'style','check','Units','pixels',...
    'string','Show in subsets','Fontsize',12,'Position',[575 490 123 23],...
    'Value',0,'Callback',@check_showSubset_callback);
edit_showSubset = uicontrol('Parent',hdl_gui,'style','edit','Units','pixels','string','4','Fontsize',14,...
    'Enable','off','Callback',@edit_showSubset_callback); 
    set(edit_showSubset,'Position',[690 485 45 26]);
Nsubset = str2double(get(edit_showSubset,'str'));
slider_showSubset = uicontrol('Parent',hdl_gui,'style','slider','Units','pixels','Enable','off',...
    'Min',1,'Max',Nsubset,'Value',1,'SliderStep',[1/Nsubset 1/Nsubset],'Callback',@displayAllR); 
    set(slider_showSubset,'Position',[745 485 101 22]);
ax_allR = axes('Units','pixels','Position',[565 516 310 500],'Xlim',[0 1],'Ylim',[0  1],'Box','on',...
    'Visible','off','clipping','off'); 
    
    displayAllR;
    
% Trajectory and spike locations
ax_allSpikeLoc = axes('Units','pixels','Position',[890 795 220 220],'Xlim',[0 1],'Ylim',[0  1],'Box','on',...
    'Visible','on','clipping','off');
    
    plot(downx,downy,'Color',[0.5 0.5 0.5],'linewidth',6); 
        title('Spike locations for all cells','Fontweight','normal','Fontsize',12); axis square; axis off; 
        hold on
        cmap = colormap(ax_allSpikeLoc,jet(Numcells));
        for j = 1:1:Numcells
            ind1 = find(activespikes(j,:)>0);
            scatter(activex(ind1),activey(ind1),0.75,cmap(j,:),'filled');
        end
        hold off

ax_cellSpikeLoc = axes('Units','pixels','Position',[890 550 220 220],'Xlim',[0 1],'Ylim',[0  1],'Box','off',...
    'Visible','off','clipping','off');

    displaySpikeLoc(1);

% Time series for specified cell
% Raw calcium
ax_rawCa = axes('Units','pixels','Position',[75 385 770 70],'Xlim',[0 1],'Ylim',[0  1],'Box','off',...
    'Visible','off','clipping','off'); 

% Smoothened dR/R
ax_dR = axes('Units','pixels','Position',[75 300 770 70],'Xlim',[0 1],'Ylim',[0  1],'Box','off',...
    'Visible','off','clipping','off'); 
    
% Spikes
ax_spikes = axes('Units','pixels','Position',[75 215 770 70],'Xlim',[0 1],'Ylim',[0  1],'Box','off',...
    'Visible','off','clipping','off'); 
    
    showTS(1);
    
% Settings panel
panel_settings = uipanel('Parent',hdl_gui,'Units','pixels','Position',[38 27 190 130]);
% Cell number
str = sprintf('Cell number (1 to %g)', Numcells);
text_cellNum = uicontrol('Parent',panel_settings,'style','text','Units','pixels','string',str,'Fontsize',16,...
    'FontWeight','bold'); 
    set(text_cellNum,'Position',[10 98 171 23]);
edit_cellNum = uicontrol('Parent',panel_settings,'style','edit','Units','pixels','string','1','Fontsize',14,...
    'Callback',@edit_cellNum_callback); 
    set(edit_cellNum,'Position',[14 63 45 26]);
slider_cellNum = uicontrol('Parent',panel_settings,'style','slider','Units','pixels',...
    'Min',1,'Max',Numcells,'Value',1,'SliderStep',[1/Numcells 1/Numcells],'Callback',@slider_cellNum_callback); 
    set(slider_cellNum,'Position',[72 60 101 22]);

% Other settings: Nepochs, Nbins, PFsmoothFac, Vthr
text_Nepochs = uicontrol('Parent',panel_settings,'style','text','Units','pixels','string','No. of epochs:',...
    'Enable','off','HorizontalAlignment','left','Fontsize',12);
    set(text_Nepochs,'Position',[10 25 82 23]);
text_Vthr = uicontrol('Parent',panel_settings,'style','text','Units','pixels','string','Speed threshold:',...
    'Enable','off','HorizontalAlignment','left','Fontsize',12);
    set(text_Vthr,'Position',[10 5 101 23]);

edit_Nepochs = uicontrol('Parent',panel_settings,'style','edit','Units','pixels','string','1','Fontsize',11,...
    'Enable','off','Callback',@edit_Nepochs_callback); 
    set(edit_Nepochs,'Position',[110 35 44 19]);
edit_Vthr = uicontrol('Parent',panel_settings,'style','edit','Units','pixels','string','1','Fontsize',11,...
    'Enable','off','Callback',@edit_Vthr_callback); 
    set(edit_Vthr,'Position',[110 10 44 19]);
    
% Cell images
% GCaMP image
ax_gcamp = axes('Units','pixels','Position',[295 45 110 110],'Xlim',[0 1],'Ylim',[0  1],'Box','on',...
    'Visible','on','clipping','off'); 

% mRuby image
ax_mruby = axes('Units','pixels','Position',[425 45 110 110],'Xlim',[0 1],'Ylim',[0  1],'Box','on',...
    'Visible','on','clipping','off'); 
    
    displayGcampMruby(1);

% PF maps for entire recording and for each cell
% ax_allPF = 
axes('Units','pixels','Position',[900 65 265 425],'Xlim',[0 1],'Ylim',[0  1],'Box','on',...
    'Visible','on','clipping','off'); 
    
    imagesc(sorted_placeMap); 
    xticks(Nbins/6:Nbins/6:Nbins); xticklabels(degperbin*(Nbins/6:Nbins/6:Nbins)); 
    xlabel('Position (degrees)'); yticks([]); ylabel('Cells');
    title('Normalised place field maps','Fontweight','normal','Fontsize',12);

ax_cellPF = axes('Units','pixels','Position',[580 60 265 15],'Xlim',[0 1],'Ylim',[0  1],'Box','on',...
    'Visible','off','clipping','off'); 
    
ax_ep1cellPF = axes('Units','pixels','Position',[580 80 265 15],'Xlim',[0 1],'Ylim',[0  1],'Box','on',...
    'Visible','off','clipping','off'); 

ax_ep2cellPF = axes('Units','pixels','Position',[580 100 265 15],'Xlim',[0 1],'Ylim',[0  1],'Box','on',...
    'Visible','off','clipping','off'); 

ax_ep3cellPF = axes('Units','pixels','Position',[580 120 265 15],'Xlim',[0 1],'Ylim',[0  1],'Box','on',...
    'Visible','off','clipping','off'); 

ax_ep4cellPF = axes('Units','pixels','Position',[580 140 265 15],'Xlim',[0 1],'Ylim',[0  1],'Box','on',...
    'Visible','off','clipping','off'); 

    displayCellPF(1);


% Epoch 1-4 panel
panel_epochs = uipanel('Parent',hdl_gui,'Units','pixels','Position',[1180 7 578 1027],'Visible','off');
% Epoch 1
ax_epoch1 = axes('Parent',panel_epochs,'Units','pixels','Position',[32 569 251 425],'Xlim',[0 1],'Ylim',[0  1],...
    'Box','on','Visible','on','clipping','off'); 
check_ep1overall = uicontrol('Parent',panel_epochs,'style','check','Units','pixels',...
    'string','overall','Fontsize',12,'Position',[105 515 72 23],...
    'Value',1,'Callback',@check_ep1overall_callback);
check_ep1own = uicontrol('Parent',panel_epochs,'style','check','Units','pixels',...
    'string','own','Fontsize',12,'Position',[182 515 60 23],...
    'Value',0,'Callback',@check_ep1own_callback);

% Epoch 2
ax_epoch2 = axes('Parent',panel_epochs,'Units','pixels','Position',[309 569 251 425],'Xlim',[0 1],'Ylim',[0  1],...
    'Box','on','Visible','on','clipping','off'); 
check_ep2overall = uicontrol('Parent',panel_epochs,'style','check','Units','pixels',...
    'string','overall','Fontsize',12,'Position',[387 515 72 23],...
    'Value',1,'Callback',@check_ep2overall_callback);
check_ep2own = uicontrol('Parent',panel_epochs,'style','check','Units','pixels',...
    'string','own','Fontsize',12,'Position',[464 515 60 23],...
    'Value',0,'Callback',@check_ep2own_callback);

% Epoch 3
ax_epoch3 = axes('Parent',panel_epochs,'Units','pixels','Position',[31 58 251 425],'Xlim',[0 1],'Ylim',[0  1],...
    'Box','on','Visible','on','clipping','off'); 
check_ep3overall = uicontrol('Parent',panel_epochs,'style','check','Units','pixels',...
    'string','overall','Fontsize',12,'Position',[105 3 72 23],...
    'Value',1,'Callback',@check_ep3overall_callback);
check_ep3own = uicontrol('Parent',panel_epochs,'style','check','Units','pixels',...
    'string','own','Fontsize',12,'Position',[182 3 60 23],...
    'Value',0,'Callback',@check_ep3own_callback);

% Epoch 4
text_epoch4 = uicontrol('Parent',panel_epochs,'style','text','Units','pixels','string','Epoch 4',...
    'Enable','on','Fontsize',12); 
    set(text_epoch4,'Position',[376 485 117 20]);
ax_epoch4 = axes('Parent',panel_epochs,'Units','pixels','Position',[309 59 251 425],'Xlim',[0 1],'Ylim',[0  1],...
    'Box','on','Visible','on','clipping','off'); 
text_ep4sorting = uicontrol('Parent',panel_epochs,'style','text','Units','pixels','string','Sorting:',...
    'Enable','on','Fontsize',12); 
    set(text_ep4sorting,'Position',[315 4 70 20]);
check_ep4overall = uicontrol('Parent',panel_epochs,'style','check','Units','pixels',...
    'string','overall','Fontsize',12,'Position',[387 4 72 23],...
    'Value',0,'Callback',@check_ep4overall_callback);
check_ep4own = uicontrol('Parent',panel_epochs,'style','check','Units','pixels',...
    'string','own','Fontsize',12,'Position',[464 4 60 23],...
    'Value',0,'Callback',@check_ep4own_callback);

%% GUI subfunctions
function displayROI(id)
    % Display neuron image overlaid with specified ROI mask
    axes(ax_whole2P);
    imagesc(mean_imratio); 
    yticks(round(size(mean_imG,2)/4:size(mean_imG,2)/4:size(mean_imG,2))); 
    yticklabels(round(pix2um*(size(mean_imG,2)/4:size(mean_imG,2)/4:size(mean_imG,2))));
    xticks([]); ylabel('Size (um)','Fontsize',11);
    title('mean(GCaMP) / mean(mRuby)','Fontweight','normal','Fontsize',12);
    colormap(ax_whole2P,jet);
    % colormap(gray);
    hold on
    if get(check_allROI,'Value') == 1
        % Outline all ROIs in blue
        for j = 1:Numcells
            plot(outline{1,1,j}{1}(:,2),outline{1,1,j}{1}(:,1),'w','Linewidth',1);
        end
    end
    % Outline ROI of current cell in green
    plot(outline{1,1,id}{1}(:,2),outline{1,1,id}{1}(:,1),'g','Linewidth',3);
    hold off
end

function showTS(id)
    % Ca time series
    plot(ax_rawCa, downt, tsG(id,:), 'Color', 'b'); axis(ax_rawCa,[0 T tsGmin tsGmax]); 
    ylabel(ax_rawCa,'F (au)','Fontweight','bold','Fontsize',12); xticks(ax_rawCa,[]);

    % delta R/R
    cellR = R(id,:);
    axes(ax_dR);
    for i = 1:length(iLoop)-1
        if mod(i,2) == 0
            color = 'r';
        else
            color = 'b';
        end
        plot(downt(iLoop(i)+1:iLoop(i+1)), cellR(iLoop(i)+1:iLoop(i+1)), 'Color',color); axis(ax_dR,[0 T Rmin Rmax]);
        hold on
    end
    hold off
    ylabel(ax_dR,'dR/R','Fontweight','bold','Fontsize',12); xticks(ax_dR,[]);

    % Spikes
    plot(ax_spikes, downt, spikes(id,:), 'Color', 'b'); axis(ax_spikes,[0 T spikesmin spikesmax]);
    ylabel(ax_spikes,'Spikes','Fontweight','bold','Fontsize',12); xlabel(ax_spikes,'Time (s)','Fontweight','bold');
end

function displaySpikeLoc(id)
    axes(ax_cellSpikeLoc);
    plot(downx,downy,'Color',[0.5 0.5 0.5],'linewidth',6); axis off; axis square;
    str = sprintf('Cell %g', id);
    title(str,'Fontweight','normal','Fontsize',11);
    hold on
    cmap = colormap(ax_cellSpikeLoc,jet(Numcells));
    ind1 = find(activespikes(id,:)>0);
    scatter(activex(ind1),activey(ind1),40,activespikes(id,ind1),'filled'); 
    c = colorbar('Location','east','Units','pixels','Position',[1120 553 16 205]);%,...
        %'Limits',[cmin cmax]);
    %c.Label.String = 'Spikes/s';
    title(c,'Spikes/s');
    hold off
end

function displayGcampMruby(id)
    bounds = outline{1,1,id}{1};
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
        if xmax > size(mean_imG,2)
            xmax = size(mean_imG,2);
        end
    ymax = max(bounds(:,1))+pad;
        if ymax > size(mean_imG,1)
            ymax = size(mean_imG,1);
        end
    xrange = xmax - xmin;
    yrange = ymax - ymin;
    if xrange < yrange
        xyrange = xrange;
    else 
        xyrange = yrange;
    end
    
    im_gcamp = imcrop(mean_imG,[xmin ymin xyrange xyrange]);
    axes(ax_gcamp);
    imagesc(im_gcamp); %colormap(ax_gcamp,gray); 
    yticks(round(xyrange/2:xyrange/2:xyrange)); yticklabels(round(pix2um*(xyrange/2:xyrange/2:xyrange)));
    ylabel('Size (um)');
    xticks([]); xlabel('GCaMP6s','Fontweight','normal','Fontsize',11);
    
    im_mruby = imcrop(mean_imR,[xmin ymin xyrange xyrange]);
    axes(ax_mruby);
    imagesc(im_mruby); %colormap(ax_mruby,gray);
    xticks([]); yticks([]); xlabel('mRuby','Fontweight','normal','Fontsize',11);
end

function displayCellPF(id)
    Nepochs = str2double(get(edit_Nepochs,'String'));
    if  Nepochs == 1
        axes(ax_ep1cellPF);
        imagesc(sorted_placeMap(id,:)); 
        xticks(Nbins/6:Nbins/6:Nbins); xticklabels(degperbin*(Nbins/6:Nbins/6:Nbins)); 
        xlabel('Position (degrees)'); yticks([]); 
        str = sprintf('Cell %g', id);
        title(str,'Fontweight','normal','Fontsize',12);
    else
        switch Nepochs
            case 2
                axes(ax_cellPF);
                imagesc(sorted_placeMap(id,:)); 
                xticks(Nbins/6:Nbins/6:Nbins); xticklabels(degperbin*(Nbins/6:Nbins/6:Nbins)); 
                xlabel('Position (degrees)'); yticks([]); 
                str = sprintf('Cell %g', id);
                ylabel('All','Position','east');
                
                axes(ax_ep1cellPF);
                imagesc(sorted_placeMap_ep(id,:,1)); 
                xticks([]); yticks([]); ylabel('1','Position','east');
                
                axes(ax_ep2cellPF);
                imagesc(sorted_placeMap_ep(id,:,2)); 
                xticks([]); yticks([]); ylabel('2','Position','east');
            case 3
                axes(ax_cellPF);
                imagesc(sorted_placeMap(id,:)); 
                xticks(Nbins/6:Nbins/6:Nbins); xticklabels(degperbin*(Nbins/6:Nbins/6:Nbins)); 
                xlabel('Position (degrees)'); yticks([]); 
                str = sprintf('Cell %g', id);
                ylabel('All','Position','east');
                
                axes(ax_ep1cellPF);
                imagesc(sorted_placeMap_ep(id,:,1)); 
                xticks([]); yticks([]); ylabel('1','Position','east');
                
                axes(ax_ep2cellPF);
                imagesc(sorted_placeMap_ep(id,:,2)); 
                xticks([]); yticks([]); ylabel('2','Position','east');
                
                axes(ax_ep3cellPF);
                imagesc(sorted_placeMap_ep(id,:,3)); 
                xticks([]); yticks([]); ylabel('3','Position','east');
            case 4
                axes(ax_cellPF);
                imagesc(sorted_placeMap(id,:)); 
                xticks(Nbins/6:Nbins/6:Nbins); xticklabels(degperbin*(Nbins/6:Nbins/6:Nbins)); 
                xlabel('Position (degrees)'); yticks([]); 
                str = sprintf('Cell %g', id);
                ylabel('All','Position','east');
                
                axes(ax_ep1cellPF);
                imagesc(sorted_placeMap_ep(id,:,1)); 
                xticks([]); yticks([]); ylabel('1','Position','east');
                
                axes(ax_ep2cellPF);
                imagesc(sorted_placeMap_ep(id,:,2)); 
                xticks([]); yticks([]); ylabel('2','Position','east');
                
                axes(ax_ep3cellPF);
                imagesc(sorted_placeMap_ep(id,:,3)); 
                xticks([]); yticks([]); ylabel('3','Position','east');
                
                axes(ax_ep4cellPF);
                imagesc(sorted_placeMap_ep(id,:,4)); 
                xticks([]); yticks([]); ylabel('4','Position','east');
        end
    end

end

function displayAllR(varargin)
    showSubset = get(check_showSubset,'Value');
    axes(ax_allR);
    if ~showSubset
        iosr.figures.multiwaveplot(downt,1:Numcells,R,'gain',8); yticks([]); xticks([]); 
        title('Smoothened dR/R','Fontweight','normal','Fontsize',12); 
    else
        Nsubset = str2double(get(edit_showSubset,'str'));
        count = round(get(slider_showSubset,'Value')); 
        bound = round( linspace(1,Numcells,Nsubset+1) );
        iosr.figures.multiwaveplot(downt,1:length(bound(count):bound(count+1)),R(bound(count):bound(count+1),:),'gain',8); 
        yticks([]); xticks([]); 
        title('Smoothened dR/R','Fontweight','normal','Fontsize',12); 
    end
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
        if get(check_showLabels,'Value') == 1
            for j = 1:Numcells
                text(centroids(j,1),centroids(j,2),num2str(j),'Color','w');
            end
        end 
    end
end

function check_showLabels_callback(varargin)
    if get(check_showLabels,'Value') == 1
        for j = 1:Numcells
            text(centroids(j,1),centroids(j,2),num2str(j),'Color','w');
        end
    else
        id = str2double(get(edit_cellNum,'string'));
        displayROI(id);    
    end
end

function check_showSubset_callback(varargin)
    if get(check_showSubset,'Value') == 0
        set(edit_showSubset,'Enable','off');
        set(slider_showSubset,'Enable','off');
    else
        set(edit_showSubset,'Enable','on');
        set(slider_showSubset,'Enable','on');
    end
    displayAllR;
end

function edit_showSubset_callback(varargin)
    Nsubset = str2double(get(edit_showSubset,'string'));

    if isnan(Nsubset)
        set(edit_showSubset,'string','4');
        displayAllR;
    else
        if ~isinteger(Nsubset)
            Nsubset = round(Nsubset);
            set(edit_showSubset,'String',num2str(Nsubset));
        end
        if Nsubset < 1 
            set(edit_cellNum,'string','2');
            displayAllR;
        elseif Nsubset > 8
            set(edit_cellNum,'string','6');
            displayAllR;
        else
            set(edit_showSubset,'String',num2str(Nsubset));
            displayAllR;
        end
    end
end

function edit_cellNum_callback(varargin)
    id = str2double(get(edit_cellNum,'string'));

    if isnan(id)
        set(edit_cellNum,'string','1');
        set(slider_cellNum,'Value',1);
        displayROI(1);
        showTS(1);
        displaySpikeLoc(1);
        displayGcampMruby(id);
        displayCellPF(1);
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
            displaySpikeLoc(1);
            displayGcampMruby(1);
            displayCellPF(1);
        elseif id > Numcells
            set(edit_cellNum,'string',num2str(Numcells));
            set(slider_cellNum,'Value',Numcells);
            displayROI(Numcells);
            showTS(Numcells);
            displaySpikeLoc(Numcells);
            displayGcampMruby(Numcells);
            displayCellPF(Numcells);
        else
            set(edit_cellNum,'String',num2str(id));
            displayROI(id);
            showTS(id);
            displaySpikeLoc(id);
            displayGcampMruby(id);
            displayCellPF(id);
            set(slider_cellNum,'Value',id);
        end
    end
end

function slider_cellNum_callback(varargin)
    id = round(get(slider_cellNum,'Value'));
    displayROI(id);
    showTS(id);
    displaySpikeLoc(id);
    displayGcampMruby(id);
    displayCellPF(id);
    set(edit_cellNum,'string',num2str(id));
end

function edit_Nepochs_callback(varargin)
    Nepochs = str2double(get(edit_Nepochs,'string'));

    if isnan(Nepochs)
        set(edit_Nepochs,'string','2');
    else
        if ~isinteger(Nepochs)
            Nepochs = round(Nepochs);
            set(edit_Nepochs,'String',num2str(Nepochs));
        end
        if Nsubset < 1 
            set(edit_Nepochs,'string','2');
        elseif Nsubset > 4
            set(edit_Nepochs,'string','6');
        end
    end
    
    if Nepochs ~= 1
        set(hdl_gui,'Position',[1087 648 1770 1040]);
        set(panel_epochs,'Visible','on');
        switch Nepochs
            case 2 
                axes(ax_epoch1); imagesc(sorted_)
            case 3
            case 4
        end
    else
        set(hdl_gui,'Position',[1087 648 1200 1040]);
        set(panel_epochs,'Visible','off')
    end
end

end