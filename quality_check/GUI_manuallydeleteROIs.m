% Written by Ann Go
% GUI for refining ROIs
% GUI shows ROI mask and timeseries extracted by CaImAn. User may scroll
% through ROIs with slider arrows or choose ROIs using mouse pointer.
% Before saving refined ROIs, GUI makes a copy of unrefined ROIs and the
% figure summaries in a folder titled 'unrefined_ROIs/'.

function GUI_manuallydeleteROIs(list, reffile)

%% GUI structures
% GUI figure
hfig = figure('Name','Manually delete ROIs','NumberTitle','off','Resize','off',...
    'Position',[1268 1044 880 820]);

% 2P image with ROIs
ax_image = axes('Units','pixels','Position',[29 231 550 550],'Box','off','Visible','on','clipping','off'); 
[ mouseid, expname, fov ] = find_mouseIDexpname(list);
im_str = strrep( [mouseid ': ' ' ' expname], '_', ' ' );
title(im_str,'Fontweight','bold','Fontsize',14);

% Time series
ax_tsG = axes('Units','pixels','Position',[29 131 833 50],'Box','off','Visible','on','clipping','off');
title('Ca time series','Fontweight','bold','Fontsize',14);
ax_df = axes('Units','pixels','Position',[29 27 834 50],'Box','off','Visible','on','clipping','off');
title('Delta f/f','Fontweight','bold','Fontsize',14);

% control buttons
% pbutton_usemouse 
uicontrol('style','pushbutton','Units','pixels','Position',[610 600 250 37],'string','Use mouse to choose ROI',...
    'Fontsize',14,'Callback',@pbutton_usemouse_callback);
% pbutton_delete 
uicontrol('style','pushbutton','Units','pixels','Position',[655 550 150 37],'string','Delete ROI',...
    'Fontsize',14,'Callback',@pbutton_delete_callback);
% pbutton_undodelete 
uicontrol('style','pushbutton','Units','pixels','Position',[655 500 150 37],'string','Undo delete ROI',...
    'Fontsize',14,'Callback',@pbutton_undodelete_callback);
% pbutton_reset 
uicontrol('style','pushbutton','Units','pixels','Position',[655 450 150 37],'string','RESET ROIs',...
    'Fontsize',14,'Callback',@pbutton_reset_callback);
% pbutton_exit 
uicontrol('style','pushbutton','Units','pixels','Position',[630 260 200 37],'string','EXIT',...
    'Fontsize',14,'Callback',@pbutton_exit_callback);

check_finalise = uicontrol('style','checkbox','Units','pixels','Position',[610 350 100 26],'string','Finalise',...
    'Fontsize',14,'Value',0,'Callback',@check_finalise_callback);
pbutton_save = uicontrol('style','pushbutton','Units','pixels','Position',[717 343 118 37],'string','SAVE',...
    'Fontsize',14,'Enable','off','Callback',@pbutton_save_callback);


%% load data
[data_locn,~,err] = load_neuroSEEmodules;
if ~isempty(err)
    beep
    cprintf('Errors',err);    
    return
end
if ~isempty(fov)
    grp_sdir = [data_locn 'Analysis/' mouseid '/' fov '/' mouseid '_' expname '/group_proc/imreg_normcorre_CaImAn/'...
            mouseid '_' expname '_imreg_ref' reffile '/'];
else
    grp_sdir = [data_locn 'Analysis/' mouseid '/' mouseid '_' expname '/group_proc/imreg_normcorre_CaImAn/'...
            mouseid '_' expname '_imreg_ref' reffile '/'];
end

% data
fname_roidata = [grp_sdir mouseid '_' expname '_ref' reffile '_segment_output.mat'];        
roidata = load(fname_roidata);
im = roidata.corr_image;
masks = roidata.masks; 
inc = 1:size(masks,3);
exc = [];
tsG = roidata.tsG(:,1:min([30000, size(roidata.tsG,2)]));   % only read the first 30000 points as very long arrays 
                                                    % take too long to load and refresh gui
df_f = roidata.df_f(:,1:min([30000, size(roidata.tsG,2)]));
t = (1:length(tsG))/30;

Nrois = size(masks,3);
outline = getOutlines(masks);

% set app data values
inc_app = inc;    % app-based indices, manipulated in the app
exc_app = exc;    % app-based indices, manipulated in the app

% Cell number
str = sprintf('ROI number (1 to %g):', Nrois);
text_roiNum = uicontrol('style','text','Units','pixels','string',str,'Fontsize',12,...
    'FontWeight','bold','Position',[600 720 171 23]); 
edit_roiNum = uicontrol('style','edit','Units','pixels','string','1','Fontsize',14,...
    'Position',[764 720 45 26],'Callback',@edit_roiNum_callback); 
slider_roiNum = uicontrol('style','slider','Units','pixels','Position',[665 685 101 22],...
    'Min',1,'Max',Nrois,'Value',1,'SliderStep',[1/Nrois 1/Nrois],'Callback',@slider_roiNum_callback); 

% Initialise GUI
displayROI(1);
showTS(1);

%% Subfunctions and Callback functions
function displayROI(id)
    % Display neuron image overlaid with specified ROI mask
    axes(ax_image);
    imagesc(im); 
    im_str = strrep( [mouseid ': ' ' ' expname], '_', ' ' );
    title(im_str,'Fontweight','bold','Fontsize',14);
    colormap(ax_image,gray); axis off; 
    hold on
    % Outline all ROIs 
    for i = 1:length(inc_app)
        k = inc_app(i);
        plot(outline{1,1,k}{1}(:,2),outline{1,1,k}{1}(:,1),'w','Linewidth',1);
    end

    % Outline ROI of current cell in green
    if ismember(id,exc_app)
        plot(outline{1,1,id}{1}(:,2),outline{1,1,id}{1}(:,1),'r','Linewidth',3);
    else
        plot(outline{1,1,id}{1}(:,2),outline{1,1,id}{1}(:,1),'g','Linewidth',3);
    end
    hold off
end

function showTS(id)
    % Ca time series
    axes(ax_tsG);
    plot(t, tsG(id,:), 'Color', 'b');  
    title('Ca time series','Fontweight','bold','Fontsize',14);

    % df/f
    axes(ax_df);
    plot(t, df_f(id,:), 'Color', 'r');  
    title('Delta f/f','Fontweight','bold','Fontsize',14);   
end

function outline = getOutlines(masks)
    for j = 1:size(masks,3)
        outline{:,:,j} = bwboundaries(masks(:,:,j));    % boundary of each ROI
    end
end

function edit_roiNum_callback(varargin)
    id = str2double(get(edit_roiNum,'string'));
    if isnan(id)
        set(edit_roiNum,'string','1');
        set(slider_roiNum,'Value',1);
        displayROI(1);
        showTS(1);
    else
        if ~isinteger(id)
            id = round(id);
            set(edit_roiNum,'String',num2str(id));
        end
        if id < 1 
            set(edit_roiNum,'string','1');
            set(slider_roiNum,'Value',1);
            displayROI(1);
            showTS(1);
        elseif id > Nrois
            set(edit_roiNum,'string',num2str(Nrois));
            set(slider_roiNum,'Value',Nrois);
            displayROI(Nrois);
            showTS(Nrois);
        else
            set(edit_roiNum,'String',num2str(id));
            set(slider_roiNum,'Value',id);
            displayROI(id);
            showTS(id);
        end
    end
end

function slider_roiNum_callback(varargin)
    id = round(get(slider_roiNum,'Value'));
    displayROI(id);
    showTS(id);
    set(edit_roiNum,'string',num2str(id));
end

function pbutton_usemouse_callback(varargin)
    axes(ax_image);
    [x,y] = ginput(1);
    
    for j = 1:size(masks,3)
        ind(j) = masks(round(y),round(x),j);
    end
    id = find(ind);
    if length(id) > 1
        id = id(1);
    end
    if length(id) == 1
        displayROI(id);
        showTS(id);
        set(edit_roiNum,'String',num2str(id));
        edit_roiNum_callback;
    end
end

function pbutton_delete_callback(varargin)
    id = str2double(get(edit_roiNum,'string'));
    inc_app = inc_app(inc_app~=id);
    exc_app = sort([exc_app, id]);
    id = id + 1;
    if id > Nrois
        id = 1;
    end
    
    set(edit_roiNum,'String',num2str(id));
    edit_roiNum_callback;
end

function pbutton_undodelete_callback(varargin)
    id = str2double(get(edit_roiNum,'string'));
    inc_app = unique(sort([inc_app, id]));
    exc_app = exc_app(exc_app~=id);
    
    displayROI(id);
    showTS(id);
end
function pbutton_reset_callback(varargin)
    inc_app = inc;              
    exc_app = exc;
    
    id = str2double(get(edit_roiNum,'string'));
    displayROI(id);
    showTS(id); 
end

function check_finalise_callback(varargin)
    checkval = get(check_finalise,'Value');
    if checkval
        set(pbutton_save,'Enable','On')
    else
        set(pbutton_save,'Enable','Off')
    end
end

function pbutton_save_callback(varargin)
    set(check_finalise,'Value',0);
    set(pbutton_save,'Enable','Off');
    
    % First make a copy of the ROI segmentation output and ROI and elim_ROI
    % figures
    copydir = [grp_sdir 'unrefined_ROIs/'];
    if ~exist( copydir, 'dir' ), mkdir( copydir ); fileattrib(copydir,'+w','g','s'); end
    copyfile(fname_roidata,...
                 [copydir mouseid '_' expname '_ref' reffile '_segment_output.mat'])
    copyfile([grp_sdir mouseid '_' expname '_ref' reffile '_elimROIs.fig'],...
             [copydir mouseid '_' expname '_ref' reffile '_elimROIs.fig'])
    copyfile([grp_sdir mouseid '_' expname '_ref' reffile '_elimROIs.png'],...
             [copydir mouseid '_' expname '_ref' reffile '_elimROIs.png'])
    copyfile([grp_sdir mouseid '_' expname '_ref' reffile '_ROIs.fig'],...
             [copydir mouseid '_' expname '_ref' reffile '_ROIs.fig'])
    copyfile([grp_sdir mouseid '_' expname '_ref' reffile '_ROIs.png'],...
             [copydir mouseid '_' expname '_ref' reffile '_ROIs.png'])
        
    % Now save refined ROI segmentation output
    masks_all = masks;
    masks_all(:,:,size(masks,3)+1:size(masks,3)+size(roidata.elim_masks,3)) = roidata.elim_masks;
    tsG_all(1:size(masks,3),:) = roidata.tsG;
    tsG_all(size(masks,3)+1:size(masks,3)+size(roidata.elim_masks,3),:) = roidata.elim_tsG;
    df_f_all(1:size(masks,3),:) = roidata.df_f;
    df_f_all(size(masks,3)+1:size(masks,3)+size(roidata.elim_masks,3),:) = roidata.elim_df_f;
    
    inc = inc_app;
    exc = [exc_app, size(masks,3)+1:size(masks,3)+size(roidata.elim_masks,3)];
    masks = masks_all(:,:,inc);     

    output.tsG = tsG_all(inc,:);        output.elim_tsG = tsG_all(exc,:);
    output.df_f = df_f_all(inc,:);      output.elim_df_f = df_f_all(exc,:);
    output.masks = masks;               output.elim_masks = masks_all(:,:,exc);
    output.corr_image = im;
    output.F0 = roidata.F0;
    output.A = roidata.A;
    output.params = roidata.params;       
    save(fname_roidata,'-struct','output');

    % Regenerate plots for masks and excluded masks
    plotopts.plot_ids = 1; % set to 1 to view the ID number of the ROIs on the plot
    fig = plotContoursOnSummaryImage(im, masks, plotopts);
    figname_pref = [grp_sdir mouseid '_' expname '_ref' reffile];
    savefig(fig, [figname_pref '_ROIs']);
    saveas(fig, [figname_pref '_ROIs'], 'png');
    close(fig);
    
    fig = plotContoursOnSummaryImage(im, output.elim_masks, plotopts);
    figname_pref = [grp_sdir mouseid '_' expname '_ref' reffile];
    savefig(fig, [figname_pref '_elimROIs']);
    saveas(fig, [figname_pref '_elimROIs'], 'png');
    close(fig);
    
    % Refresh gui with saved data
    inc = 1:size(masks,3);  inc_app = inc;              
    exc = [];               exc_app = exc;
    tsG = tsG(inc,:);
    df_f = df_f(inc,:);
    
    Nrois = size(masks,3);
    outline = getOutlines(masks);
    str = sprintf('ROI number (1 to %g):', Nrois);
    set(text_roiNum,'String',str);
    set(edit_roiNum,'String','1');
    set(slider_roiNum,'Value',1);
    displayROI(1);
    showTS(1);
end

function pbutton_exit_callback(varargin)
    close(hfig);
end
end