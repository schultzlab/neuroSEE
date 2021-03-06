% Written by Ann Go
% 
% This GUI allows user to manually refine spike output of
% neuroSEE_spkExtract.
%
% INPUTS:
%   spikes      : spike estimates
%   tsG         : raw timeseries (optional if dtsG exists)
%   dtsG        : fissa-corrected timeseries (optional if tsG exists) 
%   df_f        : dF/F 
%   ddf_f       : fissa-corrected dF/F (optional) 
%   data_locn   : repository for GCaMP data
%   file        : part of filename of 2P image in the format
%                  yyyymmdd_HH_MM_SS
%   params      : complete params 
%   corr_image  : correlation image from green channel (optional)
%   masks       : ROI masks (optional)
%
% Note for optional inputs: use empty matrix
%   e.g. [spikes, params] = GUI_manually_refine_spikes( spikes, tsG, [], df_f, [], data_locn, file, params, [], [])


function GUI_manually_refine_spikes_imreg( origspikes, tsG, dtsG, df_f, ddf_f, corr_image, masks, origparams, globalvar, savedir, savename_pref )

if nargin < 7, masks = []; end 
if nargin < 8, corr_image = zeros(512,512); end
if nargin < 9, globalvar = false; end
if nargin < 10, fsave = false; else, fsave = true; end
if ~isempty(tsG)
    N = size(tsG,1); T = size(tsG,2);
elseif ~isempty(dtsG)
    N = size(dtsG,1); T = size(dtsG,2);
elseif ~isempty(df_f)
    N = size(df_f,1); T = size(df_f,2);
end 

if globalvar, global spikes params; end    


%% GUI structures

% GUI figure
hfig_h = 581;
hfig_w = 1756;
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
      
% baseline percentage
panel_prctile = uipanel('Parent',fig_gui,'Units','pixels','Position',[562 30 161 101]);
text_prctile = uicontrol('Parent',panel_prctile,'style','text','Units','pixels','string',...
    'Baseline percentile','Fontsize',12); 
    set(text_prctile,'Position',[14 61 128 29]);
edit_prctile = uicontrol('Parent',panel_prctile,'style','edit','Units','pixels',...
    'string',num2str(origparams.spkExtract.bl_prctile),'Enable','on','Callback',@edit_prctile_callback); 
    set(edit_prctile,'Position',[10 40 40 22]);
slider_prctile = uicontrol('Parent',panel_prctile,'style','slider','Units','pixels',...
    'Min',0,'Max',100,'Value',origparams.spkExtract.bl_prctile,'SliderStep',[0.05 0.05],'Callback',@slider_prctile_callback); 
    set(slider_prctile,'Position',[60 35 88 24]);
rbutton_autoprctile = uicontrol('Parent',panel_prctile,'style','radiobutton','Units','pixels',...
    'string','auto','Value',0,'Enable','off','Fontsize',12,'Callback',@rbutton_autoprctile_callback); 
    set(rbutton_autoprctile,'Position',[13 9 59 26]);
edit_autoprctile = uicontrol('Parent',panel_prctile,'style','edit','Units','pixels',...
    'string','0','Enable','off','Callback',@edit_autoprctile_callback); 
    set(edit_autoprctile,'Position',[70 10 40 22]);

% options
panel_option = uipanel('Parent',fig_gui,'Units','pixels','Position',[745 30 228 101]);
text_spikeSNR = uicontrol('Parent',panel_option,'style','text','Units','pixels','string',...
    'spike SNR','Fontsize',12); 
    set(text_spikeSNR,'Position',[17 70 73 16]);
edit_spikeSNR = uicontrol('Parent',panel_option,'style','edit','Units','pixels',...
    'string',num2str(origparams.spkExtract.spk_SNR),'Enable','on','Callback',@edit_spkSNR_callback); 
    set(edit_spikeSNR,'Position',[137 68 58 19]);
text_decaytime = uicontrol('Parent',panel_option,'style','text','Units','pixels','string',...
    'decay time','Fontsize',12); 
    set(text_decaytime,'Position',[17 47 73 16]);
edit_decaytime = uicontrol('Parent',panel_option,'style','edit','Units','pixels',...
    'string',num2str(origparams.spkExtract.decay_time),'Enable','on','Callback',@edit_decaytime_callback); 
    set(edit_decaytime,'Position',[137 45 58 19]);
text_lamprob = uicontrol('Parent',panel_option,'style','text','Units','pixels','string',...
    'lambda probability','Fontsize',12); 
    set(text_lamprob,'Position',[17 24 115 16]);
edit_lamprob = uicontrol('Parent',panel_option,'style','edit','Units','pixels',...
    'string',num2str(origparams.spkExtract.lam_pr),'Enable','on','Callback',@edit_lamprob_callback); 
    set(edit_lamprob,'Position',[137 22 58 19]);

% Ca time series plot
ax_F = axes('Parent',fig_gui,'Units','pixels','Position',[534 445 1200 95],'Xlim',[0 1],'Ylim',[0  1],'Box','off',...
            'Visible','on','clipping','off'); 

% df/f0 
ax_dF = axes('Parent',fig_gui,'Units','pixels','Position',[534 313 1200 95],'Xlim',[0 1],'Ylim',[0  1],'Box','off',...
            'Visible','on','clipping','off');

% spikes plot
ax_spikes = axes('Parent',fig_gui,'Units','pixels','Position',[534 180 1200 95],'Xlim',[0 1],'Ylim',[0  1],'Box','off',...
            'Visible','on','clipping','off'); 

% reset, apply to all, save, quit
pbutton_reset = uicontrol('Parent',fig_gui,'style','pushbutton','Units','pixels','string','Reset','Fontsize',12,...
    'Position',[1000 88 125 33],'Enable','off','Value',0,'Callback',@pbutton_reset_callback);

tbutton_applytoall = uicontrol('Parent',fig_gui,'style','togglebutton','Units','pixels','string','Apply to all','Fontsize',12,...
    'Position',[1000 41 125 33],'Enable','off','Value',0,'Callback',@tbutton_applytoall_callback);

pbutton_save = uicontrol('Parent',fig_gui,'style','pushbutton','Units','pixels','string','SAVE','Fontsize',12,...
    'Position',[1175 61 125 44],'Enable','off','Value',0,'Callback',@pbutton_save_callback);

uicontrol('Parent',fig_gui,'style','pushbutton','Units','pixels','string','QUIT','Fontsize',12,...
    'Position',[1500 61 125 44],'Enable','on','Value',0,'Callback',@pbutton_quit_callback);

    
%% initial GUI data
if ~isempty(dtsG)
    C1 = dtsG; str_C1 = 'FISSA-corrected timeseries';
else
    C1 = tsG; str_C1 = 'Raw timeseries';
end
if ~isempty(ddf_f)
    C2 = ddf_f; str_C2 = 'FISSA-corrected dF/F';
else
    C2 = df_f; str_C2 = 'dF/F';
end

df_prctile(1:N) = origparams.spkExtract.bl_prctile;
lam_pr(1:N) = origparams.spkExtract.lam_pr;
spk_SNR(1:N) = origparams.spkExtract.spk_SNR;

setappdata(fig_gui,'curr_C2',C2);
setappdata(fig_gui,'curr_spikes',origspikes);
setappdata(fig_gui,'saveCount',0);
setappdata(fig_gui,'df_prctile',df_prctile);
setappdata(fig_gui,'lam_pr',lam_pr);
setappdata(fig_gui,'spk_SNR',spk_SNR);

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
                    plot(trace(:,2),trace(:,1),'w','Linewidth',1);
                end
            end
        end
        outline = bwboundaries(masks(:,:,id));
        if size(outline,1) > 0
            trace = outline{1};
            plot(trace(:,2),trace(:,1),'g','Linewidth',2);
        end
        hold off
    end
end

function edit_cellNum_callback(varargin)
    id = str2double(get(edit_cellNum,'string'));

    if isnan(id)
        id = 1;
    else
        if id < 1 
            id = 1;
        elseif id > N
            id = N;
        end
    end
    set(edit_cellNum,'string',num2str(id));
    set(slider_cellNum,'Value',id);
    displayROI(id); showTS(id);
    
    df_prctile = getappdata(fig_gui,'df_prctile');
    lam_pr = getappdata(fig_gui,'lam_pr');
    set(edit_prctile,'string',num2str(df_prctile(id)));
    set(slider_prctile,'value',df_prctile(id));
    set(edit_lamprob,'string',num2str(lam_pr(id)));
end

function slider_cellNum_callback(varargin)
    id = round(get(slider_cellNum,'Value'));
    set(edit_cellNum,'string',num2str(id));
    displayROI(id); showTS(id);
    
    df_prctile = getappdata(fig_gui,'df_prctile');
    lam_pr = getappdata(fig_gui,'lam_pr');
    set(edit_prctile,'string',num2str(df_prctile(id)));
    set(slider_prctile,'value',df_prctile(id));
    set(edit_lamprob,'string',num2str(lam_pr(id)));
end

function check_allROI_callback(varargin)
    id = str2double(get(edit_cellNum,'string'));
    displayROI(id);    
end

function edit_prctile_callback(varargin)
    id = str2double(get(edit_cellNum,'string'));
    val = str2double(get(edit_prctile,'string'));
    set(pbutton_reset,'Enable','on'); 
    set(tbutton_applytoall,'Enable','on'); 
    if fsave, set(pbutton_save,'Enable','on'); end
    
    df_prctile = getappdata(fig_gui,'df_prctile');
    if isnan(val)
        set(slider_prctile,'Value',df_prctile(id));
    else
        if val < 0 
            set(edit_prctile,'string','0');
            set(slider_prctile,'Value',0);
        elseif val > 100
            set(edit_prctile,'string','100');
            set(slider_prctile,'Value',100);
            set(slider_prctile,'Value',df_prctile(id));
        end
    end
    if get(tbutton_applytoall,'Value') 
        for ii = 1:N
            df_prctile(ii) = val;
        end
        calcC2_all;
    else
        df_prctile(id) = val;
        calcC2(id);
    end
    setappdata(fig_gui,'df_prctile',df_prctile);
end

function slider_prctile_callback(varargin)
    id = str2double(get(edit_cellNum,'string'));
    val = round(get(slider_prctile,'Value'));
    set(edit_prctile,'string',num2str(val));
    set(pbutton_reset,'Enable','on'); 
    set(tbutton_applytoall,'Enable','on'); 
    df_prctile = getappdata(fig_gui,'df_prctile');
    if fsave, set(pbutton_save,'Enable','on'); end
    if get(tbutton_applytoall,'Value') 
        for ii = 1:N
            df_prctile(ii) = val;
        end
        calcC2_all;
    else
        df_prctile(id) = val;
        calcC2(id);
    end
    setappdata(fig_gui,'df_prctile',df_prctile);
end

function rbutton_autoprctile_callback(varargin)
    set(pbutton_reset,'Enable','on'); 
    set(tbutton_applytoall,'Enable','on'); 
    if fsave, set(pbutton_save,'Enable','on'); end
    if get(rbutton_autoprctile,'Value')
        set(edit_prctile,'Enable','off');
        set(slider_prctile,'Enable','off');
    else
        set(edit_prctile,'Enable','on');
        set(slider_prctile,'Enable','on');
    end
    if get(tbutton_applytoall,'Value') 
        calcC2_all;
    else
        calcC2(id);
    end
end
        
function edit_spkSNR_callback(varargin)
    id = str2double(get(edit_cellNum,'string'));
    val = str2double(get(edit_spikeSNR,'string'));
    spk_SNR = getappdata(fig_gui,'spk_SNR');
    set(pbutton_reset,'Enable','on'); 
    set(tbutton_applytoall,'Enable','on'); 
    if fsave, set(pbutton_save,'Enable','on'); end
    if get(tbutton_applytoall,'Value') 
        for ii = 1:N
            spk_SNR(ii) = val;
        end
        calcSpikes_all;
    else
        spk_SNR(id) = val;
        calcSpikes(id);
    end
    setappdata(fig_gui,'spk_SNR',spk_SNR);
end

function edit_decaytime_callback(varargin)
    set(pbutton_reset,'Enable','on'); 
    if fsave, set(pbutton_save,'Enable','on'); end
    if get(tbutton_applytoall,'Value') 
        calcSpikes_all;
    else
        id = str2double(get(edit_cellNum,'string'));
        calcSpikes(id);
    end
end

function edit_lamprob_callback(varargin)
    id = str2double(get(edit_cellNum,'string'));
    val = str2double(get(edit_lamprob,'string'));
    lam_pr = getappdata(fig_gui,'lam_pr');
    set(pbutton_reset,'Enable','on'); 
    set(tbutton_applytoall,'Enable','on'); 
    if fsave, set(pbutton_save,'Enable','on'); end
    if get(tbutton_applytoall,'Value') 
        for ii = 1:N
            lam_pr(ii) = val;
        end
        calcSpikes_all;
    else
        lam_pr(id) = val;
        calcSpikes(id);
    end
    setappdata(fig_gui,'lam_pr',lam_pr);
end

function pbutton_reset_callback(varargin)
    set(pbutton_reset,'Enable','off');
    if getappdata(fig_gui,'saveCount')
        if fsave, set(pbutton_save,'Enable','on'); end
    else
        set(pbutton_save,'Enable','off');
    end
    set(edit_prctile,'Enable','on');
        set(edit_prctile,'string',num2str(origparams.spkExtract.bl_prctile));
        setappdata(fig_gui,'df_prctile',origparams.spkExtract.bl_prctile);
    set(slider_prctile,'Enable','on');
        set(slider_prctile,'Value',origparams.spkExtract.bl_prctile);
    set(rbutton_autoprctile,'Value',0);
    set(edit_spikeSNR,'String',num2str(origparams.spkExtract.spk_SNR));
    set(edit_decaytime,'String',num2str(origparams.spkExtract.decay_time));
    set(edit_lamprob,'String',num2str(origparams.spkExtract.lam_pr));
        setappdata(fig_gui,'lam_pr',origparams.spkExtract.lam_pr);
    
    setappdata(fig_gui,'curr_C2',C2);
    setappdata(fig_gui,'curr_spikes',origspikes);
    id = str2double(get(edit_cellNum,'string'));
    showTS(id);
end

function tbutton_applytoall_callback(varargin)
    if get(tbutton_applytoall,'Value')
        df_prctile = getappdata(fig_gui,'df_prctile');
        spk_SNR = getappdata(fig_gui,'spk_SNR');
        lam_pr = getappdata(fig_gui,'lam_pr');
        for ii = 1:N
            df_prctile(ii) = str2double(get(edit_prctile,'string'));
            spk_SNR(ii) = str2double(get(edit_spikeSNR,'string'));
            lam_pr(ii) = str2double(get(edit_lamprob,'string'));
        end
        calcC2_all;
        
        setappdata(fig_gui,'df_prctile',df_prctile);
        setappdata(fig_gui,'spk_SNR',spk_SNR);
        setappdata(fig_gui,'lam_pr',lam_pr);
    end
end

function pbutton_save_callback(varargin)
    setappdata(fig_gui,'saveCount',1);
    spikes = getappdata(fig_gui,'curr_spikes');
    
    if globalvar
        if ~get(rbutton_autoprctile,'Value')
            outparams.bl_prctile = str2double(get(edit_prctile,'string'));
        else
            outparams.bl_prctile = getappdata(fig_gui,'auto_prctile');
        end
        outparams.spk_SNR = str2double(get(edit_spikeSNR,'string'));
        outparams.decay_time = str2double(get(edit_decaytime,'string'));
        outparams.lam_pr = str2double(get(edit_lamprob,'string'));
        params = origparams;
        params.spkExtract = outparams;
        
        setappdata(fig_gui,'saved_params',output.params);
    end

    fname_mat = [savedir savename_pref '_spikes_mrefined.mat'];

    save(fname_mat,'spikes');
    setappdata(fig_gui,'saved_spikes',spikes);
    set(pbutton_save,'Enable','off');
end

function pbutton_quit_callback(varargin)
    close(fig_gui);
end

function calcC2(id)
    df_prctile_id = str2double(get(edit_prctile,'string'));
    fo = ones(1,T) * prctile(C2(id,:),df_prctile_id);
    C2_id= (C2(id,:) - fo); % ./ fo;

    curr_C2 = getappdata(fig_gui,'curr_C2');
    curr_C2(id,:) = C2_id;
    setappdata(fig_gui,'curr_C2',curr_C2);
    calcSpikes(id);
end

function calcSpikes(id)
    C = getappdata(fig_gui,'curr_C2');
    curr_spikes = getappdata(fig_gui,'curr_spikes');
    decay_time = str2double(get(edit_decaytime,'String'));
    lam_prob = str2double(get(edit_lamprob,'String'));
    spike_SNR = str2double(get(edit_spikeSNR,'String'));
    
    spkmin = spike_SNR*GetSn(C(id,:));
    lam = choose_lambda(exp(-1/(30.9*decay_time)),GetSn(C(id,:)),lam_prob);

    [~,spkAR,~] = deconvolveCa(C(id,:),'ar2','method','thresholded','optimize_pars',true,'maxIter',20,...
                                'window',150,'lambda',lam,'smin',spkmin);
    curr_spikes(id,:) = spkAR(:);
    
    setappdata(fig_gui,'curr_spikes',curr_spikes);
    showTS(id);
end

function calcC2_all(varargin)
    if get(rbutton_autoprctile,'Value')
        auto_prctile = zeros(1,N);
    end
    curr_C2 = zeros(size(C2));
    for i = 1:N
        if get(rbutton_autoprctile,'Value')
            [~,cdf_val] = estimate_percentile_level(C2(i,:));
            df_prctile = mean(cdf_val);
            auto_prctile(i) = df_prctile;
        else
            df_prctile = str2double(get(edit_prctile,'string'));
        end
        fo = ones(1,T) * prctile(C2(i,:),df_prctile);
        curr_C2(i,:) = (C2(i,:) - fo); % ./ fo;
    end
    if get(rbutton_autoprctile,'Value')
        setappdata(fig_gui,'auto_prctile',auto_prctile);
    end
    setappdata(fig_gui,'curr_C2',curr_C2);
    calcSpikes_all;
end

function calcSpikes_all(varargin)
    C = getappdata(fig_gui,'curr_C2');
    decay_time = str2double(get(edit_decaytime,'String'));
    lam_pr = str2double(get(edit_lamprob,'String'));
    spk_SNR = str2double(get(edit_spikeSNR,'String'));
    
    curr_spikes = zeros(N,T);        
    for i = 1:N
        spkmin = spk_SNR*GetSn(C(i,:));
        lam = choose_lambda(exp(-1/(30.9*decay_time)),GetSn(C(i,:)),lam_pr);
        
        [~,spkAR,~] = deconvolveCa(C(i,:),'ar2','method','thresholded','optimize_pars',true,'maxIter',20,...
                                    'window',150,'lambda',lam,'smin',spkmin);
        curr_spikes(i,:) = spkAR(:);

        h = waitbar(i/N);
    end
    close(h);
    setappdata(fig_gui,'curr_spikes',curr_spikes);
    id = str2double(get(edit_cellNum,'string'));
    showTS(id);
end

function showTS(id)
    % Ca time series
    axes(ax_F);
    plot(C1(id,:)); title(str_C1);

    % df/f or dR/R
    curr_C2 = getappdata(fig_gui,'curr_C2');
    axes(ax_dF);
    plot(curr_C2(id,:)); title(str_C2);

    % spikes
    curr_spikes = getappdata(fig_gui,'curr_spikes');
    axes(ax_spikes);
    plot(curr_spikes(id,:)); title('Spikes');
    
    df_prctile = getappdata(fig_gui,'auto_prctile');
    if get(rbutton_autoprctile,'Value')
        set(edit_autoprctile,'String',num2str(df_prctile(id)));
    end
end

end