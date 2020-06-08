% Written by Ann Go
% 
% This GUI plots
% (1) the decontaminated time series OR if it does not exist, the raw time
%       series
% (2) the decontaminated df_f OR if it does not exist, df_f, OR in the case
%       of older files, it plots dR_R
% (3) the extracted spikes from using both oasisAR1 and oasisAR2 for
% comparison
% It also lets user see how spike output changes with different settings.
% 
% Note for missing inputs: use empty matrix
%   e.g. GUI_extractSpikes_oasisAR1vsAR2( tsG, df_f, [], [], [] )
% Required inputs:
% (1) either tsG or dtsG
% (2) one of the following: df_f, ddf_f, R


function GUI_extractSpikes_oasisAR1vsAR2( tsG, df_f, dtsG, ddf_f, R )
N = size(tsG,1); T = size(tsG,2);

%% GUI structures

% GUI figure
hfig_h = 674;
hfig_w = 1001;
hdl_gui = figure('Name','oasisAR1 vs AR2','NumberTitle','off',...%'Resize','off',...
'Position',[1680-hfig_w,1050-hfig_h,hfig_w,hfig_h]);

% cell number 
panel_cellNum = uipanel('Parent',hdl_gui,'Units','pixels','Position',[50 23 187 101]);
str = sprintf('Cell number: 1 to %g', N);
text_cellNum = uicontrol('Parent',panel_cellNum,'style','text','Units','pixels','string',str,'Fontsize',14); 
    set(text_cellNum,'Position',[4 56 169 28]);
edit_cellNum = uicontrol('Parent',panel_cellNum,'style','edit','Units','pixels','string','1','Callback',@edit_cellNum_callback); 
    set(edit_cellNum,'Position',[10 21 40 28]);
slider_cellNum = uicontrol('Parent',panel_cellNum,'style','slider','Units','pixels',...
    'Min',1,'Max',N,'Value',1,'SliderStep',[1/N 1/N],'Callback',@slider_cellNum_callback); 
    set(slider_cellNum,'Position',[64 18 108 24]);

% baseline percentage
panel_prctile = uipanel('Parent',hdl_gui,'Units','pixels','Position',[259 23 161 101]);
text_prctile = uicontrol('Parent',panel_prctile,'style','text','Units','pixels','string',...
    'Baseline percentile','Fontsize',12); 
    set(text_prctile,'Position',[14 61 128 29]);
edit_prctile = uicontrol('Parent',panel_prctile,'style','edit','Units','pixels',...
    'string','0','Enable','on','Callback',@edit_prctile_callback); 
    set(edit_prctile,'Position',[10 40 40 22]);
slider_prctile = uicontrol('Parent',panel_prctile,'style','slider','Units','pixels',...
    'Min',0,'Max',100,'Value',0,'SliderStep',[0.05 0.05],'Callback',@slider_prctile_callback); 
    set(slider_prctile,'Position',[60 35 88 24]);
rbutton_autoprctile = uicontrol('Parent',panel_prctile,'style','radiobutton','Units','pixels',...
    'string','auto','Value',0,'Enable','on','Fontsize',12,'Callback',@rbutton_autoprctile_callback); 
    set(rbutton_autoprctile,'Position',[13 9 59 26]);
edit_autoprctile = uicontrol('Parent',panel_prctile,'style','edit','Units','pixels',...
    'string','0','Enable','off','Callback',@edit_autoprctile_callback); 
    set(edit_autoprctile,'Position',[70 10 40 22]);

% options
panel_option = uipanel('Parent',hdl_gui,'Units','pixels','Position',[450 23 228 101]);
text_spikeSNR = uicontrol('Parent',panel_option,'style','text','Units','pixels','string',...
    'spike SNR','Fontsize',12); 
    set(text_spikeSNR,'Position',[17 70 73 16]);
edit_spikeSNR = uicontrol('Parent',panel_option,'style','edit','Units','pixels',...
    'string','0.5','Enable','on','Callback',@edit_options_callback); 
    set(edit_spikeSNR,'Position',[137 68 58 19]);
text_decaytime = uicontrol('Parent',panel_option,'style','text','Units','pixels','string',...
    'decay time','Fontsize',12); 
    set(text_decaytime,'Position',[17 47 73 16]);
edit_decaytime = uicontrol('Parent',panel_option,'style','edit','Units','pixels',...
    'string','0.4','Enable','on','Callback',@edit_options_callback); 
    set(edit_decaytime,'Position',[137 45 58 19]);
text_lamprob = uicontrol('Parent',panel_option,'style','text','Units','pixels','string',...
    'lambda probability','Fontsize',12); 
    set(text_lamprob,'Position',[17 24 115 16]);
edit_lamprob = uicontrol('Parent',panel_option,'style','edit','Units','pixels',...
    'string','0.99','Enable','on','Callback',@edit_options_callback); 
    set(edit_lamprob,'Position',[137 22 58 19]);

% Ca time series plot
ax_F = axes('Parent',hdl_gui,'Units','pixels','Position',[50 553 901 90],'Xlim',[0 1],'Ylim',[0  1],'Box','off',...
            'Visible','on','clipping','off'); 

% df/f or dR/R plot
ax_dF = axes('Parent',hdl_gui,'Units','pixels','Position',[50 424 901 90],'Xlim',[0 1],'Ylim',[0  1],'Box','off',...
            'Visible','on','clipping','off');

% oasisAR1: spikes plot
ax_AR1 = axes('Parent',hdl_gui,'Units','pixels','Position',[50 293 901 90],'Xlim',[0 1],'Ylim',[0  1],'Box','off',...
            'Visible','on','clipping','off'); 

% oasisAR2: spikes plot
ax_AR2 = axes('Parent',hdl_gui,'Units','pixels','Position',[50 162 901 90],'Xlim',[0 1],'Ylim',[0  1],'Box','off',...
            'Visible','on','clipping','off'); 

    
%% initial GUI data
if ~isempty(dtsG)
    C1 = dtsG; str_C1 = 'FISSA-corrected timeseries';
else
    C1 = tsG; str_C1 = 'Raw timeseries';
end
if ~isempty(ddf_f)
    C2 = ddf_f; str_C2 = 'FISSA-corrected dF/F';
else
    if ~isempty(df_f)
        C2 = df_f; str_C2 = 'dF/F';
    else
        C2 = R; str_C2 = 'dR/R';
    end
end
setappdata(hdl_gui,'curr_C2',C2);
calcC2;


%% GUI Subfunctions
function edit_cellNum_callback(varargin)
    id = str2double(get(edit_cellNum,'string'));

    if isnan(id)
        set(edit_cellNum,'string','1');
        set(slider_cellNum,'Value',1);
        showTS(1);
    else
        if id < 1 
            set(edit_cellNum,'string','1');
            set(slider_cellNum,'Value',1);
            showTS(1);
        elseif id > N
            set(edit_cellNum,'string',num2str(N));
            set(slider_cellNum,'Value',N);
            showTS(N);
        else
            set(slider_cellNum,'Value',id);
            showTS(id);
        end
    end
end

function slider_cellNum_callback(varargin)
    id = round(get(slider_cellNum,'Value'));
    set(edit_cellNum,'string',num2str(id));
    showTS(id);
end

function edit_prctile_callback(varargin)
    df_prctile = str2double(get(edit_prctile,'string'));

    if isnan(df_prctile)
        set(edit_prctile,'string','85');
        set(slider_prctile,'Value',85);
        calcC2;
    else
        if df_prctile < 0 
            set(edit_prctile,'string','0');
            set(slider_prctile,'Value',0);
            calcC2;
        elseif df_prctile > 100
            set(edit_prctile,'string','100');
            set(slider_prctile,'Value',100);
            calcC2;
        else
            set(slider_prctile,'Value',df_prctile);
            calcC2;
        end
    end
end

function slider_prctile_callback(varargin)
    df_prctile = round(get(slider_prctile,'Value'));
    set(edit_prctile,'string',num2str(df_prctile));
    calcC2;
end

function rbutton_autoprctile_callback(varargin)
    if get(rbutton_autoprctile,'Value')
        set(edit_prctile,'Enable','off');
        set(slider_prctile,'Enable','off');
        calcC2;
    else
        set(edit_prctile,'Enable','on');
        set(slider_prctile,'Enable','on');
        calcC2;
    end
end

function edit_options_callback(varargin)
    calcSpikes;
end

function calcC2(varargin)
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
        setappdata(hdl_gui,'auto_prctile',auto_prctile);
    end
    setappdata(hdl_gui,'curr_C2',curr_C2);
    calcSpikes;
end

function calcSpikes(varargin)
    C = getappdata(hdl_gui,'curr_C2');
    decay_time = str2double(get(edit_decaytime,'String'));
    lam_pr = str2double(get(edit_lamprob,'String'));
    spk_SNR = str2double(get(edit_spikeSNR,'String'));
    
    spikesAR1 = zeros(N,T);
    spikesAR2 = zeros(N,T);        
    for i = 1:N
        spkmin = spk_SNR*GetSn(C(i,:));
        lam = choose_lambda(exp(-1/(30.9*decay_time)),GetSn(C(i,:)),lam_pr);
        
        [~,spkAR1,~] = deconvolveCa(C(i,:),'ar1','method','thresholded','optimize_pars',true,'maxIter',20,...
                                    'window',150,'lambda',lam,'smin',spkmin);
        spikesAR1(i,:) = spkAR1(:);
        
        [~,spkAR2,~] = deconvolveCa(C(i,:),'ar2','method','thresholded','optimize_pars',true,'maxIter',20,...
                                    'window',150,'lambda',lam,'smin',spkmin);
        spikesAR2(i,:) = spkAR2(:);

        h = waitbar(i/N);
    end
    close(h);
    setappdata(hdl_gui,'curr_spikes_AR1',spikesAR1);
    setappdata(hdl_gui,'curr_spikes_AR2',spikesAR2);
    id = str2double(get(edit_cellNum,'string'));
    showTS(id);
end

function showTS(id)
    % Ca time series
    axes(ax_F);
    plot(C1(id,:)); title(str_C1);

    % df/f or dR/R
    curr_C2 = getappdata(hdl_gui,'curr_C2');
    axes(ax_dF);
    plot(curr_C2(id,:)); title(str_C2);

    % oasisAR1 spikes
    spikesAR1 = getappdata(hdl_gui,'curr_spikes_AR1');
    axes(ax_AR1);
    plot(spikesAR1(id,:)); title('Spikes (oasisAR1)');
    
    % oasisAR2 spikes
    spikesAR2 = getappdata(hdl_gui,'curr_spikes_AR2');
    axes(ax_AR2);
    plot(spikesAR2(id,:)); title('Spikes (oasisAR2)');
    
    df_prctile = getappdata(hdl_gui,'auto_prctile');
    if get(rbutton_autoprctile,'Value')
        set(edit_autoprctile,'String',num2str(df_prctile(id)));
    end
end

end