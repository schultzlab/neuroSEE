function GUI_viewImReg(file, mcorr_method, ref_array)

hdl_gui = figure('Name','Image registration result','NumberTitle','off','Resize','off',...
    'Position',[1268 1044 880 644]);

ax_image = axes('Units','pixels','Position',[20 23 600 600],'Xlim',[0 1],'Ylim',[0  1],'Box','off',...
    'Visible','off','clipping','off'); 
text_regref = uicontrol('Parent',hdl_gui,'style','text','Units','pixels','string','Registration reference','Fontsize',16,...
    'FontWeight','bold'); 
    set(text_regref,'Position',[645 594 201 31]);
list_regref = uicontrol('Parent',hdl_gui,'style','list','Units','pixels','string',ref_array,'Fontsize',12,...
    'Max',2,'Min',1,'Value',1,'Callback',@list_regref_callback); 
    set(list_regref,'Position',[650 493 201 101]);
text_colour = uicontrol('Parent',hdl_gui,'style','text','Units','pixels','string','Show reference in','Fontsize',16,...
    'FontWeight','bold'); 
    set(text_colour,'Position',[635 438 201 31]);
check_red = uicontrol('Parent',hdl_gui,'style','checkbox','Units','pixels','string','red','Fontsize',12,...
    'ForegroundColor','r','Value',1,'Callback',@check_red_callback); 
    set(check_red,'Position',[668 412 101 32]);
check_green = uicontrol('Parent',hdl_gui,'style','checkbox','Units','pixels','string','green','Fontsize',12,...
    'ForegroundColor',[0.4 0.6 0.2],'Value',0,'Callback',@check_green_callback); 
    set(check_green,'Position',[668 390 101 32]);

% load image file
[data_locn,comp,err] = load_neuroSEEmodules;
if ~isempty(err)
    beep
    cprintf('Errors',err);    
    return
end
fdir = [ data_locn 'Data/' file(1:8) '/Processed/' file '/imreg_' mcorr_method '_ref' ];
file_g = zeros(512,512,size(ref_array,2));
for n = 1:size(ref_array,1)
    fname_mat_file = [fdir ref_array(n,:) '/' file '_imreg_ref' ref_array(n,:) '_output.mat'];
    green = load(fname_mat_file,'green');
    file_g(:,:,n) = green.green.meanregframe;
end

% load ref files
template_g = zeros(512,512,size(ref_array,2));
for n = 1:size(ref_array,1)
    tdir = [ data_locn 'Data/' ref_array(n,1:8) '/Processed/' ref_array(n,:) '/mcorr_' mcorr_method '/' ];
    fname_mat_template = [tdir ref_array(n,:) '_mcorr_output.mat'];
    green = load(fname_mat_template,'green');
    template_g(:,:,n) = green.green.meanregframe;
end

showImages(1);
% Callback functions
    function showImages(refnum)
        checkred = get(check_red,'Value');
        axes(ax_image);
        if checkred
            C1 = imfuse( file_g(:,:,refnum), template_g(:,:,refnum), 'falsecolor', 'Scaling', 'joint', 'ColorChannels', [2 1 0]);
        else
            C1 = imfuse( file_g(:,:,refnum), template_g(:,:,refnum), 'falsecolor', 'Scaling', 'joint', 'ColorChannels', [1 2 0]);
        end
        imshow(C1);  
    end

    function list_regref_callback(varargin)
        refnum = get(list_regref,'Value');
        showImages(refnum)
    end

    function check_red_callback(varargin)
        checkred = get(check_red,'Value');
        if checkred
            set(check_green,'Value',0);
        else
            set(check_green,'Value',1);
        end
        refnum = get(list_regref,'Value');
        showImages(refnum)
    end

    function check_green_callback(varargin)
        checkgreen = get(check_green,'Value');
        if checkgreen
            set(check_red,'Value',0);
        else
            set(check_red,'Value',1);
        end
        refnum = get(list_regref,'Value');
        showImages(refnum)
    end
end