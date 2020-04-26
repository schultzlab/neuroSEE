% Written by Ann Go
%
% OUTPUTS
%   tsG     : raw timeseries
%   dtsG    : decontaminated raw timeseries
%   ddf_f   : (decontaminated) df_f

function [tsG, dtsG, ddf_f, params] = neuroSEE_neuropilDecon( masks, data_locn, file, params, force, list, reffile )

if nargin<7, reffile = []; end
if nargin<6, list = []; end
if nargin<5, force = 0; end

mcorr_method = params.methods.mcorr_method;
segment_method = params.methods.segment_method;

if isempty(list)
    tiffile = [data_locn,'Data/',file(1:8),'/Processed/',file,'/mcorr_',mcorr_method,'/',file,'_2P_XYT_green_mcorr.tif'];
    fissadir = [data_locn,'Data/',file(1:8),'/Processed/',file,'/mcorr_',mcorr_method,'/',segment_method,'/FISSA/'];
    
    fname_mat = [fissadir file '_fissa_output.mat'];
    fname_mat_temp = [fissadir 'FISSAout/matlab.mat'];
    fname_fig1 = [fissadir file '_fissa_result.fig'];
    fname_fig2 = [fissadir file '_fissa_df_f.fig'];
else
    [ mouseid, expname ] = find_mouseIDexpname(list);
    if strcmpi(file, reffile)
        tiffile = [data_locn 'Data/' file(1:8) '/Processed/' file '/mcorr_' mcorr_method '/' file '_2P_XYT_green_mcorr.tif'];
        fissadir = [data_locn 'Data/' file(1:8) '/Processed/' file '/mcorr_' mcorr_method '/' ...
                    segment_method '_' mouseid '_' expname '/FISSA/'];
    else
        imreg_method = params.methods.imreg_method;
        if strcmpi(imreg_method, mcorr_method)
            tifdir = [data_locn 'Data/' file(1:8) '/Processed/' file '/imreg_' imreg_method '_ref' reffile '/'];
        else
            tifdir = [data_locn 'Data/' file(1:8) '/Processed/' file '/imreg_' imreg_method '_ref' reffile '_' mcorr_method '/'];
        end
        tiffile = [tifdir file '_2P_XYT_green_imreg_ref' reffile '.tif'];
        fissadir = [tifdir segment_method '_' mouseid '_' expname '/FISSA/'];
    end
    fname_mat = [fissadir file '_' mouseid '_' expname '_ref' reffile '_fissa_output.mat'];
    fname_mat_temp = [fissadir 'FISSAout/matlab.mat'];
    fname_fig1 = [fissadir file '_' mouseid '_' expname '_ref' reffile '_fissa_result.fig'];
    fname_fig2 = [fissadir file '_' mouseid '_' expname '_ref' reffile '_fissa_df_f.fig'];
end

prevstr = [];
if force || ~exist(fname_mat,'file')
    if isempty(list)
        str = sprintf('%s: Doing FISSA correction...\n',file);
    else
        str = sprintf('%s: Doing FISSA correction...\n',[mouseid '_' expname '_' file]);
    end
    refreshdisp(str, prevstr);
    prevstr = str;
    
    if force || and( ~exist(fname_mat,'file'), ~exist(fname_mat_temp,'file') )
        runFISSA( masks, tiffile, fissadir );
    end
    raw = load(fname_mat_temp,'raw');
    result = load(fname_mat_temp,'result');
    
    % Convert raw timeseries cell array structure to a matrix
    tsG = zeros(size(masks,3),size(raw.raw.cell0.trial0,2));
    for i = 1:numel(fieldnames(raw.raw))
        tsG(i,:) = raw.raw.(['cell' num2str(i-1)]).trial0(1,:);
    end
    
    % Convert decontaminated timeseries cell array structure to a matrix
    dtsG = zeros(size(masks,3),size(result.result.cell0.trial0,2));
    for i = 1:numel(fieldnames(result.result))
        dtsG(i,:) = result.result.(['cell' num2str(i-1)]).trial0(1,:);
    end

    % Calculate df_f
    ddf_f = zeros(size(dtsG));
    ddf_prctile = params.fissa.ddf_prctile;
    for i = 1:size(dtsG,1)
        x = lowpass( medfilt1(dtsG(i,:),params.fissa.ddf_medfilt1), 1, params.fr );
        fo = ones(size(x)) * prctile(x,ddf_prctile);
        while fo == 0
            fo = ones(size(x)) * prctile(x,ddf_prctile+5);
            ddf_prctile = ddf_prctile+5;
        end
        ddf_f(i,:) = (x - fo) ./ fo;
    end
    
    % Save output
    output.tsG = tsG;
    output.dtsG = dtsG;
    output.ddf_f = ddf_f;
    output.params = params.fissa;
    if ~exist( fissadir, 'dir' ), mkdir( fissadir ); end
    save(fname_mat,'-struct','output');
    if isempty(list)
        str = sprintf('%s: FISSA correction done\n',file);
    else
        str = sprintf('%s: FISSA correction done\n',[mouseid '_' expname  '_' file]);
    end
    refreshdisp(str, prevstr);
    
    % plot 
    makeplot(dtsG, ddf_f);
else
    % If it exists, load it 
    fissa_output = load(fname_mat);
    if isfield(fissa_output,'tsG'), tsG = fissa_output.tsG; end
    dtsG = fissa_output.dtsG;
    ddf_f = fissa_output.ddf_f;
    params.fissa = fissa_output.params;

    if ~exist(fname_fig1,'file') || ~exist(fname_fig2,'file')
        makeplot(dtsG, ddf_f);
    end
    if isempty(list)
        fprintf('%s: Neuropil decontamination output found and loaded\n',file);
    else
        fprintf('%s: Neuropil decontamination output found and loaded\n',[mouseid '_' expname  '_' file]);
    end
end

function makeplot(dtsG, ddf_f)
    % raw timeseries
    fig = figure;
    iosr.figures.multiwaveplot(1:size(dtsG,2),1:size(dtsG,1),dtsG,'gain',5); yticks([]); xticks([]); 
    title('Fissa-corrected raw timeseries','Fontweight','normal','Fontsize',12); 
    savefig(fig, fname_fig1(1:end-4));
    saveas(fig, fname_fig1(1:end-4),'png');
    close(fig);

    % dF/F
    fig = figure;
    iosr.figures.multiwaveplot(1:size(ddf_f,2),1:size(ddf_f,1),ddf_f,'gain',5); yticks([]); xticks([]); 
    title('Fissa-corrected dF/F','Fontweight','normal','Fontsize',12); 
    savefig(fig, fname_fig2(1:end-4));
    saveas(fig, fname_fig2(1:end-4),'png');
    close(fig);
end

end

