% Written by Ann Go
%
% INPUTS
%   conc_runs  : flag if rois were segmented from concatenated files from
%               different runs e.g. fam1fam2fam1-fam1 but rois are
%               for fam1fam2fam1. DO NOT flag for fam1fam2fam1 files since in
%               this case it is understood that the rois are from the
%               concatenated runs.
% OUTPUTS
%   tsG     : raw timeseries
%   dtsG    : decontaminated raw timeseries
%   ddf_f   : (decontaminated) df_f


function [dtsG, ddf_f, params] = neuroSEE_neuropilDecon( masks, data_locn, params, fileorlist, reffile, conc_runs, grp_sdir, force )

if nargin<8, force = 0; end
if nargin<7, list = []; end
if nargin<6, conc_runs = false; end
if nargin<5, reffile = []; end

mcorr_method = params.methods.mcorr_method;
segment_method = params.methods.segment_method;

% determine whether fileorlist is a file or list
if ~strncmp(fileorlist,'list',4)
    file = fileorlist; list = [];
else
    list = fileorlist; file = [];
end

if isempty(list) % file
    tiff = [data_locn,'Data/',file(1:8),'/Processed/',file,'/mcorr_',mcorr_method,'/',file,'_2P_XYT_green_mcorr.tif'];
    fissadir = [data_locn,'Data/',file(1:8),'/Processed/',file,'/mcorr_',mcorr_method,'/',segment_method,'/FISSA/'];
    
    fname_mat = [fissadir file '_fissa_output.mat'];
    fname_mat_temp = [fissadir 'FISSAout/matlab.mat'];
    fname_fig1 = [fissadir file '_fissa_result.fig'];
    fname_fig2 = [fissadir file '_fissa_df_f.fig'];
else % list
    [ mouseid, expname ] = find_mouseIDexpname(list);
    if conc_runs
        concrunsname = find_concrunsname( list );
        expname = concrunsname;
    end
    listfile = [data_locn 'Digital_Logbook/lists_imaging/' list];
    files = extractFilenamesFromTxtfile( listfile );
    Nfiles = size(files,1);
    tiff = [];
    for n = 1:N
        if strcmpi(file, reffile)
            tiffile = [data_locn 'Data/' file(1:8) '/Processed/' file '/mcorr_' mcorr_method '/' file '_2P_XYT_green_mcorr.tif'];
        else
            tifdir = [data_locn 'Data/' file(1:8) '/Processed/' file '/imreg_' mcorr_method '_ref' reffile '/'];
            tiffile = [tifdir file '_2P_XYT_green_imreg_ref' reffile '.tif'];
        end
        if isempty(tiff)
            tiff = tiffile;
        else
            tiff = [tiff ',' tiffile];
        end
    end
    fissadir = [grp_sdir '/' str_fissa '/'];
    fname_mat = [fissadir mouseid '_' expname '_ref' reffile '_fissa_output.mat'];
    fname_mat_temp = [fissadir 'FISSAout/matlab.mat'];
    fname_fig1 = [fissadir mouseid '_' expname '_ref' reffile '_fissa_result.fig'];
    fname_fig2 = [fissadir mouseid '_' expname '_ref' reffile '_fissa_df_f.fig'];
end

prevstr = [];
if force || ~exist(fname_mat,'file')
    if isempty(list)
        str = sprintf('%s: Doing FISSA correction...\n',file);
    else
        str = sprintf('%s: Doing FISSA correction...\n',[mouseid '_' expname]);
    end
    refreshdisp(str, prevstr);
    prevstr = str;
    
    if force || and( ~exist(fname_mat,'file'), ~exist(fname_mat_temp,'file') )
        runFISSA( masks, tiff, fissadir );
    end
    
    result = load(fname_mat_temp,'result');
    df_result = load(fname_mat_temp,'deltaf_result');

    % Convert decontaminated timeseries cell array structure to a matrix
    T = size(result.result.cell0.trial0,2);
    dtsG = zeros(size(masks,3),T);
    for i = 1:numel(fieldnames(result.result))
        for j = 1:numel(fieldnames(result.result.cell0))
            dtsG(i,(j-1)*T+1:j*T) = result.result.(['cell' num2str(i-1)]).(['trial' num2str(j-1)])(1,:);
        end
    end

    ddf_f = zeros(size(masks,3),T);
    for i = 1:numel(fieldnames(df_result.deltaf_result))
        for j = 1:numel(fieldnames(df_result.deltaf_result.cell0))
            ddf_f(i,(j-1)*T+1:j*T) = lowpass( medfilt1(df_result.deltaf_result.(['cell' num2str(i-1)]).(['trial' num2str(j-1)])(1,:), params.fissa.ddf_medfilt1), 1, 30.9 );
        end
    end

    % Save output
    output.dtsG = dtsG;
    output.ddf_f = ddf_f;
    output.params = params.fissa;
    if isempty(list)
        str = sprintf('%s: Saving fissa output\n',file);
    else
        str = sprintf('%s: Saving fissa output\n',[mouseid '_' expname]);
    end
    refreshdisp(str, prevstr);
    
    if ~exist( fissadir, 'dir' ), mkdir( fissadir ); fileattrib(fissadir,'+w','g','s'); end
    save(fname_mat,'-struct','output');
    
    if isempty(list)
        str = sprintf('%s: FISSA correction done\n',file);
    else
        str = sprintf('%s: FISSA correction done\n',[mouseid '_' expname]);
    end
    refreshdisp(str, prevstr);
        
    % plot timeseries of 50 cells (for visibility)
    multiplot_ts(dtsG, fname_fig1(1:end-4), 'Fissa-corrected raw timeseries', true, 50);
    multiplot_ts(ddf_f, fname_fig2(1:end-4), 'Fissa-corrected dF/F', true, 50);
else
    % If it exists, load it 
    str = sprintf('%s: Loading fissa data\n', [mouseid '_' expname]);
    cprintf(str);
    fissa_output = load(fname_mat);
    dtsG = fissa_output.dtsG;
    ddf_f = fissa_output.ddf_f;
    params.fissa = fissa_output.params;

    newstr = sprintf('%s: Fissa data loaded\n', [mouseid '_' expname]);
    refreshdisp(newstr, str)
            
    if ~exist(fname_fig1,'file') || ~exist(fname_fig2,'file')
        multiplot_ts(dtsG, fname_fig1(1:end-4), 'Fissa-corrected raw timeseries');
        multiplot_ts(ddf_f, fname_fig2(1:end-4), 'Fissa-corrected dF/F');
    end
    if isempty(list)
        fprintf('%s: Neuropil decontamination output found and loaded\n',file);
    else
        fprintf('%s: Neuropil decontamination output found and loaded\n',[mouseid '_' expname]);
    end
end

end

