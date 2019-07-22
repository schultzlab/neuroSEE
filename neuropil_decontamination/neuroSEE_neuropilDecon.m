% Written by Ann Go
%
% OUTPUTS
%   tsG     : raw timeseries
%   dtsG    : decontaminated raw timeseries
%   ddf_f   : (decontaminated) df_f

function [tsG, dtsG, ddf_f, params] = neuroSEE_neuropilDecon( masks, data_locn, file, params, force )

if nargin<5, force = 0;      end

mcorr_method = params.methods.mcorr_method;
segment_method = params.methods.segment_method;
tiffile = [data_locn,'Data/',file(1:8),'/Processed/',file,'/mcorr_',mcorr_method,'/',file,'_2P_XYT_green_mcorr.tif'];
fissadir = [data_locn,'Data/',file(1:8),'/Processed/',file,'/mcorr_',mcorr_method,'/',segment_method,'/FISSA/'];
if ~exist( fissadir, 'dir' ), mkdir( fissadir ); end

fname_mat = [fissadir file '_fissa_output.mat'];
fname_mat_temp = [fissadir 'FISSAout/matlab.mat'];

prevstr = [];
if force || ~exist(fname_mat,'file')
    str = sprintf('%s: Doing FISSA correction...\n',file);
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
    for i = 1:size(dtsG,1)
        x = lowpass( medfilt1(dtsG(i,:),params.fissa.ddf_medfilt1), 1, params.fr );
        fo = ones(size(x)) * prctile(x,params.fissa.ddf_prctile);
        ddf_f(i,:) = (x - fo) ./ fo;
    end
    
    output.tsG = tsG;
    output.dtsG = dtsG;
    output.ddf_f = ddf_f;
    output.params = params.fissa;
    save(fname_mat,'-struct','output');
    str = sprintf('%s: FISSA correction done\n',file);
    refreshdisp(str, prevstr);

else
    % If it exists, load it 
    fissa_output = load(fname_mat);
    tsG = fissa_output.tsG;
    dtsG = fissa_output.dtsG;
    ddf_f = fissa_output.ddf_f;
    params.fissa = fissa_output.params;

    fprintf('%s: Neuropil decontamination output found and loaded\n',file);
end

