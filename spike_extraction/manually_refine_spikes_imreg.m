% in /neuroSEE/pipelines
test = false;                   % flag to use one of smaller files in test folder)
[data_locn,comp,err] = load_neuroSEEmodules(test);
if ~isempty(err)
    beep
    cprintf('Errors',err);    
    return
end
params = load('default_params.mat');
file = '20190406_20_14_42';
reffile = '20190406_20_27_07';
expname = 'm82_open_day2';

addpath(genpath('../spike_extraction'));

if ~strcmpi(file,reffile)
    fissadir = [data_locn 'Data/' file(1:8) '/Processed/' file '/mcorr_normcorre_ref' reffile '/CaImAn_' expname '/FISSA/'];
    load([fissadir file '_' expname '_fissa_output.mat' ]);
    load([fissadir file '_' expname '_spikes.mat' ]);
    GUI_manually_refine_spikes_imreg_notgeneral( spikes, [], dtsG, [], ddf_f, data_locn, file, reffile, expname, params, [], [] )
else
    fissadir = [data_locn 'Data/' file(1:8) '/Processed/' file '/mcorr_normcorre/CaImAn_' expname '/FISSA/'];
    load([fissadir file '_' expname '_fissa_output.mat' ]);
    load([fissadir file '_' expname '_spikes.mat' ]);
    GUI_manually_refine_spikes_imreg_notgeneral( spikes, [], dtsG, [], ddf_f, data_locn, file, reffile, expname, params, [], [] )
end
