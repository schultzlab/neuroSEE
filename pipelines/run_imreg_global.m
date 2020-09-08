file = '20181016_09_09_43';
templateglob = '20181015_09_37_54';
imregr_params = load('CtoB_r.mat');
imregnr_params = load('CtoB_nr.mat');

[data_locn,~,err] = load_neuroSEEmodules;
if ~isempty(err)
    beep
    cprintf('Errors',err);    
    return
end

imreg_global( file, templateglob, imregr_params, imregnr_params )
