file = '20181013_10_43_47';
templateglob = '20181015_09_37_54';
templateloc = '20181013_10_53_51';
imregr_params = load('AtoB_r.mat');
shifts = imregr_params.shifts;
    A = struct('A', repmat(shifts, 7420, 1));
    imregr_params.shifts = A.A;
imregnr_params = load('AtoB_nr.mat');
    shifts = imregnr_params.shifts;
    A = struct('A', repmat(shifts, 7420, 1));
    imregnr_params.shifts = A.A;
    
[data_locn,~,err] = load_neuroSEEmodules;
if ~isempty(err)
    beep
    cprintf('Errors',err);    
    return
end

imreg_global( file, templateglob, imregr_params, imregnr_params, templateloc );
