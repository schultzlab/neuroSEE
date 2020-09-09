function frun_imreg_global_batch(array_id, list, templateglob, templateloc, templateloccode, force)
if nargin<6, force = []; end

imregr_params = load([templateloccode 'toB_r.mat']);
    shifts = imregr_params.shifts;
    A = struct('A', repmat(shifts, 7420, 1));
    imregr_params.shifts = A.A;

imregnr_params = load([templateloccode 'toB_nr.mat']);
    shifts = imregnr_params.shifts;
    A = struct('A', repmat(shifts, 7420, 1));
    imregnr_params.shifts = A.A;

[~,~,err] = load_neuroSEEmodules;
if ~isempty(err)
    beep
    cprintf('Errors',err);    
    return
end

imreg_global_batch( array_id, list, templateglob, imregr_params, imregnr_params, templateloc, force );