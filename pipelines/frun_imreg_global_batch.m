function frun_imreg_global_batch(array_id, list, templateglob, templateloc, shiftsfile_prefix, force)
if nargin<6, force = []; end

imregr_params = load([shiftsfile_prefix '_r.mat']);
    shifts = imregr_params.params{array_id}.shifts;
    A = struct('A', repmat(shifts, 7420, 1));
    imregr_params.params{array_id}.shifts = A.A;

imregnr_params = load([shiftsfile_prefix '_nr.mat']);
    shifts = imregnr_params.params{array_id}.shifts;
    A = struct('A', repmat(shifts, 7420, 1));
    imregnr_params.params{array_id}.shifts = A.A;

[~,~,err] = load_neuroSEEmodules;
if ~isempty(err)
    beep
    cprintf('Errors',err);    
    return
end

imreg_global_batch( array_id, list, templateglob, imregr_params.params{array_id}, imregnr_params.params{array_id}, templateloc, force );