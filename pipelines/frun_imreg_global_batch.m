function frun_imreg_global_batch(array_id, list, templateglob, templateloc, )
list = 'list_m62_fam1fam2_s2-fam1.txt';
templateglob = '20181015_09_37_54';
imregr_params = load('CtoB_r.mat');
imregnr_params = load('CtoB_nr.mat');
imreg_global( array_id, list, templateglob, imregr_params, imregnr_params, templateloc, false )