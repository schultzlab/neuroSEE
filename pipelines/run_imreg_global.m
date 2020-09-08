file = '20181016_09_09_43';
templateglob = '20181015_09_37_54';
imregr_params = load('CtoB-r.mat');
imregnr_params = load('CtoB-nr.mat');
imreg_global( file, templateglob, imregr_params, imregnr_params )