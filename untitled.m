frun_pipeline_imreg( 'list_m62_fov2_fam1fam2-fam2.txt', '20181013_10_53_51', true, [0;0;0;0;0;0], [1;1;1;1;1;1], 5, 88, 99, 0.05, 0.5, 3.0, true )
list = 'list_m62_fov2_fam1fam2-fam1.txt'; 
reffile = '20181013_10_53_51'; 
dofissa = true; 
force = [0;0;0;1;0;0];
dostep = [1;1;1;1;1;1]; 
tsub = 5; 
bl_prctile = 85; 
prctile_thr = 99; 
pfactivet_thr = 0.03; 
activetrials_thr = 0.3; 
min_SNR = 3.0; 
conc_env = true;

