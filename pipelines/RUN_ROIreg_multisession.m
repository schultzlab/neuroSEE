list = 'explist_m82_open_s1-2.txt';
ref_array = ['20190404_17_41_23';...
             '20190406_20_56_30'];
bl_prctile_array = [85; 85];
sessionind_array = [1; 3];
pfactivet_thr = 0.02;
force = false; 
figclose = false;
frun_ROIreg_multisession( 'explist_m82_open_s1-2.txt', ['20190404_17_41_23';'20190406_20_56_30'], [75; 75], 0.02, [1; 3], false, true )


% list = 'explist_m62_fam1_s1-4_complete.txt';
% ref_array = ['20181015_09_37_54';...
%              '20181013_10_53_51';...
%              '20181017_09_57_03';...
%              '20181016_09_09_43'];
% bl_array = [85; 90; 78; 78];
% sessionind_array = [2; 1; 4; 3];
% force = false; 
% figclose = false;
% 
% r = load('AtoB_r.mat','options');
% params_mc(1).r = r.options;
% nr = load('AtoB_nr.mat','options');
% params_mc(1).nr = nr.options;
% 
% r = load('DtoB_r.mat','options');
% params_mc(2).r = r.options;
% nr = load('DtoB_nr.mat','options');
% params_mc(2).nr = nr.options;
% 
% r = load('CtoB_r.mat','options');
% params_mc(3).r = r.options;
% nr = load('CtoB_nr.mat','options');
% params_mc(3).nr = nr.options;
% 
% params_mc_array = params_mc;
% 
% options = params.ROIreg;
% options_mc = params_mc_array;
% figname_pref = [];
% 
% params_mc_array(2).nr = params_mc_array(1).nr;