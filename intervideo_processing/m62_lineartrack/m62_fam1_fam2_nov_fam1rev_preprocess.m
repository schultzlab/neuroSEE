clear; close all;

%% Basic setup
addpath(genpath('../behaviour'));
addpath(genpath('../intervideo_processing'));
addpath(genpath('../motion_correction'));
addpath(genpath('../PF_mapping'));
addpath(genpath('../pipelines'));
addpath(genpath('../ROI_segmentation'));
addpath(genpath('../spike_extraction'));
addpath(genpath('../utilities'));

% Data location
data_locn = '/Volumes/thefarm2/live/CrazyEights/AD_2PCa/';
if ~exist(data_locn,'dir')
    data_locn = '/Volumes/RDS/project/thefarm2/live/CrazyEights/AD_2PCa/';
end
if ~exist(data_locn,'dir')
    data_locn = '/rds/general/user/mgo/projects/thefarm2/live/CrazyEights/AD_2PCa/';
end

data_locn1 = [data_locn 'Data/20181013/Processed/'];
data_locn2 = [data_locn 'Data/20181015/Processed/'];
data_locn3 = [data_locn 'Data/20181016/Processed/'];
data_locn4 = [data_locn 'Data/20181017/Processed/'];


%% Load images and ompare red images to see which ones are most similar
load(strcat(data_locn1,'20181015_09_26_48/20181015_09_26_48_2P_mcorr_output.mat'))
fam1d2_0926.green = green;
fam1d2_0926.red = red;

load(strcat(data_locn1,'20181015_09_26_48/20181015_09_26_48_2P_mcorr_output.mat'))
fam1d2_0926.green = green;
fam1d2_0926.red = red;

load(strcat(data_locn1,'20181015_09_33_26/20181015_09_33_26_2P_mcorr_output.mat'))
fam1d2_0933.green = green;
fam1d2_0933.red = red;

load(strcat(data_locn1,'20181015_09_37_54/20181015_09_37_54_2P_mcorr_output.mat'))
fam1d2_0937.green = green;
fam1d2_0937.red = red;

load(strcat(data_locn1,'20181015_09_42_24/20181015_09_42_24_2P_mcorr_output.mat'))
fam1d2_0942.green = green;
fam1d2_0942.red = red;

load(strcat(data_locn1,'20181015_09_46_53/20181015_09_46_53_2P_mcorr_output.mat'))
fam1d2_0946.green = green;
fam1d2_0946.red = red;

load(strcat(data_locn1,'20181015_09_51_30/20181015_09_51_30_2P_mcorr_output.mat'))
fam1d2_0951.green = green;
fam1d2_0951.red = red;

load(strcat(data_locn1,'20181015_09_59_57/20181015_09_59_57_2P_mcorr_output.mat'))
novd2_0959.green = green;
novd2_0959.red = red;

load(strcat(data_locn1,'20181015_10_04_26/20181015_10_04_26_2P_mcorr_output.mat'))
novd2_1004.green = green;
novd2_1004.red = red;

load(strcat(data_locn1,'20181015_10_08_57/20181015_10_08_57_2P_mcorr_output.mat'))
novd2_1008.green = green;
novd2_1008.red = red;

load(strcat(data_locn2,'20181016_09_09_43/20181016_09_09_43_2P_mcorr_output.mat'))
fam1d3_0909.green = green;
fam1d3_0909.red = red;

load(strcat(data_locn2,'20181016_09_14_03/20181016_09_14_03_2P_mcorr_output.mat'))
fam1d3_0914.green = green;
fam1d3_0914.red = red;

load(strcat(data_locn2,'20181016_09_18_26/20181016_09_18_26_2P_mcorr_output.mat'))
fam1d3_0918.green = green;
fam1d3_0918.red = red;

load(strcat(data_locn2,'20181016_09_22_45/20181016_09_22_45_2P_mcorr_output.mat'))
fam1d3_0922.green = green;
fam1d3_0922.red = red;

load(strcat(data_locn2,'20181016_09_44_06/20181016_09_44_06_2P_mcorr_output.mat'))
fam1d3_0944.green = green;
fam1d3_0944.red = red;

load(strcat(data_locn2,'20181016_09_49_06/20181016_09_49_06_2P_mcorr_output.mat'))
fam1d3_0949.green = green;
fam1d3_0949.red = red;

load(strcat(data_locn2,'20181016_09_57_43/20181016_09_57_43_2P_mcorr_output.mat'))
fam1revd3_0957.green = green;
fam1revd3_0957.red = red;

load(strcat(data_locn2,'20181016_10_07_07/20181016_10_07_07_2P_mcorr_output.mat'))
fam1revd3_1007.green = green;
fam1revd3_1007.red = red;

load(strcat(data_locn2,'20181016_10_11_35/20181016_10_11_35_2P_mcorr_output.mat'))
fam1revd3_1011.green = green;
fam1revd3_1011.red = red;

load(strcat(data_locn2,'20181016_10_15_57/20181016_10_15_57_2P_mcorr_output.mat'))
fam1revd3_1015.green = green;
fam1revd3_1015.red = red;
clear green red params

% figure;
% subplot(451);imagesc(fam1d2_0926.red.meanregframe); title('Fam1 d2 09:26 (8)'); caxis([2000 2400]); axis off; axis square;
% subplot(452);imagesc(fam1d2_0933.red.meanregframe); title('Fam1 d2 09:33 (12)'); caxis([2000 2400]); axis off; axis square;
% subplot(453);imagesc(fam1d2_0937.red.meanregframe); title('Fam1 d2 09:37 (11)'); caxis([2000 2400]); axis off; axis square;
% subplot(454);imagesc(fam1d2_0942.red.meanregframe); title('Fam1 d2 09:42 (9)'); caxis([2000 2400]); axis off; axis square;
% subplot(455);imagesc(fam1d2_0946.red.meanregframe); title('Fam1 d2 09:46 (17)'); caxis([2000 2400]); axis off; axis square;
% subplot(456);imagesc(fam1d2_0946.red.meanregframe); title('Fam1 d2 09:51 (14)'); caxis([2000 2400]); axis off; axis square;
% subplot(457);imagesc(fam1d2_0959.red.meanregframe); title('Nov d2 09:59 (27)'); caxis([2000 2400]); axis off; axis square;
% subplot(458);imagesc(fam1d2_1004.red.meanregframe); title('Nov d2 10:04 (26)'); caxis([2000 2400]); axis off; axis square;
% subplot(459);imagesc(fam1d2_1008.red.meanregframe); title('Nov d2 10:08 (28)'); caxis([2000 2400]); axis off; axis square;
% 
% subplot(4,5,11);imagesc(fam1d3_0909.red.meanregframe); title('Fam1 d3 09:09 (21)'); caxis([2000 2400]); axis off; axis square;
% subplot(4,5,12);imagesc(fam1d3_0914.red.meanregframe); title('Fam1 d3 09:14 (17)'); caxis([2000 2400]); axis off; axis square;
% subplot(4,5,13);imagesc(fam1d3_0918.red.meanregframe); title('Fam1 d3 09:18 (21)'); caxis([2000 2400]); axis off; axis square;
% subplot(4,5,14);imagesc(fam1d3_0922.red.meanregframe); title('Fam1 d3 09:22 (19)'); caxis([2000 2400]); axis off; axis square;
% subplot(4,5,15);imagesc(fam1d3_0944.red.meanregframe); title('Fam1 d3 09:44 (13)'); caxis([2000 2400]); axis off; axis square;
% subplot(4,5,16);imagesc(fam1d3_0949.red.meanregframe); title('Fam1 d3 09:49 (18)'); caxis([2000 2400]); axis off; axis square;
% subplot(4,5,17);imagesc(fam1revd3_0957.red.meanregframe); title('Fam1rev d3 09:57 (20)'); caxis([2000 2400]); axis off; axis square;
% subplot(4,5,18);imagesc(fam1revd3_1007.red.meanregframe); title('Fam1rev d3 10:07 (35)'); caxis([2000 2400]); axis off; axis square;
% subplot(4,5,19);imagesc(fam1revd3_1011.red.meanregframe); title('Fam1rev d3 10:11 (46)'); caxis([2000 2400]); axis off; axis square;
% subplot(4,5,20);imagesc(fam1revd3_1015.red.meanregframe); title('Fam1rev d3 10:15 (39)'); caxis([2000 2400]); axis off; axis square;
 

%% Choose a reference image and create a superset of video frames for ROI segmentation

% reference image: fam1d2_0946

% novd2_1004
cprintf('novd2_1004: calculating shift...\n')
[ novd2_1004.shift, novd2_1004.red.globalreg_meanframe ] = globalregisterImage(fam1d2_0946.red.meanregframe, novd2_1004.red.meanregframe, 1 );
cprintf('novd2_1004: reading tif...\n')
novd2_1004_imG = read_file( [data_locn1 '20181015_10_04_26/20181015_10_04_26_2P_XYT_green_mcorr.tif'] );
novd2_1004_imR = read_file( [data_locn1 '20181015_10_04_26/20181015_10_04_26_2P_XYT_red_mcorr.tif'] );

cprintf('novd2_1004: registering images...\n')
[ novd2_1004_imGglobalreg, novd2_1004_imRglobalreg ] = ...
    globalregisterStack( novd2_1004_imG, novd2_1004_imR, novd2_1004.shift(1) , novd2_1004.shift(2) );
clear novd2_1004_imG novd2_1004_imR

cprintf('novd2_1004: saving registered tif images...\n')
fname_tif_novd2_1004_imGglobalreg = [data_locn1 '20181015_10_04_26/20181015_10_04_26_2P_XYT_green_mcorr_globalreg.tif'];
    writeTifStack( novd2_1004_imGglobalreg,fname_tif_novd2_1004_imGglobalreg );
fname_tif_novd2_1004_imRglobalreg = [data_locn1 '20181015_10_04_26/20181015_10_04_26_2P_XYT_red_mcorr_globalreg.tif'];
    writeTifStack( novd2_1004_imRglobalreg,fname_tif_novd2_1004_imRglobalreg );
clear novd2_1004_imGglobalreg novd2_1004_imRglobalreg fname_tif_novd2_1004_imGglobalreg fname_tif_novd2_1004_imRglobalreg
            
% fam1d3_0922
cprintf('fam1d3_0922: calculating shift...\n')
[ fam1d3_0922.shift, fam1d3_0922.red.globalreg_meanframe ] = globalregisterImage(fam1d2_0946.red.meanregframe, fam1d3_0922.red.meanregframe, 1 );
cprintf('fam1d3_0922: reading tif...\n')
fam1d3_0922_imG = read_file( [data_locn2 '20181016_09_22_45/20181016_09_22_45_2P_XYT_green_mcorr.tif'] );
fam1d3_0922_imR = read_file( [data_locn2 '20181016_09_22_45/20181016_09_22_45_2P_XYT_red_mcorr.tif'] );

cprintf('fam1d3_0922: registering images...\n')
[ fam1d3_0922_imGglobalreg, fam1d3_0922_imRglobalreg ] = ...
    globalregisterStack( fam1d3_0922_imG, fam1d3_0922_imR, fam1d3_0922.shift(1) , fam1d3_0922.shift(2) );
clear fam1d3_0922_imG fam1d3_0922_imR

cprintf('fam1d3_0922: saving registered tif images...\n')
fname_tif_fam1d3_0922_imGglobalreg = [data_locn2 '20181016_09_22_45/20181016_09_22_45_2P_XYT_green_mcorr_globalreg.tif'];
    writeTifStack( fam1d3_0922_imGglobalreg,fname_tif_fam1d3_0922_imGglobalreg );
fname_tif_fam1d3_0922_imRglobalreg = [data_locn2 '20181016_09_22_45/20181016_09_22_45_2P_XYT_red_mcorr_globalreg.tif'];
    writeTifStack( fam1d3_0922_imRglobalreg,fname_tif_fam1d3_0922_imRglobalreg );
clear fam1d3_0922_imGglobalreg fam1d3_0922_imRglobalreg fname_tif_fam1d3_0922_imGglobalreg fname_tif_fam1d3_0922_imRglobalreg

% fam1revd3_1007
cprintf('fam1revd3_1007: calculating shift...\n')
[ fam1revd3_1007.shift, fam1revd3_1007.red.globalreg_meanframe ] = globalregisterImage(fam1d2_0946.red.meanregframe, fam1revd3_1007.red.meanregframe, 1 );
cprintf('fam1revd3_1007: reading tif...\n')
fam1revd3_1007_imG = read_file( [data_locn2 '20181016_10_07_07/20181016_10_07_07_2P_XYT_green_mcorr.tif'] );
fam1revd3_1007_imR = read_file( [data_locn2 '20181016_10_07_07/20181016_10_07_07_2P_XYT_red_mcorr.tif'] );

cprintf('fam1revd3_1007: registering images...\n')
[ fam1revd3_1007_imGglobalreg, fam1revd3_1007_imRglobalreg ] = ...
    globalregisterStack( fam1revd3_1007_imG, fam1revd3_1007_imR, fam1revd3_1007.shift(1) , fam1revd3_1007.shift(2) );
clear fam1revd3_1007_imG fam1revd3_1007_imR

cprintf('fam1revd3_1007: saving registered tif images...\n')
fname_tif_fam1revd3_1007_imGglobalreg = [data_locn2 '20181016_10_07_07/20181016_10_07_07_2P_XYT_green_mcorr_globalreg.tif'];
    writeTifStack( fam1revd3_1007_imGglobalreg,fname_tif_fam1revd3_1007_imGglobalreg );
fname_tif_fam1revd3_1007_imRglobalreg = [data_locn2 '20181016_10_07_07/20181016_10_07_07_2P_XYT_red_mcorr_globalreg.tif'];
    writeTifStack( fam1revd3_1007_imRglobalreg,fname_tif_fam1revd3_1007_imRglobalreg );
clear fam1revd3_1007_imGglobalreg fam1revd3_1007_imRglobalreg fname_tif_fam1revd3_1007_imGglobalreg fname_tif_fam1revd3_1007_imRglobalreg

clear data_locn data_locn1 data_locn2
save('fam1_nov_fam1rev_20181015to16_m62.mat')
    
