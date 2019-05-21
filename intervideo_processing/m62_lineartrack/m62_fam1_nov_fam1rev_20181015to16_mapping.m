clear; close all;
tic

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

filedir1 = [data_locn 'Data/20181015/Processed/'];
filedir2 = [data_locn 'Data/20181016/Processed/'];
filedir3 = [data_locn 'Summaries/m62/fam1_nov_fam1rev/'];

file_ref = '20181015_09_46_53';
% Files to be registered
files_fam1d2 = ['20181015_09_33_26';...
                '20181015_09_37_54';...
                '20181015_09_51_30'];
files_novd2  = ['20181015_09_59_57';...
                '20181015_10_04_26'];
files_fam1d3 = ['20181016_09_22_45';...
                '20181016_09_44_06';...
                '20181016_09_49_06'];
files_fam1revd3 = ['20181016_09_57_43';...
                '20181016_10_07_07'];
     
% Read the masks
load([filedir3 'subsampledSplicedStacks_segment_output.mat'])
clear params
Nmasks = size(masks,3);


%% Load images
load('m62_fam1_nov_fam1rev_20181015to16_preprocess.mat')

% fam1 day2
fam1d2(1).green = fam1d2_0933.green;
fam1d2(1).red   = fam1d2_0933.red;

fam1d2(2).green = fam1d2_0937.green;
fam1d2(2).red   = fam1d2_0937.red;

fam1d2(3).green = fam1d2_0951.green;
fam1d2(3).red   = fam1d2_0951.red;

% nov day2
novd2(1).green = novd2_0959.green;
novd2(1).red   = novd2_0959.red;

novd2(2).green = novd2_1004.green;
novd2(2).red   = novd2_1004.red;

% fam1 day3
fam1d3(1).green = fam1d3_0922.green;
fam1d3(1).red   = fam1d3_0922.red;

fam1d3(2).green = fam1d3_0944.green;
fam1d3(2).red   = fam1d3_0944.red;

fam1d3(3).green = fam1d3_0949.green;
fam1d3(3).red   = fam1d3_0949.red;

% fam1rev day3
fam1revd3(1).green = fam1revd3_0957.green;
fam1revd3(1).red   = fam1revd3_0957.red;

fam1revd3(2).green = fam1revd3_1007.green;
fam1revd3(2).red   = fam1revd3_1007.red;


%% Reference image
str = sprintf('Reading reference image: %s\n', file_ref);
cprintf(str)

ref = fam1d2_0946;
ref.tsG = cell_tsG;
ref.tsR = cell_tsR;
clear cell_tsG cell_tsR

ref.R = ratiometric_Ca( ref.tsG, ref.tsR, 11 );
ref.spikes = nndORoasis(ref.R, 2, 0.94, 2.4);


%% Register other images and calculate tsG, tsR, R, spikes

% fam1 day2
for n = 1:size(files_fam1d2,1)
    str = sprintf('Processing %s ...\n', files_fam1d2(n,:));
    cprintf(str)
%     load([filedir1 files_fam1d2(n,:) '/' files_fam1d2(n,:) '_2P_mcorr_output.mat'])
%     fam1d2(n).green = green;
%     fam1d2(n).red = red;
%     clear green red 

    fam1d2_imG = read_file( [filedir1 files_fam1d2(n,:) '/' files_fam1d2(n,:) '_2P_XYT_green_mcorr.tif'] );
    fam1d2_imR = read_file( [filedir1 files_fam1d2(n,:) '/' files_fam1d2(n,:) '_2P_XYT_red_mcorr.tif'] );

    [ fam1d2(n).shift, fam1d2(n).red.globalreg_meanframe ] = globalregisterImage(ref.red.meanregframe, fam1d2(n).red.meanregframe, 1 );
    [ fam1d2_imGglobalreg, fam1d2_imRglobalreg ] = ...
        globalregisterStack( fam1d2_imG, fam1d2_imR, fam1d2(n).shift(1) , fam1d2(n).shift(2) );
    clear fam1d2_imG fam1d2_imR
    
    fname_tif_fam1d2_imGglobalreg = [filedir1 files_fam1d2(n,:) '/' files_fam1d2(n,:) '_2P_XYT_green_mcorr_globalreg.tif'];
        writeTifStack( fam1d2_imGglobalreg, fname_tif_fam1d2_imGglobalreg );
    fname_tif_fam1d2_imRglobalreg = [filedir1 files_fam1d2(n,:) '/' files_fam1d2(n,:) '_2P_XYT_red_mcorr_globalreg.tif'];
        writeTifStack( fam1d2_imRglobalreg, fname_tif_fam1d2_imRglobalreg );
    clear fname_tif_fam1d2_imGglobalreg fname_tif_fam1d2_imRglobalreg
        
    for i = 1:Nmasks
        maskind = masks(:,:,i);
        for j = 1:size(fam1d2_imGglobalreg,3)
            imG_reshaped = reshape( fam1d2_imGglobalreg(:,:,j), 512*512, 1);
            fam1d2(n).tsG( i, j ) = mean( imG_reshaped(maskind) );
            imR_reshaped = reshape( fam1d2_imRglobalreg(:,:,j), 512*512, 1);
            fam1d2(n).tsR( i, j ) = mean( imR_reshaped(maskind) );
        end
    end
    clear fam1d2_imGglobalreg fam1d2_imRglobalreg imG_reshaped imR_reshaped maskind i j

    fam1d2(n).R = ratiometric_Ca( fam1d2(n).tsG, fam1d2(n).tsR, 11 );
    fam1d2(n).spikes = nndORoasis(fam1d2(n).R, 2, 0.94, 2.4);
end

% nov day2
for n = 1:size(files_novd2,1)
    str = sprintf('Processing %s ...\n', files_novd2(n,:));
    cprintf(str)
%     load([filedir1 files_novd2(n,:) '/' files_novd2(n,:) '_2P_mcorr_output.mat'])
%     novd2(n).green = green;
%     novd2(n).red = red;
%     clear green red params

    novd2_imG = read_file( [filedir1 files_novd2(n,:) '/' files_novd2(n,:) '_2P_XYT_green_mcorr.tif'] );
    novd2_imR = read_file( [filedir1 files_novd2(n,:) '/' files_novd2(n,:) '_2P_XYT_red_mcorr.tif'] );

    [ novd2(n).shift, novd2(n).red.globalreg_meanframe ] = globalregisterImage(ref.red.meanregframe, novd2(n).red.meanregframe, 1 );
    [ novd2_imGglobalreg, novd2_imRglobalreg ] = ...
        globalregisterStack( novd2_imG, novd2_imR, novd2(n).shift(1), novd2(n).shift(2) );
    clear novd2_imG novd2_imR
    
    fname_tif_novd2_imGglobalreg = [filedir1 files_novd2(n,:) '/' files_novd2(n,:) '_2P_XYT_green_mcorr_globalreg.tif'];
        writeTifStack( novd2_imGglobalreg, fname_tif_novd2_imGglobalreg );
    fname_tif_novd2_imRglobalreg = [filedir1 files_novd2(n,:) '/' files_novd2(n,:) '_2P_XYT_red_mcorr_globalreg.tif'];
        writeTifStack( novd2_imRglobalreg, fname_tif_novd2_imRglobalreg );
    clear fname_tif_novd2_imGglobalreg fname_tif_novd2_imRglobalreg
        
    for i = 1:Nmasks
        maskind = masks(:,:,i);
        for j = 1:size(novd2_imGglobalreg,3)
            imG_reshaped = reshape( novd2_imGglobalreg(:,:,j), 512*512, 1);
            novd2(n).tsG( i, j ) = mean( imG_reshaped(maskind) );
            imR_reshaped = reshape( novd2_imRglobalreg(:,:,j), 512*512, 1);
            novd2(n).tsR( i, j ) = mean( imR_reshaped(maskind) );
        end
    end
    clear novd2_imGglobalreg novd2_imRglobalreg imG_reshaped imR_reshaped maskind i j

    novd2(n).R = ratiometric_Ca( novd2(n).tsG, novd2(n).tsR, 11 );
    novd2(n).spikes = nndORoasis(novd2(n).R, 2, 0.94, 2.4);
end

% fam1 day3
for n = 1:size(files_fam1d3,1)
    str = sprintf('Processing %s ...\n', files_fam1d3(n,:));
    cprintf(str)
%     load([filedir2 files_fam1d3(n,:) '/' files_fam1d3(n,:) '_2P_mcorr_output.mat'])
%     fam1d3(n).green = green;
%     fam1d3(n).red = red;
%     clear green red 

    fam1d3_imG = read_file( [filedir2 files_fam1d3(n,:) '/' files_fam1d3(n,:) '_2P_XYT_green_mcorr.tif'] );
    fam1d3_imR = read_file( [filedir2 files_fam1d3(n,:) '/' files_fam1d3(n,:) '_2P_XYT_red_mcorr.tif'] );

    [ fam1d3(n).shift, fam1d3(n).red.globalreg_meanframe ] = globalregisterImage(ref.red.meanregframe, fam1d3(n).red.meanregframe, 1 );
    [ fam1d3_imGglobalreg, fam1d3_imRglobalreg ] = ...
        globalregisterStack( fam1d3_imG, fam1d3_imR, fam1d3(n).shift(1) , fam1d3(n).shift(2) );
    clear fam1d3_imG fam1d3_imR
    
    fname_tif_fam1d3_imGglobalreg = [filedir2 files_fam1d3(n,:) '/' files_fam1d3(n,:) '_2P_XYT_green_mcorr_globalreg.tif'];
        writeTifStack( fam1d3_imGglobalreg, fname_tif_fam1d3_imGglobalreg );
    fname_tif_fam1d3_imRglobalreg = [filedir2 files_fam1d3(n,:) '/' files_fam1d3(n,:) '_2P_XYT_red_mcorr_globalreg.tif'];
        writeTifStack( fam1d3_imRglobalreg, fname_tif_fam1d3_imRglobalreg );
    clear fname_tif_fam1d3_imGglobalreg fname_tif_fam1d3_imRglobalreg
        
    for i = 1:Nmasks
        maskind = masks(:,:,i);
        for j = 1:size(fam1d3_imGglobalreg,3)
            imG_reshaped = reshape( fam1d3_imGglobalreg(:,:,j), 512*512, 1);
            fam1d3(n).tsG( i, j ) = mean( imG_reshaped(maskind) );
            imR_reshaped = reshape( fam1d3_imRglobalreg(:,:,j), 512*512, 1);
            fam1d3(n).tsR( i, j ) = mean( imR_reshaped(maskind) );
        end
    end
    clear fam1d3_imGglobalreg fam1d3_imRglobalreg imG_reshaped imR_reshaped maskind i j

    fam1d3(n).R = ratiometric_Ca( fam1d3(n).tsG, fam1d3(n).tsR, 11 );
    fam1d3(n).spikes = nndORoasis(fam1d3(n).R, 2, 0.94, 2.4);
end

% fam1rev day3
for n = 1:size(files_fam1revd3,1)
    str = sprintf('Processing %s ...\n', files_fam1revd3(n,:));
    cprintf(str)
%     load([filedir2 files_fam1revd3(n,:) '/' files_fam1revd3(n,:) '_2P_mcorr_output.mat'])
%     fam1revd3(n).green = green;
%     fam1revd3(n).red = red;
%     clear green red 

    fam1revd3_imG = read_file( [filedir2 files_fam1revd3(n,:) '/' files_fam1revd3(n,:) '_2P_XYT_green_mcorr.tif'] );
    fam1revd3_imR = read_file( [filedir2 files_fam1revd3(n,:) '/' files_fam1revd3(n,:) '_2P_XYT_red_mcorr.tif'] );

    [ fam1revd3(n).shift, fam1revd3(n).red.globalreg_meanframe ] = globalregisterImage(ref.red.meanregframe, fam1revd3(n).red.meanregframe, 1 );
    [ fam1revd3_imGglobalreg, fam1revd3_imRglobalreg ] = ...
        globalregisterStack( fam1revd3_imG, fam1revd3_imR, fam1revd3(n).shift(1) , fam1revd3(n).shift(2) );
    clear fam1revd3_imG fam1revd3_imR
    
    fname_tif_fam1revd3_imGglobalreg = [filedir2 files_fam1revd3(n,:) '/' files_fam1revd3(n,:) '_2P_XYT_green_mcorr_globalreg.tif'];
        writeTifStack( fam1revd3_imGglobalreg, fname_tif_fam1revd3_imGglobalreg );
    fname_tif_fam1revd3_imRglobalreg = [filedir2 files_fam1revd3(n,:) '/' files_fam1revd3(n,:) '_2P_XYT_red_mcorr_globalreg.tif'];
        writeTifStack( fam1revd3_imRglobalreg, fname_tif_fam1revd3_imRglobalreg );
    clear fname_tif_fam1revd3_imGglobalreg fname_tif_fam1revd3_imRglobalreg
        
    for i = 1:Nmasks
        maskind = masks(:,:,i);
        for j = 1:size(fam1revd3_imGglobalreg,3)
            imG_reshaped = reshape( fam1revd3_imGglobalreg(:,:,j), 512*512, 1);
            fam1revd3(n).tsG( i, j ) = mean( imG_reshaped(maskind) );
            imR_reshaped = reshape( fam1revd3_imRglobalreg(:,:,j), 512*512, 1);
            fam1revd3(n).tsR( i, j ) = mean( imR_reshaped(maskind) );
        end
    end
    clear fam1revd3_imGglobalreg fam1revd3_imRglobalreg imG_reshaped imR_reshaped maskind i j

    fam1revd3(n).R = ratiometric_Ca( fam1revd3(n).tsG, fam1revd3(n).tsR, 11 );
    fam1revd3(n).spikes = nndORoasis(fam1revd3(n).R, 2, 0.94, 2.4);
end
clear Nmasks data_locn filedir1 filedir2 filedir3 files_fam1d2 files_novd2 files_fam1d3 files_fam1revd3
clear file_ref str n i j mean_imratio


%% Save output
fam1d2_0946 = ref;

fam1d2_0933 = fam1d2(1);
fam1d2_0937 = fam1d2(2);
fam1d2_0951 = fam1d2(3);
clear fam1d2

novd2_0959 = novd2(1);
novd2_1004 = novd2(2);
clear novd2

fam1d3_0922 = fam1d3(1);
fam1d3_0944 = fam1d3(2);
fam1d3_0949 = fam1d3(3);
clear fam1d3

fam1revd3_0957 = fam1revd3(1);
fam1revd3_1007 = fam1revd3(2);
clear fam1revd3

save('m62_fam1_nov_fam1rev_20181015to16_segment.mat')

toc


    
