clear; close all;
tic

%% Basic setup
test = 1;       % set to 1 if testing, this will use one of smaller files in ../test
default = 0;    % set to 1 to use default parameters
display = 0;    % set to 1 to display results (all results are saved in 
                %    individual file directories with summary pdfs regardless)
force = 0;      % set to 1 to overwrite saved processed files. This will 
                %    force pipeline to redo all steps incl. raw to tif
                %    conversion. If force = 0, processing step is skipped
                %    if output of said step already exists in individual
                %    file directory.

% gcp;           % start parallel pool
addpath(genpath('../behaviour'));
addpath(genpath('../intervideo_processing'));
addpath(genpath('../motion_correction'));
addpath(genpath('../PF_mapping'));
addpath(genpath('../ROI_segmentation'));
addpath(genpath('../spike_extraction'));
addpath(genpath('../utilities'));
addpath(genpath('../pipelines'));

% Data location
data_locn = '/Volumes/thefarm2/live/CrazyEights/AD_2PCa/';
if ~exist(data_locn,'dir')
    data_locn = '/Volumes/RDS/project/thefarm2/live/CrazyEights/AD_2PCa/';
end
if ~exist(data_locn,'dir')
    data_locn = '/rds/general/user/mc6017/projects/thefarm2/live/CrazyEights/AD_2PCa/';
end

if test
%     data_locn = '/Users/mc6017/Documents';
     data_locn = [data_locn 'Data_forTesting/'];
end

    
%% Read files


%  file ='20181016_09_22_45'; 


% file = '20181016_09_09_43';

%% Set parameters

if default
    params = load('default_params.mat');
end

if ~default
    % Specify own parameters
    % ROI segmentation
        params.cellrad = 8;            % expected radius of a cell (pixels)
        params.maxcells = 200;          % estimated number of cells in FOV
end

if test
    params.maxcells = 60;  
end
tic;
%% Load motion corrected tif files 
file = '20181016_09_09_43';
[imG, imR] = load_imagefile(data_locn,file,force,'_mcorr');

imG = imG(:,:,1:5:size(imG,3));
imR = imR(:,:,1:5:size(imR,3));

%% Load motion corrected tif files 
file = '20181016_09_14_03';

[imG_2, imR_2]  = load_imagefile(data_locn,file,force,'_mcorr');

imG_2 = imG_2(:,:,1:5:size(imG_2,3));
imR_2 = imR_2(:,:,1:5:size(imR_2,3));

imG = cat(3,imG,imG_2);
imR = cat(3,imR,imR_2);
mean_r = mean(imR,3);
clear imG_2; clear imR_2;


%% Load motion corrected tif files 
file = '20181016_09_18_26';
[imG_3, imR_3] = load_imagefile(data_locn,file,force,'_mcorr');

imG_3 = imG_3(:,:,1:5:size(imG_3,3));
imR_3 = imR_3(:,:,1:5:size(imR_3,3));

imG = cat(3,imG,imG_3);
imR = cat(3,imR,imR_3);
clear imG_3; clear imR_3;
%%
file = '20181016_09_22_45';
[imG_3, imR_3] = load_imagefile(data_locn,file,force,'_mcorr');

imG_3 = imG_3(:,:,1:5:size(imG_3,3));
imR_3 = imR_3(:,:,1:5:size(imR_3,3));

imG = cat(3,imG,imG_3);
imR = cat(3,imR,imR_3);
clear imG_3; clear imR_3;
%%
file = '20181016_09_44_06';
[imG_3, imR_3] = load_imagefile(data_locn,file,force,'_mcorr');

imG_3 = imG_3(:,:,1:5:size(imG_3,3));
imR_3 = imR_3(:,:,1:5:size(imR_3,3));

imG = cat(3,imG,imG_3);
imR = cat(3,imR,imR_3);
clear imG_3; clear imR_3;
%%
file = '20181016_09_49_06';
[imG_3, imR_3] = load_imagefile(data_locn,file,force,'_mcorr');

imG_3 = imG_3(:,:,1:5:size(imG_3,3));
imR_3 = imR_3(:,:,1:5:size(imR_3,3));

imG = cat(3,imG,imG_3);
imR = cat(3,imR,imR_3);
clear imG_3; clear imR_3;
%%
file = '20181016_09_57_43';
[imG_3, imR_3] = load_imagefile(data_locn,file,force,'_mcorr');


imG_3 = imG_3(:,:,1:5:size(imG_3,3));
imR_3 = imR_3(:,:,1:5:size(imR_3,3));

imG = cat(3,imG,imG_3);
imR = cat(3,imR,imR_3);
clear imG_3; clear imR_3;
%%
file = '20181016_10_07_07';
[imG_3, imR_3] = load_imagefile(data_locn,file,force,'_mcorr');

imG_3 = imG_3(:,:,1:5:size(imG_3,3));
imR_3 = imR_3(:,:,1:5:size(imR_3,3));

imG = cat(3,imG,imG_3);
imR = cat(3,imR,imR_3);
clear imG_3; clear imR_3;
%%
file = '20181016_10_11_35';
[imG_3, imR_3] = load_imagefile(data_locn,file,force,'_mcorr');

imG_3 = imG_3(:,:,1:5:size(imG_3,3));
imR_3 = imR_3(:,:,1:5:size(imR_3,3));


imG = cat(3,imG,imG_3);
imR = cat(3,imR,imR_3);
clear imG_3; clear imR_3;
%%
file = '20181016_10_15_57';
[imG_3, imR_3] = load_imagefile(data_locn,file,force,'_mcorr');

imG_3 = imG_3(:,:,1:5:size(imG_3,3));
imR_3 = imR_3(:,:,1:5:size(imR_3,3));


imG = cat(3,imG,imG_3);
imR = cat(3,imR,imR_3);

clear imG_3; clear imR_3;
t = toc;
str = sprintf('All 10 tiff files loaded in %2.2f minutes',t/60);
SendSlackNotification('https://hooks.slack.com/services/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',str, '#general','@manfredi.castelli17', [], [], []);
%% Use ABLE to extract ROIs and raw time series
% Saved: image with ROIs (fig, pdf), mat with fields {tsG, tsR, masks, mean_imratio, params}

% [tsG, tsR, masks, mean_imratio, params] = neuroSEE_segment_Manfredi( imG(:,:,1:5:size(imG,3)), imR(:,:,1:5:size(imG,3)), mean(imR,3), ...
%                                                              data_locn, file, params, force );
[tsG, tsR, masks, mean_imratio, params] = neuroSEE_segment_Manfredi( imG, imR, mean(imR,3), ...
                                                             data_locn, file, params, force );
[R, spikes, params] = neuroSEE_extractSpikes( tsG, tsR, data_locn, file, params, force );
GUI_viewROIsSpikes( mean_imratio, masks, tsG, tsR, spikes );
    


t = toc;
str = sprintf('%s: Processing done in %g hrs', file, round(t/3600,2));
cprintf(str)