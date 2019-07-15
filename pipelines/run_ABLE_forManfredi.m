clear; close all;
tic

%% USER: Set basic settings
                    % Set to
test = 0;           % [0,1] 1: debug mode (this will use one of smaller files in test folder)
default = 1;        % [0,1] 1: default parameters
force = [0;...      % [0,1] 1: load raw images instead of tif, i.e. force raw to tif conversion of images
         0;...      % [0,1] 1: redo motion correction even if motion corrected images exist
         0];        % [0,1] 1: force pipeline to redo all steps after motion correction
                    %       0, processing step is skipped if output of said
                    %          step already exists in individual file directory.
mcorr_method = 1;   % [1,2] 1: CaImAn NoRMCorre method, 2: fft-rigid method (Katie's)
segment_method = 2; % [1,2] 1: ABLE, 2: CaImAn
fissa_yn = 1;       % [0,1] 1: implement FISSA, 0: skip step



%% Load module folders and define data directory

data_locn = load_neuroSEEmodules(test);


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

    
%% Read files


%  file ='20181016_09_22_45'; 


% file = '20181016_09_09_43';


file = '20190406_18_51_09'; 
[imG, imR] = load_imagefile(data_locn,file,force,'_mcorr');
% %% Load motion corrected tif files 
% file = '20181016_09_09_43';
% [imG, imR] = load_imagefile(data_locn,file,force,'_mcorr');
% 
imG = imG(:,:,1:5:size(imG,3));
imR = imR(:,:,1:5:size(imR,3));
% 
%% Load motion corrected tif files 
file = '20190406_18_52_49';

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

file = '20190406_18_59_41';
[imG_3, imR_3] = load_imagefile(data_locn,file,force,'_mcorr');


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

% file = '20181016_09_22_45';
% [imG_3, imR_3] = load_imagefile(data_locn,file,force,'_mcorr');
% 
% imG_3 = imG_3(:,:,1:5:size(imG_3,3));
% imR_3 = imR_3(:,:,1:5:size(imR_3,3));
% 
% imG = cat(3,imG,imG_3);
% imR = cat(3,imR,imR_3);
% clear imG_3; clear imR_3;
% %%
% file = '20181016_09_44_06';
% [imG_3, imR_3] = load_imagefile(data_locn,file,force,'_mcorr');
% 
% imG_3 = imG_3(:,:,1:5:size(imG_3,3));
% imR_3 = imR_3(:,:,1:5:size(imR_3,3));
% 
% imG = cat(3,imG,imG_3);
% imR = cat(3,imR,imR_3);
% clear imG_3; clear imR_3;
% %%
% file = '20181016_09_49_06';
% [imG_3, imR_3] = load_imagefile(data_locn,file,force,'_mcorr');
% 
% imG_3 = imG_3(:,:,1:5:size(imG_3,3));
% imR_3 = imR_3(:,:,1:5:size(imR_3,3));
% 
% imG = cat(3,imG,imG_3);
% imR = cat(3,imR,imR_3);
% clear imG_3; clear imR_3;
% %%
% file = '20181016_09_57_43';
% [imG_3, imR_3] = load_imagefile(data_locn,file,force,'_mcorr');
% 
% 
% imG_3 = imG_3(:,:,1:5:size(imG_3,3));
% imR_3 = imR_3(:,:,1:5:size(imR_3,3));
% 
% imG = cat(3,imG,imG_3);
% imR = cat(3,imR,imR_3);
% clear imG_3; clear imR_3;
%%
% file = '20181016_10_07_07';
% file = '20190404_17_29_21';
% [imG_3, imR_3] = load_imagefile(data_locn,file,force,'_mcorr');
% 
% imG_3 = imG_3(:,:,1:5:size(imG_3,3));
% imR_3 = imR_3(:,:,1:5:size(imR_3,3));
% 
% imG = cat(3,imG,imG_3);
% imR = cat(3,imR,imR_3);
% % clear imG_3; clear imR_3;
% %%
% file = '20181016_10_11_35';
% [imG_3, imR_3] = load_imagefile(data_locn,file,force,'_mcorr');
% 
% imG_3 = imG_3(:,:,1:5:size(imG_3,3));
% imR_3 = imR_3(:,:,1:5:size(imR_3,3));
% 
% 
% imG = cat(3,imG,imG_3);
% imR = cat(3,imR,imR_3);
% % clear imG_3; clear imR_3;
% %%
% file = '20181016_10_15_57';
% [imG_3, imR_3] = load_imagefile(data_locn,file,force,'_mcorr');
% 
% imG_3 = imG_3(:,:,1:5:size(imG_3,3));
% imR_3 = imR_3(:,:,1:5:size(imR_3,3));
% 
% 
% imG = cat(3,imG,imG_3);
% imR = cat(3,imR,imR_3);

% clear imG_3; clear imR_3;
t = toc;
str = sprintf('MAC->All 10 tiff files loaded in %2.2f minutes',t/60);
SendSlackNotification('https://hooks.slack.com/services/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',str, '#general','@manfredi.castelli17', [], [], []);

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