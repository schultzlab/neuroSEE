%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   NeuroSEE: An automated Neuronal Source Extraction
%             and Exploration toolbox
%   
%   Author: Seigfred Prado   
%   Supervisor: Simon Schultz
%   Acknowledgment: Stephanie Reynolds, Pier Luigi Dragotti
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [cell_tsG, cell_tsR] = pipeline_neuroSEE_batch(filename,sumImage,cellRadius,numCells)

if nargin<5; numCells=150; end
if nargin<4; cellRadius=10; end

disp('started')
%%%%%%%%% Load data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
filename1 = filename;

% Get the information of the imaging data file
info = imfinfo(filename1,'tiff');
num_frames = length(info); 
% Initialise the arrays for the green and red channels
green_channel = zeros(info(1).Height,info(1).Width,num_frames/2,'single');
red_channel = zeros(info(1).Height,info(1).Width,num_frames/2,'single');

% Load the images and save the values in the array
for i = 1:2:num_frames
    green_channel(:,:,(i+1)/2) = imread(filename1,'tiff',i);
    red_channel(:,:,(i+1)/2) = imread(filename1,'tiff',i+1);
end
% Make sure there are no zeros
green_channel(find(green_channel<1)) = 1;
red_channel(find(red_channel<1)) = 1;

dim = size(green_channel); % Dimensions of each channel
N = dim(3); % Number of frames per channel

% Add the folder where the other scripts are located
addpath([pwd, 'segmentation/dependencies/'])
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%% Compute the summary images %%%%%%%%%%%%%%%%%%%%%%
%
disp('sumImage in')
filename2 = sumImage;

meanIm_red = mean(red_channel,3);
meanIm_green = mean(green_channel,3);
corrIm_red = crossCorr(red_channel);
corrIm_green = crossCorr(green_channel);
% The summary image can be from the mean or correlation above or 
% a composite image of red and green channels as below:
summary_image = imread(filename2); 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% Initialisation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Set the initial images
disp('initial images set')
green_initial_image = medfilt2(corrIm_green, [5 5]);
red_initial_image = meanIm_red;
%
% Set alpha parameter

radius = cellRadius; % expected cell radius in pixels; usually set to 10

alpha = 0.05; % starts with a low value then increase that during tuning

init_opt.blur_radius = 4; % default is 1; radius of the blurring applied
                          % to the input summary image
                          
% Secondary Metric
% Divide the red over the green average image so that the only remaining
% signal is from neuronal nuclei, thus eliminating field illumination
% inhomogeneity and non-cellular structures like blood vessels
init_opt.secondary_metric = meanIm_red./meanIm_green;
retuned_opt.secondary_metric = init_opt.secondary_metric;
init_opt.second_alpha = 0.05; % same as alpha

% Initialise by getting the candidate ROI seeds
phi_0 = initialise(green_initial_image, radius, alpha, init_opt);

exp_ROIs = numCells; 

retuned_alpha = alpha;
initial_masks_num = size(phi_0,3);
retuned_masks_num = initial_masks_num;
while (retuned_masks_num > exp_ROIs)
    retuned_alpha = retuned_alpha + 0.05;
    retuned_opt.second_alpha = retuned_alpha;
    disp(['Current Alpha is ', num2str(retuned_alpha)]);
    %{
    if retuned_alpha > 0.95;
        break
    end
    %}
    phi = initialise(green_initial_image, radius, retuned_alpha, retuned_opt);
    retuned_masks_num = size(phi,3);
end
%
phi_0=phi;

%
% Set lambda parameter

% Initialise the level set algorithm parameters
seg_opt.lambda              = 200; % this depends on the data
seg_opt.mergeCorr           = 0.95; % correlation coefficient threshold
                                    % above which two neigbouring ROIs will
                                    % be merged
seg_opt.mergeDuring         = 1;

% Tune lambda based on the data
[lambda, phi] = lambda_tune(phi_0, green_channel, red_channel, radius,...
                                 seg_opt, corrIm_green, meanIm_red);
seg_opt.lambda = lambda;

tic;
%

%%%%%%%%%% Segmentation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
seg_opt.corrIm = corrIm_green;
seg_opt.meanIm = meanIm_red;
%seg_opt = rmfield(seg_opt,'plot_progress');
[masks, cell_tsG, nhbd_tsG, cell_tsR, nhbd_tsR] = segment_cells(phi_0, ...
    green_channel, red_channel, radius, seg_opt);
runtime = toc
mask_num = size(masks,3); % Number of detected ROIs

%{
% Display initial results
clear opts;
% Plot masks on summary image
opts.plot_ids = 1; %Set to 1 if you want to view the ID number of the ROIs
plotContoursOnSummaryImage(summary_image, masks, opts);
%}

filename2 = filename(1:20);
save(filename2, '-v7.3');
clear all;
end

