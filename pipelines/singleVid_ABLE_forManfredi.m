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
    data_locn = '/rds/general/user/mgo/projects/thefarm2/live/CrazyEights/AD_2PCa/';
end

if test
    data_locn = [data_locn 'Data_forTesting/'];
end

    
%% Read files

file = '20181016_10_11_35'; 
% file = '20190406_20_38_41';

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

%% Load motion corrected tif files 

[imG, imR] = load_imagefile(data_locn,file,force,'_mcorr');

%% Use ABLE to extract ROIs and raw time series
% Saved: image with ROIs (fig, pdf), mat with fields {tsG, tsR, masks, mean_imratio, params}

[tsG, tsR, masks, mean_imratio, params] = neuroSEE_segment( imG, imR, mean(mean(imR,3)), ...
                                                            data_locn, file, params, force );


t = toc;
str = sprintf('%s: Processing done in %g hrs', file, round(t/3600,2));
cprintf(str)