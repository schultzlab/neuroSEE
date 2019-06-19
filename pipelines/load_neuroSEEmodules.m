function data_locn = load_neuroSEEmodules(test)

% Add module folders
addpath(genpath('../behaviour'));
addpath(genpath('../motion_correction'));
addpath(genpath('../neuropil_decontamination'));
addpath(genpath('../PF_mapping'));
addpath(genpath('../ROI_registration'));
addpath(genpath('../ROI_segmentation'));
addpath(genpath('../spike_extraction'));
addpath(genpath('../utilities'));
addpath(genpath('../pipelines'));

% Data location
% Works in Ann's Macs
data_locn = '/Volumes/thefarm2/live/CrazyEights/AD_2PCa/';
if ~exist(data_locn,'dir')
    data_locn = '/Volumes/RDS/project/thefarm2/live/CrazyEights/AD_2PCa/';
end
% Works in linuxbox
if ~exist(data_locn,'dir')
    data_locn = '/rds/general/user/mgo/projects/thefarm2/live/CrazyEights/AD_2PCa/';
end
% Works in HPC
if ~exist(data_locn,'dir')
    data_locn = '/rds/general/project/thefarm2/live/CrazyEights/AD_2PCa/';
end

if test
    data_locn = [data_locn 'Data_forTesting/'];
end
