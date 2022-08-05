% Written by Ann Go
%
% This function adds module folders required in the pipelines and returns
% the data folder location and any errors

function [data_locn,comp,err] = load_neuroSEEmodules(test)

if nargin < 1, test = false; end

%% Add module folders
addpath(genpath('../behaviour'));
addpath(genpath('../motion_correction'));
addpath(genpath('../neuropil_decontamination'));
addpath(genpath('../PF_mapping'));
addpath(genpath('../pipelines'));
addpath(genpath('../quality_control'));
addpath(genpath('../ROI_registration'));
addpath(genpath('../ROI_segmentation'));
addpath(genpath('../spike_extraction'));
addpath(genpath('../utilities'));

err = [];

%% Data location
mydir  = pwd;
ind   = strfind(mydir,'/');
user = mydir(ind(2)+1:ind(3)-1);

% Works in local Macs
data_locn = '/Volumes/thefarm2/live/CrazyEights/AD_2PCa/';
comp = 'mac';
if ~exist(data_locn,'dir')
    data_locn = '/Volumes/RDS/project/thefarm2/live/CrazyEights/AD_2PCa/';
end
if ~exist(data_locn,'dir')
    data_locn = ['/Volumes/' user '/projects/thefarm2/live/CrazyEights/AD_2PCa/'];
end

% Works in linuxbox
if ~exist(data_locn,'dir')
    data_locn = ['/rds/general/user/' user '/projects/thefarm2/live/CrazyEights/AD_2PCa/'];
    comp = 'linuxbox';
end
% Works in HPC
if ~exist(data_locn,'dir')
    data_locn = '/rds/general/project/thefarm2/live/CrazyEights/AD_2PCa/';
    comp = 'hpc';
end

% If data_locn can't be found, thefarm2 might not be mounted
if ~exist(data_locn,'dir')
    data_locn = []; comp = [];
    err = '"data_locn" not found. "thefarm2" may not be mounted.\n';
end
if isempty(err)    
    if test
        data_locn = [data_locn 'Data_forTesting/'];
    end
end
