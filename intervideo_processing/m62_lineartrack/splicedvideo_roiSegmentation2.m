clear; close all;
tic

%% Basic setup
addpath(genpath('../behaviour'));
addpath(genpath('../intervideo_processing'));
addpath(genpath('../motion_correction'));
addpath(genpath('../PF_mapping'));
addpath(genpath('../ROI_segmentation'));
addpath(genpath('../spike_extraction'));
addpath(genpath('../utilities'));
addpath(genpath('../pipelines'));

% Data location
data_locn = '/Volumes/thefarm2/live/CrazyEights/AD_2PCa/Summaries/m62/fam1_nov_fam1rev/';
if ~exist(data_locn,'dir')
    data_locn = '/Volumes/RDS/project/thefarm2/live/CrazyEights/AD_2PCa/Summaries/m62/fam1_nov_fam1rev/';
end
if ~exist(data_locn,'dir')
    data_locn = '/rds/general/user/mgo/projects/thefarm2/live/CrazyEights/AD_2PCa/Summaries/m62/fam1_nov_fam1rev/';
end

%% file
file = 'first2000SplicedStacks';
file_gr = [data_locn file '_gr.tif'];
file_red = [data_locn file '_red.tif'];
str = sprintf('Reading %s _gr.tif\n', file);
cprintf(str)
imG = read_file( file_gr );
str = sprintf('Reading %s _red.tif\n', file);
cprintf(str)
imR = read_file( file_red );

params = load('default_params.mat');

% Use ABLE to extract ROIs and raw time series
cellrad = params.cellrad;
maxcells = params.maxcells;

cell_tsR = zeros(1,size(imR,3));
[cell_tsG, masks, mean_imratio] = ABLE( imG, mean(imR,3), file, cellrad, maxcells );
for i = 1:size(masks,3)
    ind = find(masks(:,:,i));
    for j = 1:size(imR,3)
        imR_reshaped = reshape(imR(:,:,j),size(imR,1)*size(imR,2),1);
        cell_tsR(i,j) = mean(imR_reshaped(ind));
    end
end

% Save output
% Plot masks on summary image and save plot
plotopts.plot_ids = 1; % set to 1 to view the ID number of the ROIs on the plot
fig = plotContoursOnSummaryImage(mean_imratio, masks, plotopts);
fname_fig = [data_locn file];
savefig(fig,fname_fig);
saveas(fig,fname_fig,'pdf');
close(fig);

% Save masks in a mat file
fname_out = [data_locn file '_segment_output.mat'];
str = sprintf('Saving %s ...\n', fname_out);
cprintf(str);
save(fname_out,'cell_tsG','cell_tsR','masks','mean_imratio','params')

toc