function[masks] = initialiseROISegmentation(metric,red_metric,radius, alpha,...
                                            options,tune_red,tune_green)
%
% Author:      Stephanie Reynolds
% Date:        25/09/2017
% Supervisors: Pier Luigi Dragotti, Simon R Schultz
% Overview:    This function is used as the initialisation for the segmentation
%              algorithm. Peaks in the 2D summary image(s) are identified 
%              as candidate ROIs. Peaks are found with an built-in MATLAB 
%              image processing function 'imextendedmax'. Peaks with 
%              relative height (with respect to their neighbours) that are 
%              less than alpha x sigma are suppressed. Here, sigma is
%              the standard deviation of the summary image and alpha is a 
%              tuning parameter. 
% Reference:   Reynolds et al. (2016) ABLE: an activity-based level set 
%              segmentation algorithm for two-photon calcium imaging data
%
%Acknowledgment:   Stephanie Reynolds, Pier Luigi Dragotti
%%%%%%%%%%%%%%%   INPUTS    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% metric                     MxN summary image of video, usually the
%                            pixelwise cross-correlation (see crossCorr.m)
% radius                     radius of a cell
% alpha                      tuning parameter, peaks below alpha*sigma will be
%                            suppressed
% options                    A variable of type struct. In the following we
%                            describe the fields of this struct. 
% options.blur_radius        [Default: 1] If present, this is the radius of
%                            blurring applied to the summary images. 
%                            blurred with radius options.blur_radius.
% options.secondary_metric   M x N array corresponding to a summary image,
%                            e.g. the mean image. If this field is present, 
%                            initialisation is  performed on both the first 
%                            argument 'metric' and a second summary image. The value
%                            options.second_alpha (the value of alpha for the
%                            second summary image) must also be present. 
% 
%  (Edit by Manfredi)
%  tune                      filtering parameter for segmenting image.
%
%%%%%%%%%%%%%%%   OUTPUTS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% masks                      MxNxK array, where K is the number of ROIs found.
%                            In each sheet (masks(:,:,ii)): -1 indicates that a
%                            pixel is inside the i-th ROI, +1 indicates that the
%                            pixel is outside the i-th ROI.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dim     = size(metric);
maxSize = round(pi * radius^2 * 1.8);

%% Loading pretrained convolutional neural network
load('cnn_2.mat');

%% Segmenting and Making measuremnts of objects detected

[segmentedR, ~,~,~,stats_red] = segment_Manfredi(red_metric,tune_red);

[segmentedG, ~,~,~,stats_green] = segment_Manfredi((metric - 0.05),tune_green);

% Deleting ROI detected  at boundaries due to correlation artifacts
boundaryDetect = [];
count = 1;

for x = 1 : size(stats_red,1)
    if (stats_red(x).Centroid(1) < 3) ||(stats_red(x).Centroid(1) > 509)
        
        boundaryDetect(count) = x;
        count = count +1;
    end
end
stats_red(boundaryDetect) = [];

boundaryDetect = [];
count = 1;
for x = 1 : size(stats_red,1)
    if (stats_red(x).Centroid(2) < 15) ||(stats_red(x).Centroid(2) > 500)
        
        boundaryDetect(count) = x;
        count = count +1;
    end
end
stats_red(boundaryDetect) = [];


h = figure('Name','Green');imshow(metric);hold on; label_centroid(stats_green,'g.');label_centroid(stats_red,'r.');hold off;
savefig(h,'DistributionOfCentroids.fig','compact');



%% Classify each ROI found from red channel 
notCell = seedClassification((metric ),stats_red);
% eliminating roi that classifier discarded 
stats_red(notCell) = [];
% %% Classify each ROI found from green channel 
% notGreen = seedClassification(metric,stats_green);
% stats_green(notGreen) = [];


h = figure('Name','Final centroids after classification');imshow((metric - 0.05));hold on; label_centroid([stats_red;stats_green],'g.');hold off;
savefig(h,'Segmented.fig','compact');

% h = figure('Name','Green');imshow(segmentedG);hold on; label_centroid(stats_green,'g.');label_centroid(stats_red,'r.');hold off;
% savefig(h,'Segmented.fig','compact');

[~,red_Roi_Merged] = merge_Manfredi(stats_red,stats_green);
% eliminating roi that merged 
stats_red(red_Roi_Merged) = [];




%% Creating masks 
obj_num    = size(stats_red,1);
masks      = zeros(dim(1), dim(2), obj_num);
masksG     = zeros(dim(1),dim(2),obj_num);

% msave()inSize = min([stats_red.Area]);

for ii  = 1:obj_num
    mask             = zeros(dim);
    mask(stats_red(ii).PixelIdxList(:)) = 1;
    masks(:,:,ii)    = mask;
end

for ii  = 1:size(stats_green,1)
    maskG             = zeros(dim);
    maskG(stats_green(ii).PixelIdxList(:)) = 1;
    masksG(:,:,ii)    = maskG;
end

masks = cat(3,masksG,masks);

%%BInarizing each mask since concatenation as affected morphology o f pixel

for ii  = 1:size(masks,3)
    masks(:,:,ii) = imbinarize(masks(:,:,ii));
%     masks(:,:,ii) = medfilt2(masks(:,:,ii),[3 3]);
end
% % Remove any that are too large 
% nnzs  = squeeze(sum(sum(masks,1),2));
% masks(:,:,nnzs > maxSize)     = [];
% Final ROIs

masks = -1*masks + ~masks;
h = figure('Name','Green');imagesc(mean(masks,3));hold on; label_centroid(stats_green,'g.');label_centroid(stats_red,'r.');hold off;
savefig(h,'masks+centr.fig','compact');
end









