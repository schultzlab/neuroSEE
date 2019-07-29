function[masks] = initialiseROISegmentation_Manfredi(metric,red_metric,ratio,radius, alpha,...
                                            options,tune_red,tune_green)
%
% Author:      Manfredi Castellil Stephanie Reynolds
% Date:        25/09/2017
% Acknowledgment: Stephanie Reynolds
% Supervisors: Pier Luigi Dragotti, Simon R Schultz
% Overview:    This function is used as the initialisation for the segmentation
%              algorithm. The function computes the segmenttion of the correlation image
%              from green channel; and of red channel mean image.
%              To gets more seeds the red channel seeds are compared to the ones found in the mean_imratio
%              and they are merged(mean_imratio -by inspection- helped us indentify cell overlap).
%              Finally green and red ROI candidates  are selected based on geometry features
%              (Area and Eccenctricity)/ The final candidates are used for creating var masks
%
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
%  tune                      filtering parameter for segmenting image.
%
%%%%%%%%%%%%%%%   OUTPUTS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% masks                      MxNxK array, where K is the number of ROIs found.
%                            In each sheet (masks(:,:,ii)): -1 indicates that a
%                            pixel is inside the i-th ROI, +1 indicates that the
%                            pixel is outside the i-th ROI.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dim     = size(metric);
maxSize = round(pi * (radius*1.7)^2 );
minSize = round(pi * (0.55*radius)^2);

%% Loading pretrained convolutional neural network - still trying to improve it (89.8 % accurate)
%load('cnn_2.mat');

%% Segmenting and Making measuremnts of objects detected
[~, ~,stats_red] = segment_Manfredi(red_metric+alpha ,tune_red,alpha);
[~, ~,stats_green] = segment_Manfredi(metric - alpha,tune_green,alpha);
[~, ~,stats_ratio] = segment_Manfredi(imadjust(ratio),tune_red,alpha);

% merging ratio and red seeds
[~,Roi_Merged] = merge_Manfredi(stats_ratio,stats_red);
% eliminating roi that merged
stats_ratio(Roi_Merged) = [];

stats_red = [stats_red;stats_ratio];
clear stat_ratio;

% Deleting ROI detected  at boundaries due to correlation artifacts
boundaryDetect = [];
count = 1;

for x = 1 : size(stats_red,1)
    if (stats_red(x).Centroid(1) < 7) ||(stats_red(x).Centroid(1) > 507)

        boundaryDetect(count) = x;
        count = count +1;
    end
end
stats_red(boundaryDetect) = [];

boundaryDetect = [];
count = 1;
for x = 1 : size(stats_red,1)
    if (stats_red(x).Centroid(2) < 6) ||(stats_red(x).Centroid(2) > 506)

        boundaryDetect(count) = x;
        count = count +1;
    end
end
stats_red(boundaryDetect) = [];

%% Classify each ROI found from red channel
%notCell = seedClassification((metric ),stats_red);
% eliminating roi that classifier discarded

%stats_red(notCell) = [];


%% Selecting ROI based on geometrical features
idx = find([stats_red.Area] > maxSize | [stats_red.Area] < minSize | [stats_red.Eccentricity] > 0.9 );
stats_red(idx) = [];

idx = find([stats_green.Area] > maxSize | [stats_green.Area] < minSize |[stats_green.Eccentricity] > 0.92);
stats_green(idx) = [];

[~,red_Roi_Merged] = merge_Manfredi(stats_red,stats_green);
% eliminating roi that merged
stats_red(red_Roi_Merged) = [];

%% Creating masks
obj_num    = size(stats_red,1);
masks      = zeros(dim(1), dim(2), obj_num);
masksG     = zeros(dim(1),dim(2),size(stats_green,1));


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

%BInarizing each mask since concatenation as affected morphology o f pixel

%for ii  = 1:size(masks,3)
%    masks(:,:,ii) = imbinarize(masks(:,:,ii));
%     masks(:,:,ii) = medfilt2(masks(:,:,ii),[3 3]);
%% % Remove any that are too large
% nnzs  = squeeze(sum(sum(masks,1),2));
% masks(:,:,nnzs > maxSize)     = [];
%% Final ROIs

masks = -1*masks + ~masks;
%h = figure('Name','Green');imagesc(mean(masks,3));hold on; label_centroid(stats_green,'g.');label_centroid(stats_red,'r.');hold off;

end
