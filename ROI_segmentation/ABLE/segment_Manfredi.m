function [segmented_image,D,stat] = segment_Manfredi(img,tune,alpha)

% Author:           Manfredi Castelli
% Date:             25/06/2019
% Supervisors:      Simon R Schultz

% Overview:         This function contains the segmentation algorithm trough Watershed trasform.
% Acknowledgment:   Stephanie Reynolds, Pier Luigi Dragotti
% Reference:        Reynolds et al. (2016) ABLE: an activity-based level set
%                   segmentation algorithm for two-photon calcium imaging data
%
% Overview:         Function detects seeds in two different ways due to differences in
%                   histogram of red and green images.
%                   Red segmetnation is based on background removal;
%                   Green requires only a binarization with an adaptive threshold.
%                   Then, both images are processed by computing the distance trasform
%                   hence,detecing the minima(and then reinverting image). Finally segmetation is done trough
%                   Watershed  trasform. Objects detected are then measured by function regionprops
%
%%%%%%%%%%%%%%%   INPUTS    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% img:            Image in form of MxN array to be segmented.
% tune:           Filtering option:
%                 -1: red channel;
%                 -2 ovelay of red and green channel;
%                 -3: correlated green channel image;
% alpha :         tuning parameter for distance trasform
%
%%%%%%%%%%%%%%%   OUTPUTS    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Segmented_image:  segmented image of input img;
% D:                distance transform output of img;
% stat:             measurement and features of img and its objects.
%
h  = sqrt(nanvar(img(:)));

%1) convert mat to grayscale image->binarize img->remove noise
 if tune == 1 %red channel image
        cut = 90;

        se = strel('disk',8);
        background = imopen(img,se);
        I2 = img - background; %background removal
        img = imadjust(I2);
        img = imbinarize(img, 'adaptive');

        img = bwareafilt(img,[cut 650]); %selction blobs that match area limits specified
        %img = imclearborder(img);
        img = imfill(img,'holes');

    elseif tune == 2 % ovelray
            cut = 60;

    else
%         correlated green image has lots of noise around border.
        cut = 100;

%         img = medfilt2(img,[3,3]);
    img = imbinarize(img, 'adaptive');
    img = bwareaopen(img, cut);
    img = imclearborder(img,8);
    img = imfill(img,'holes');

    end

% bwdist uses a fast algorithm to compute the true Euclidean distanc
% transform. it is the distance at each pixel to the nearest nonzero pixel.
    D = -bwdist(~img);
    mask = imextendedmin(D,h*alpha);
    D2 = imimposemin(D,mask);

    % applying watershed and segmenting image
    Ld2 = watershed(D2);
    segmented_image = img;
    segmented_image(Ld2 == 0) = 0;

%     filtering noisy segmented image for well defined boundaries

    if tune == 1
        segmented_image = bwareaopen(segmented_image, 110,4);
    else
        segmented_image = bwareaopen(segmented_image, 90,4);
    end

    segmented_image = medfilt2(segmented_image,[2 2]);

    %% smoothing edges
    E = edge(segmented_image,'canny');
    Ed = imdilate(E,strel('disk',7));
    %Filtered image
    Ifilt = imfilter(segmented_image,fspecial('gaussian'));
    %Use Ed as logical index into I to and replace with Ifilt
    segmented_image(Ed) = Ifilt(Ed);
    

%     labelling each object in image and extracting features
    Ilabel = bwlabel(segmented_image);
    stat = regionprops(Ilabel,  'all');
end
