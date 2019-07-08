function [segmented_image,D,N,B,stat] = segment_Manfredi(img,tune)

% Author:           Manfredi Castelli
% Date:             25/06/2019
% Supervisors:      Simon R Schultz

% Overview:         This function contains the segmentation algorithm trough Watershed trasform.
% Acknowledgment:   Stephanie Reynolds, Pier Luigi Dragotti
% Reference:        Reynolds et al. (2016) ABLE: an activity-based level set 
%                   segmentation algorithm for two-photon calcium imaging data
%
% 
%%%%%%%%%%%%%%%   INPUTS    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% img:            Image in form of MxN array to be segmented. 
% tune:           Filtering option: 
%                 -1: red channel;
%                 -2 ovelay of red and green channel;
%                 -3: correlated green channel image;
% 
%
%%%%%%%%%%%%%%%   OUTPUTS    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Segmented_image:  segmented image of input img;
% D:                distance transform output of img;
% N,B:              boundaries stats;
% stat:             measurement and features of img and its objects.
% 


%1) convert mat to grayscale image->binarize img->remove noise
img = mat2gray(img);
img = imadjust(img);
%  cut is the mininum number of pixels of each cell 
    if tune == 1 %red channel image
        cut = 100;
        
    elseif tune == 2 % ovelray
            cut = 60;
            
    else 
%         correlated green image has lots of noise around border.
        cut = 140;
        img = imclearborder(img,4);
        
    end

         
%     img = rgb2gray(img);
%     img = medfilt2(img);
    
% binarizing image with adaptive threshold
    img = imbinarize(img, 'adaptive');
    img = bwareaopen(img, cut);
    img = imfill(img,'holes');


% bwdist uses a fast algorithm to compute the true Euclidean distanc
% transform. it is the distance at each pixel to the nearest nonzero pixel.
    D = -bwdist(~img);
    
    mask = imextendedmin(D,2);

    % applying watershed and segmenting image
    D2 = imimposemin(D,mask);
    Ld2 = watershed(D2);
    segmented_image = img;
    segmented_image(Ld2 == 0) = 0;
    
%     filtering noisy segmented image for well defined boundaries
    segmented_image = medfilt2(segmented_image,[3,3]);
    segmented_image = imclearborder(segmented_image,4);
    
     if tune == 3
        segmented_image = bwareaopen(segmented_image, 130);
     
     else
         segmented_image = bwareaopen(segmented_image, 70);
     end
         

    [B,~,N,~] = bwboundaries(img,'noholes');

%     labelling each object in image and extracting features
    Ilabel = bwlabel(segmented_image);
    stat = regionprops(Ilabel,  'all');
end

