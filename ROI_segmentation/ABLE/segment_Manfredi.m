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
 if tune == 1 %red channel image
        cut = 140;
        
        se = strel('disk',8);
        background = imopen(img,se);       
        I2 = img - background;
        img = imadjust(I2); 
       
        img = imbinarize(img, 'adaptive');
        img = bwareafilt(img,[140 600]);
        img = imclearborder(img);
        img = imfill(img,'holes');
        
    elseif tune == 2 % ovelray
            cut = 60;
            
    else 
%         correlated green image has lots of noise around border.
        cut = 140;
       
%         img = medfilt2(img,[3,3]);
    img = imbinarize(img, 'adaptive');
    img = bwareaopen(img, cut);
    img = imclearborder(img);
    img = imfill(img,'holes');
        
    end

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
    
     if tune == 3
         segmented_image = bwareaopen(segmented_image, 120);
     
     else
         segmented_image = bwareaopen(segmented_image, 110);
     end
         

    [B,~,N,~] = bwboundaries(img,'noholes');

%     labelling each object in image and extracting features
    Ilabel = bwlabel(segmented_image);
    stat = regionprops(Ilabel,  'all');
end

