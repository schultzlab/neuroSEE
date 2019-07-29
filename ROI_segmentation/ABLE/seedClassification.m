function indexOfNotCell = seedClassification(img,stats_red)
% written by Manfredi to be used when CNN classifier improves
%
%


    load('cnn_2.mat');

    img = imadjust(img)+ 0.2;

    for i = 1: size(stats_red,1)
        title = sprintf('%d.png',i);
        fullFileName = fullfile('/Users/mc6017/Desktop/UROP/test',title);

         x = stats_red(i).Centroid(1) - 12;
         y = stats_red(i).Centroid(2)  - 12;

        bounding_box = [x y 2*12 2*12];

        crop_section = imcrop(img,bounding_box);
        crop_section = imresize(crop_section,[50 50]);

         %Gray scale normalization
         crop_section = mat2gray(crop_section);
         crop_section = im2double(crop_section);
%          crop_section = uint8(255*mat2gray(crop_section));

        imwrite(crop_section,fullFileName);

    end

    imds = imageDatastore(fullfile('/Users/mc6017/Desktop/UROP/test'));
    [YPred,probs] = classify(net,imds);

    notCell = (YPred ~= 'Cell') ;
    indexOfNotCell = find(notCell);


    for i = 1: size(stats_red,1)
        fprintf('%d: %s\n',i,YPred(i));
    end
end
