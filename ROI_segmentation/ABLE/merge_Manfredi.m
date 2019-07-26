
function [centroids,roiR_merged] = merge_Manfredi(red_cent,green_centr)

% written by Manfredi  Castelli 06/2019
%%%%%%%%%%%%%%%   INPUTS    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% red_cent                   matrix containg centroids of each seed from red image
%
% green_cent                 matrix containg centroids of each seed from green image
%%%%%%%%%%%%%%%   OUTPUTS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% centroids                 array containg final centroids after merging
% roiR_merged               index of all seeds from 1st input that merged
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    roiR_merged = [];
    num_green = size(green_centr,1);
    num_red = size(red_cent,1);

    centroids = zeros(num_green+num_red,2);
    labelled = zeros(num_green,2);
    green_all_centr = zeros(num_green,2);

   % finding min radius equivalent so that min distance for merging depends on enviromented analysed
    min_radius_red = min(mean([red_cent.MajorAxisLength red_cent.MinorAxisLength],2)/2);
    min_radius_green = min(mean([green_centr.MajorAxisLength green_centr.MinorAxisLength],2)/2);

    min_radius = min(min_radius_red,min_radius_green) * 1.8; % multiplying give us some margin

    for x = 1 : num_red
        centroids(x,1) = red_cent(x).Centroid(1);
        centroids(x,2) = red_cent(x).Centroid(2);

    end

    for x = 1 : num_green
        green_all_centr(x,1) = green_centr(x).Centroid(1);
        green_all_centr(x,2) = green_centr(x).Centroid(2);

    end


    i = 0;
    j = 0;
    count = 0;

    %% Calculating seed - wise distance between the two images;
    for i = 1 : num_red

        a = zeros(num_green,2);
        x = centroids(i,1);
        y = centroids(i,2);

        for j = 1 : num_green

            distance = sqrt((x-green_centr(j).Centroid(1)).^2 + (y-green_centr(j).Centroid(2)).^2);
%             distances(j) = distance;

            a(1,1)  = green_centr(j).Centroid(1);
            a(1,2) = green_centr(j).Centroid(2);

            % check whether centroid has already merged with other seeds
            not_present = ~ismember(a,labelled,'rows');

            if (distance <= min_radius) && not_present(1)   % if seed hasn't merged already and distance is small has to merge
                fprintf('Red ROI %d and green ROI %d merged.\n',i,j);

                count = count +1;
                roiR_merged(count) = i;

                labelled(count,1) = green_centr(j).Centroid(1);
                labelled(count,2) = green_centr(j).Centroid(2);

            end


        end


    end

    %vector operation to mark final centroids
    c = ~ismember(green_all_centr,labelled,'rows');
    c  = c .* green_all_centr;
    centroids = unique([centroids;c],'rows');
    centroids(1,:) =[];

    final_numberG = size(centroids,1)-num_red;


    green_percent = (100*size(green_centr,1))/size(centroids,1);
    red_percent = (100*size(red_cent,1))/size(centroids,1);

%     fprintf('\n\nMAC->%d ROI in total merged.',count);
%     fprintf('\n\nMAC->Final number of ROI: %d\n\n',size(centroids,1));
%     fprintf('Total percentage of centroids detected by green channel: %2.2f%%\n',green_percent);
%     fprintf('Total percentage of centroids detected by red channel: %2.2f%%\n',red_percent);
    
end
