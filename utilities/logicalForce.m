function force = logicalForce(force)

% force(1): force motion correction even if motion corrected images exist
% force(2): force roi segmentation
% force(3): force neuropil decontamination
% force(4): force spike extraction
% force(5): force tracking data extraction
% force(6): force roi registration across images
% force(7): force spike/tracking data consolidation across images
% force(8): force place field mapping

% Allow only combinations of values that make sense
if force(1)
    force([2:4,6:8]) = true; % because redoing motion correction step affects all succeeding steps except step 5
elseif force(2)
    force([3:4,6:8]) = true; 
elseif force(3)
    force([4,6:8]) = true; 
elseif force(4)
    force(6:8) = true;
elseif force(5)
    force(7:8) = true;
elseif force(6)
    force(7:8) = true;
elseif force(7)
    force(8) = true;
end
