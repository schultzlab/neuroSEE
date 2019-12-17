function force = logicalForce(force)

% These steps are included in single file processing (image registration pipeline)
% force(1): force motion correction even if motion corrected images exist
% force(2): force roi segmentation
% force(3): force neuropil decontamination
% force(4): force spike extraction
% force(5): force tracking data extraction
% force(6): force place field mapping

% These steps are included in multiple file processing (image registration pipeline)
% force(1): force image registration of individual files
% force(2): force group roi segmentation
% force(3): force group neuropil decontamination
% force(4): force group spike extraction
% force(5): force group tracking data extraction
% force(6): force group place field mapping

% Allow only combinations of values that make sense
if force(1)
    force([2:4,6]) = true; 
elseif force(2)
    force([3:4,6]) = true; 
elseif force(3)
    force([4,6]) = true; 
elseif force(4)
    force(6) = true;
elseif force(5)
    force(6) = true;
end

