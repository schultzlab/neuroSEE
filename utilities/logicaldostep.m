function dostep = logicaldostep(dostep)

% These steps are included in single file processing (image registration pipeline)
% dostep(1): dostep motion correction even if motion corrected images exist
% dostep(2): dostep roi segmentation
% dostep(3): dostep neuropil decontamination
% dostep(4): dostep spike extraction
% dostep(5): dostep tracking data extraction
% dostep(6): dostep place field mapping

% These steps are included in multiple file processing (image registration pipeline)
% dostep(1): dostep image registration of individual files
% dostep(2): dostep group roi segmentation
% dostep(3): dostep group neuropil decontamination
% dostep(4): dostep group spike extraction
% dostep(5): dostep group tracking data extraction
% dostep(6): dostep group place field mapping

% Allow only combinations of values that make sense
if dostep(6)
    dostep(1:5) = true; 
elseif dostep(5)
    dostep(1:4) = true; 
elseif dostep(4)
    dostep(1:3) = true; 
elseif dostep(3)
    dostep([1,2]) = true; 
elseif dostep(2)
    dostep(1) = true;
end

