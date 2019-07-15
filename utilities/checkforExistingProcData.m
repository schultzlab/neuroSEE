% Written by Ann Go
% Function to determine whether file has been processed using same
% parameters as specified by user

function check = checkforExistingProcData(data_locn, file, mcorr_method, segment_method, dofissa)
    if dofissa
        str_fissa = '_FISSA';
    else
        str_fissa = 'noFISSA';
    end

    check = zeros(1,7);
    dir_proc = [data_locn 'Data/' file(1:8) '/Processed/' file '/'];
    % 1) Check for existing motion corrected tif files and mat file in Processed folder
    dir_mcorr = [dir_proc 'mcorr_' mcorr_method '/'];
    if all( [exist([dir_mcorr file '_2P_XYT_green_mcorr.tif'],'file'),...
             exist([dir_mcorr file '_2P_XYT_red_mcorr.tif'],'file'),...
             exist([dir_mcorr file '_mcorr_output.mat'],'file')] )
        check(1) = 1;
    end
    
    % 2) Check for existing roi segmentation output
    dir_segment = [dir_mcorr segment_method '/'];
    if exist([dir_segment file '_segment_output.mat'],'file')
        check(2) = 1;
    end
    
    % 3) Check for existing neuropil decontamination output
    dir_fissa = [dir_segment 'FISSA/'];
    if exist([dir_fissa file '_fissa_output.mat'],'file')
        check(3) = 1;
    end
    
    % 4) Check for existing spike estimation output
    if exist([dir_fissa file '_spike_output.mat'],'file')
        check(4) = 1;
    end

    % 5) Check for tracking data
    dir_track = [dir_proc 'behaviour/'];
    matfiles = dir(fullfile(dir_track,['*.','mat']));
    if numel(matfiles) > 0
        check(5) = 1;
    end
    
    % 6) Check for existing PF mapping output
    if exist([dir_fissa file '_PFmapping_output.mat'],'file')
        check(6) = 1;
    end

    % 7) Check if mat file for all proc data for the file exists
    matfiles = dir(fullfile(dir_proc,['*.','mat']));
    if numel(matfiles) > 0 
        for i = 1:numel(matfiles)
            name = matfiles(i).name;
            if contains(name,'allData')
                if contains(name,mcorr_method)
                    if contains(name,segment_method)
                        if contains(name,str_fissa)
                            check(7) = 1;
                        end
                    end
                end
            end
        end
    end
end