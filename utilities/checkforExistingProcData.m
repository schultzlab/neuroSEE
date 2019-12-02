% Written by Ann Go
% Function to determine whether file has been processed using same
% parameters as specified by user

function check = checkforExistingProcData(data_locn, text, params)
    mcorr_method = params.methods.mcorr_method;
    segment_method = params.methods.segment_method;
    dofissa = params.methods.dofissa;
    
    if dofissa
        str_fissa = 'FISSA';
    else
        str_fissa = 'noFISSA';
    end

    if strcmp(text(end-2:end),'txt')
        list = text; tlist = true;
    else
        file = text; tlist = false;
    end
    
    if tlist
        [mouseid,expname] = find_mouseIDexpname(list);
        listfile = [data_locn 'Digital_Logbook/lists/' list];
        files = extractFilenamesFromTxtfile(listfile);
        reffile = files(1,:);
        
        dir_proc = [data_locn 'Analysis/' mouseid '/environment_PFmaps/' mouseid '_' expname '_ref_' reffile '/'...
                    mcorr_method '_' segment_method '_' str_fissa '/'];
        
        check = zeros(1,4);
        if exist(dir_proc,'dir')
            % 1) Check for existing registered rois across all image files
            % on the list
            if exist([dir_proc  mouseid '_' expname '_ref_' reffile '_registered_rois.mat'],'file')
                check(1) = 1;
            end
            
            % 2) Check for existing concatenated spiking and tracking data
            if exist([dir_proc  mouseid '_' expname '_ref_' reffile '_spikes_tracking_data.mat'],'file')
                check(2) = 1;
            end
            
            % 3) Check for existing PF mapping output
            if exist([dir_proc  mouseid '_' expname '_ref_' reffile '_PFmap_output.mat'],'file')
                check(3) = 1;
            end
            
            % 4) Check if mat file for all proc data for the file exists
            if exist([dir_proc  mouseid '_' expname '_ref_' reffile '/'...
                    mcorr_method '_' segment_method '_' str_fissa '_allData.mat'],'file')
                check(4) = 1;
            end
        end
    else 
        dir_proc = [data_locn 'Data/' file(1:8) '/Processed/' file '/'];
        
        check = zeros(1,7);
        if exist(dir_proc,'dir')
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
            if exist([dir_fissa file '_spikes_output.mat'],'file')
                check(4) = 1;
            end

            % 5) Check for tracking data
            dir_track = [dir_proc 'behaviour/'];
            matfiles = dir(fullfile(dir_track,['*.','mat']));
            if numel(matfiles) > 0
                check(5) = 1;
            end

            % 6) Check for existing PF mapping output
            if exist([dir_fissa 'PFmaps/' file '_PFmap_output.mat'],'file')
                check(6) = 1;
            end

            % 7) Check if mat file for all proc data for the file exists
            if exist([dir_proc file '_' mcorr_method '_' segment_method '_' str_fissa '_allData.mat'],'file')
                check(6) = 1;
            end
        end
    end
    
end