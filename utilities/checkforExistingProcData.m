% Written by Ann Go
%
% Function to determine whether file has been processed using same
% parameters as specified by user
% This can be used to check for processed data for individual files OR for
% processed data for different experiments of a mouse.

function check = checkforExistingProcData(data_locn, text, params, reffile)
    
    %if nargin<4, see line 39
    if nargin<3 
        params = neuroSEE_setparams;
    end
    mcorr_method = params.methods.mcorr_method;
    segment_method = params.methods.segment_method;
    dofissa = params.methods.dofissa;
    bl_prctile = params.spkExtract.bl_prctile;
    
    if dofissa
        str_fissa = 'FISSA';
    else
        str_fissa = 'noFISSA';
    end
    
    % checking processed data for an individual file or an experiment?
    if strcmp(text(end-2:end),'txt')
        list = text; tlist = true;
    else
        file = text; tlist = false;
    end
    
    % if an experiment (i.e. list of files)
    if tlist
        [mouseid,expname,fov] = find_mouseIDexpname(list);
        if nargin<4 
            listfile = [data_locn 'Digital_Logbook/lists/' list];
            files = extractFilenamesFromTxtfile(listfile);
            reffile = files(1,:); 
        end
        
        groupreg_method = params.methods.groupreg_method;
        if ~isempty(fov)
            dir_proc = [data_locn 'Analysis/' mouseid '/' fov '/' mouseid '_' expname '/group_proc/'...
                    groupreg_method '_' mcorr_method '_' segment_method '/'...
                    mouseid '_' expname '_imreg_ref' reffile '/'];
        else
            dir_proc = [data_locn 'Analysis/' mouseid '/' mouseid '_' expname '/group_proc/'...
                    groupreg_method '_' mcorr_method '_' segment_method '/'...
                    mouseid '_' expname '_imreg_ref' reffile '/'];
        end
        
        check = zeros(1,6);
        if exist(dir_proc,'dir')
            % 1) Check for existing collective roi segmentation output 
            if exist([dir_proc  mouseid '_' expname '_ref' reffile '_segment_output.mat'],'file')
                check(1) = 1;
            end
            
            % 2) Check for existing consolidated fissa data 
            if exist([dir_proc '/fissa/' mouseid '_' expname '_ref' reffile '_fissa_output.mat'],'file')
                check(2) = 1;
            end
            
            % 3) Check for existing consolidated spike data 
            if exist([dir_proc '/' str_fissa '/bl_prctile' num2str(bl_prctile) '/' mouseid '_' expname '_ref' reffile '_spikes.mat'],'file')
                check(3) = 1;
            end
            
            % 4) Check for existing consolidated tracking data 
            if exist([data_locn 'Analysis/' mouseid '/' mouseid '_' expname '/group_proc/'...
                      mouseid '_' expname '_downTrackdata.mat'],'file')
                check(4) = 1;
            end
            
            % 5) Check for existing collective PF mapping output
            if exist([dir_proc '/' str_fissa '/bl_prctile' num2str(bl_prctile) '/' mouseid '_' expname '_ref' reffile '_PFmap_output.mat'],'file')
                check(5) = 1;
            end
            
            % 6) Check if mat file for all proc data for the file exists
            if exist([dir_proc  mouseid '_' expname '_ref' reffile '_' mcorr_method '_' segment_method '_' str_fissa '_allData_blprctile' num2str(bl_prctile) '.mat'],'file')
                check(6) = 1;
            end
        end
    else % individual file 
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
            if exist([dir_segment '/' str_fissa '/bl_prctile' num2str(bl_prctile) '/' file '_spikes_output.mat'],'file') || ...
               exist([dir_segment '/' str_fissa '/' file '_spikes.mat'],'file') 
                check(4) = 1;
            end

            % 5) Check for tracking data
            if exist([dir_proc 'behaviour/' file '_downTrackdata.mat'],'file')
                check(5) = 1;
            end

            % 6) Check for existing PF mapping output
            if exist([dir_segment '/' str_fissa '/bl_prctile' num2str(bl_prctile) '/' file '_PFmap_output.mat'],'file')
                check(6) = 1;
            end

            % 7) Check if mat file for all proc data for the file exists
            if exist([dir_proc file '_' mcorr_method '_' segment_method '_' str_fissa '_allData_blprctile' num2str(bl_prctile) '.mat'],'file')
                check(7) = 1;
            end
        end
    end
    
end