% Written by Ann Go
%
% Function to determine whether file has been processed using same
% parameters as specified by user
% This can be used to check for processed data for individual files OR for
% processed data for different experiments of a mouse.

function check = checkforExistingProcData(data_locn, text, params_methods, reffile)

    %if nargin<4, see line 39
    if nargin<3 
        mcorr_method = 'normcorre-nr';
        segment_method = 'CaImAn';
        dofissa = true;
    else
        mcorr_method = params_methods.mcorr_method;
        segment_method = params_methods.segment_method;
        dofissa = params_methods.dofissa;
    end
    
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
        [mouseid,expname] = find_mouseIDexpname(list);
        listfile = [data_locn 'Digital_Logbook/lists/' list];
        files = extractFilenamesFromTxtfile(listfile);
        if nargin<4, reffile = files(1,:); end
        
        dir_proc = [data_locn 'Analysis/' mouseid '/' mouseid '_' expname '/group_proc/'...
                    mcorr_method '_' segment_method '_' str_fissa '/'...
                    mouseid '_' expname '_imreg_ref' reffile '/'];
        
        check = zeros(1,4);
        if exist(dir_proc,'dir')
            % 1) Check for existing collective roi segmentation output 
            if exist([dir_proc  mouseid '_' expname '_ref' reffile '_segment_output.mat'],'file')
                check(1) = 1;
            end
            
            % 2) Check for existing consolidated data for fissa, spikes, tracking
            if exist([dir_proc  mouseid '_' expname '_ref' reffile '_fissa_spike_track_data.mat'],'file')
                check(2) = 1;
            end
            
            % 3) Check for existing collective PF mapping output
            if exist([dir_proc  mouseid '_' expname '_ref' reffile '_PFmap_output.mat'],'file')
                check(3) = 1;
            end
            
            % 4) Check if mat file for all proc data for the file exists
            if exist([dir_proc  mouseid '_' expname '_ref' reffile '_allData.mat'],'file')
                check(4) = 1;
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
            if exist([dir_fissa file '_spikes_output.mat'],'file')
                check(4) = 1;
            end

            % 5) Check for tracking data
            if exist([dir_proc 'behaviour/' file '_downTrackdata.mat'],'file')
                check(5) = 1;
            end

            % 6) Check for existing PF mapping output
            if exist([dir_fissa 'PFmaps/' file '_PFmap_output.mat'],'file')
                check(6) = 1;
            end

            % 7) Check if mat file for all proc data for the file exists
            if exist([dir_proc file '_' mcorr_method '_' segment_method '_' str_fissa '_allData.mat'],'file')
                check(7) = 1;
            end
        end
    end
    
end