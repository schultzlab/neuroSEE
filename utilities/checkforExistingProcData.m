% Written by Ann Go
% Function to determine whether file has been processed using same
% parameters as specified by user

function yn = checkforExistingProcData(data_locn, file, mcorr_method, segment_method, fissa_yn)
    if fissa_yn == 1
        str_fissa = '_FISSA';
    else
        str_fissa = 'noFISSA';
    end
    
    yn = 0;
    filedir = [data_locn 'Data/' file(1:8) '/Processed/' file '/'];
    matfiles = dir(fullfile(filedir,['*.','mat']));
    if numel(matfiles) > 0 
        for i = 1:numel(matfiles)
            name = matfiles(i).name;
            if contains(name,'allData')
                if contains(name,mcorr_method)
                    if contains(name,segment_method)
                        if contains(name,str_fissa)
                            yn = 1;
                        end
                    end
                end
            end
        end
    end
end