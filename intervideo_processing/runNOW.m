% nov day2
for n = 1:size(files_novd2,1)
    str = sprintf('Processing %s ...\n', files_novd2(n,:));
    cprintf(str)
%     load([filedir1 files_novd2(n,:) '/' files_novd2(n,:) '_2P_mcorr_output.mat'])
%     novd2(n).green = green;
%     novd2(n).red = red;
%     clear green red params

    novd2_imG = read_file( [filedir1 files_novd2(n,:) '/' files_novd2(n,:) '_2P_XYT_green_mcorr.tif'] );
    novd2_imR = read_file( [filedir1 files_novd2(n,:) '/' files_novd2(n,:) '_2P_XYT_red_mcorr.tif'] );

    [ novd2(n).shift, novd2(n).red.globalreg_meanframe ] = globalregisterImage(ref.red.meanregframe, novd2(n).red.meanregframe, 1 );
    [ novd2_imGglobalreg, novd2_imRglobalreg ] = ...
        globalregisterStack( novd2_imG, novd2_imR, novd2(n).shift(1), novd2(n).shift(2) );
    clear novd2_imG novd2_imR
    
    fname_tif_novd2_imGglobalreg = [filedir1 files_novd2(n,:) '/' files_novd2(n,:) '_2P_XYT_green_mcorr_globalreg.tif'];
        writeTifStack( novd2_imGglobalreg, fname_tif_novd2_imGglobalreg );
    fname_tif_novd2_imRglobalreg = [filedir1 files_novd2(n,:) '/' files_novd2(n,:) '_2P_XYT_red_mcorr_globalreg.tif'];
        writeTifStack( novd2_imRglobalreg, fname_tif_novd2_imRglobalreg );
    clear fname_tif_novd2_imGglobalreg fname_tif_novd2_imRglobalreg
        
    for i = 1:Nmasks
        maskind = masks(:,:,i);
        for j = 1:size(novd2_imGglobalreg,3)
            imG_reshaped = reshape( novd2_imGglobalreg(:,:,j), 512*512, 1);
            novd2(n).tsG( i, j ) = mean( imG_reshaped(maskind) );
            imR_reshaped = reshape( novd2_imRglobalreg(:,:,j), 512*512, 1);
            novd2(n).tsR( i, j ) = mean( imR_reshaped(maskind) );
        end
    end
    clear novd2_imGglobalreg novd2_imRglobalreg imG_reshaped imR_reshaped maskind i j

    novd2(n).R = ratiometric_Ca( novd2(n).tsG, novd2(n).tsR, 11 );
    novd2(n).spikes = nndORoasis(novd2(n).R, 2, 0.94, 2.4);
end

% fam1 day3
for n = 1:size(files_fam1d3,1)
    str = sprintf('Processing %s ...\n', files_fam1d3(n,:));
    cprintf(str)
%     load([filedir2 files_fam1d3(n,:) '/' files_fam1d3(n,:) '_2P_mcorr_output.mat'])
%     fam1d3(n).green = green;
%     fam1d3(n).red = red;
%     clear green red 

    fam1d3_imG = read_file( [filedir2 files_fam1d3(n,:) '/' files_fam1d3(n,:) '_2P_XYT_green_mcorr.tif'] );
    fam1d3_imR = read_file( [filedir2 files_fam1d3(n,:) '/' files_fam1d3(n,:) '_2P_XYT_red_mcorr.tif'] );

    [ fam1d3(n).shift, fam1d3(n).red.globalreg_meanframe ] = globalregisterImage(ref.red.meanregframe, fam1d3(n).red.meanregframe, 1 );
    [ fam1d3_imGglobalreg, fam1d3_imRglobalreg ] = ...
        globalregisterStack( fam1d3_imG, fam1d3_imR, fam1d3(n).shift(1) , fam1d3(n).shift(2) );
    clear fam1d3_imG fam1d3_imR
    
    fname_tif_fam1d3_imGglobalreg = [filedir2 files_fam1d3(n,:) '/' files_fam1d3(n,:) '_2P_XYT_green_mcorr_globalreg.tif'];
        writeTifStack( fam1d3_imGglobalreg, fname_tif_fam1d3_imGglobalreg );
    fname_tif_fam1d3_imRglobalreg = [filedir2 files_fam1d3(n,:) '/' files_fam1d3(n,:) '_2P_XYT_red_mcorr_globalreg.tif'];
        writeTifStack( fam1d3_imRglobalreg, fname_tif_fam1d3_imRglobalreg );
    clear fname_tif_fam1d3_imGglobalreg fname_tif_fam1d3_imRglobalreg
        
    for i = 1:Nmasks
        maskind = masks(:,:,i);
        for j = 1:size(fam1d3_imGglobalreg,3)
            imG_reshaped = reshape( fam1d3_imGglobalreg(:,:,j), 512*512, 1);
            fam1d3(n).tsG( i, j ) = mean( imG_reshaped(maskind) );
            imR_reshaped = reshape( fam1d3_imRglobalreg(:,:,j), 512*512, 1);
            fam1d3(n).tsR( i, j ) = mean( imR_reshaped(maskind) );
        end
    end
    clear fam1d3_imGglobalreg fam1d3_imRglobalreg imG_reshaped imR_reshaped maskind i j

    fam1d3(n).R = ratiometric_Ca( fam1d3(n).tsG, fam1d3(n).tsR, 11 );
    fam1d3(n).spikes = nndORoasis(fam1d3(n).R, 2, 0.94, 2.4);
end

% fam1rev day3
for n = 1:size(files_fam1revd3,1)
    str = sprintf('Processing %s ...\n', files_fam1revd3(n,:));
    cprintf(str)
%     load([filedir2 files_fam1revd3(n,:) '/' files_fam1revd3(n,:) '_2P_mcorr_output.mat'])
%     fam1revd3(n).green = green;
%     fam1revd3(n).red = red;
%     clear green red 

    fam1revd3_imG = read_file( [filedir2 files_fam1revd3(n,:) '/' files_fam1revd3(n,:) '_2P_XYT_green_mcorr.tif'] );
    fam1revd3_imR = read_file( [filedir2 files_fam1revd3(n,:) '/' files_fam1revd3(n,:) '_2P_XYT_red_mcorr.tif'] );

    [ fam1revd3(n).shift, fam1revd3(n).red.globalreg_meanframe ] = globalregisterImage(ref.red.meanregframe, fam1revd3(n).red.meanregframe, 1 );
    [ fam1revd3_imGglobalreg, fam1revd3_imRglobalreg ] = ...
        globalregisterStack( fam1revd3_imG, fam1revd3_imR, fam1revd3(n).shift(1) , fam1revd3(n).shift(2) );
    clear fam1revd3_imG fam1revd3_imR
    
    fname_tif_fam1revd3_imGglobalreg = [filedir2 files_fam1revd3(n,:) '/' files_fam1revd3(n,:) '_2P_XYT_green_mcorr_globalreg.tif'];
        writeTifStack( fam1revd3_imGglobalreg, fname_tif_fam1revd3_imGglobalreg );
    fname_tif_fam1revd3_imRglobalreg = [filedir2 files_fam1revd3(n,:) '/' files_fam1revd3(n,:) '_2P_XYT_red_mcorr_globalreg.tif'];
        writeTifStack( fam1revd3_imRglobalreg, fname_tif_fam1revd3_imRglobalreg );
    clear fname_tif_fam1revd3_imGglobalreg fname_tif_fam1revd3_imRglobalreg
        
    for i = 1:Nmasks
        maskind = masks(:,:,i);
        for j = 1:size(fam1revd3_imGglobalreg,3)
            imG_reshaped = reshape( fam1revd3_imGglobalreg(:,:,j), 512*512, 1);
            fam1revd3(n).tsG( i, j ) = mean( imG_reshaped(maskind) );
            imR_reshaped = reshape( fam1revd3_imRglobalreg(:,:,j), 512*512, 1);
            fam1revd3(n).tsR( i, j ) = mean( imR_reshaped(maskind) );
        end
    end
    clear fam1revd3_imGglobalreg fam1revd3_imRglobalreg imG_reshaped imR_reshaped maskind i j

    fam1revd3(n).R = ratiometric_Ca( fam1revd3(n).tsG, fam1revd3(n).tsR, 11 );
    fam1revd3(n).spikes = nndORoasis(fam1revd3(n).R, 2, 0.94, 2.4);
end
clear Nmasks data_locn filedir1 filedir2 filedir3 files_fam1d2 files_novd2 files_fam1d3 files_fam1revd3
clear file_ref str n i j mean_imratio


%% Save output
fam1d2_0946 = ref;

fam1d2_0933 = fam1d2(1);
fam1d2_0937 = fam1d2(2);
fam1d2_0951 = fam1d2(3);
clear fam1d2

novd2_0959 = novd2(1);
novd2_1004 = novd2(2);
clear novd2

fam1d3_0922 = fam1d3(1);
fam1d3_0944 = fam1d3(2);
fam1d3_0949 = fam1d3(3);
clear fam1d3

fam1revd3_0957 = fam1revd3(1);
fam1revd3_1007 = fam1revd3(2);
clear fam1revd3

save('m62_fam1_nov_fam1rev_20181015to16_segment.mat')