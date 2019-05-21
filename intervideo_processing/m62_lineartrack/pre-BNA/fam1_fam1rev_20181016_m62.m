data_locn = '/Volumes/RDS/project/thefarm2/live/CrazyEights/AD_2PCa/Data/20181016/Processed/';

%% Compare red images to see which ones are most similar
figure;
subplot(251);imagesc(fam1_0909.red.meanregframe); title('Fam1 09:09'); caxis([2000 2400]); axis off; axis square;
subplot(252);imagesc(fam1_0914.red.meanregframe); title('Fam1 09:14'); caxis([2000 2400]); axis off; axis square;
subplot(253);imagesc(fam1_0918.red.meanregframe); title('Fam1 09:18'); caxis([2000 2400]); axis off; axis square;
subplot(254);imagesc(fam1_0922.red.meanregframe); title('Fam1 09:22'); caxis([2000 2400]); axis off; axis square;
subplot(255);imagesc(fam1_0944.red.meanregframe); title('Fam1 09:44'); caxis([2000 2400]); axis off; axis square;
subplot(256);imagesc(fam1_0949.red.meanregframe); title('Fam1 09:49'); caxis([2000 2400]); axis off; axis square;
subplot(257);imagesc(fam1rev_0957.red.meanregframe); title('Fam1rev 09:57'); caxis([2000 2400]); axis off; axis square;
subplot(258);imagesc(fam1rev_1007.red.meanregframe); title('Fam1rev 10:07'); caxis([2000 2400]); axis off; axis square;
subplot(259);imagesc(fam1rev_1011.red.meanregframe); title('Fam1rev 10:11'); caxis([2000 2400]); axis off; axis square;
subplot(2,5,10);imagesc(fam1rev_1015.red.meanregframe); title('Fam1rev 10:15'); caxis([2000 2400]); axis off; axis square;
 

%% Compare ROIS to determine which image to use as reference
figure;
subplot(251); imagesc(fam1_0922.mean_imratio); colormap(jet); axis square; title('Fam1 09:22'); caxis([1 6]);
subplot(252); imagesc(fam1_0944.mean_imratio); colormap(jet); axis square; title('Fam1 09:44'); caxis([1 6]); 
subplot(253); imagesc(fam1_0949.mean_imratio); colormap(jet); axis square; title('Fam1 09:49'); caxis([1 6]); 
subplot(254); imagesc(fam1rev_0957.mean_imratio); colormap(jet); axis square; title('Fam1rev 09:57'); caxis([1 6]); 
subplot(255); imagesc(fam1rev_1007.mean_imratio); colormap(jet); axis square; title('Fam1rev 10:07'); caxis([1 6]); 

subplot(256); imagesc(fam1_0922.mean_imratio); caxis([1 6]); colormap(jet); axis square; 
    hold on
        for j = 1:size(fam1_0922.masks,3)
            outline = bwboundaries(fam1_0922.masks(:,:,j));
            trace = outline{1};
            plot(trace(:,2),trace(:,1),'w','Linewidth',2);
        end
    hold off;
subplot(257); imagesc(fam1_0944.mean_imratio); caxis([1 6]); title('Fam1 09:44'); colormap(jet); axis square; 
    hold on
        for j = 1:size(fam1_0944.masks,3)
            outline = bwboundaries(fam1_0944.masks(:,:,j));
            trace = outline{1};
            plot(trace(:,2),trace(:,1),'w','Linewidth',2);
        end
    hold off
subplot(258); imagesc(fam1_0949.mean_imratio); caxis([1 6]); title('Fam1 09:49'); colormap(jet); axis square; 
    hold on
        for j = 1:size(fam1_0949.masks,3)
            outline = bwboundaries(fam1_0949.masks(:,:,j));
            trace = outline{1};
            plot(trace(:,2),trace(:,1),'w','Linewidth',2);
        end
    hold off
subplot(259); imagesc(fam1rev_0957.mean_imratio); caxis([1 6]); colormap(jet); axis square; 
    hold on
        for j = 1:size(fam1rev_0957.masks,3)
            outline = bwboundaries(fam1rev_0957.masks(:,:,j));
            trace = outline{1};
            plot(trace(:,2),trace(:,1),'w','Linewidth',2);
        end
    hold off
subplot(2,5,10); imagesc(fam1rev_1007.mean_imratio); caxis([1 6]); title('Fam1rev 10:07'); colormap(jet); axis square; 
    hold on
        for j = 1:size(fam1rev_1007.masks,3)
            outline = bwboundaries(fam1rev_1007.masks(:,:,j));
            trace = outline{1};
            plot(trace(:,2),trace(:,1),'w','Linewidth',2);
        end
    hold off

    
%% Refine ROIs for reference image
% recalculate R & spikes because I changed smoothing (so that it's the last step in the calculation of R)
fam1rev_1007.R = ratiometric_Ca( fam1rev_1007.tsG, fam1rev_1007.tsR, 11 );
fam1rev_1007.spikes = nndORoasis(fam1rev_1007.R, 2, 0.94, 2.4);
GUI_viewROIsSpikes_withModeChoice(fam1rev_1007.mean_imratio,fam1rev_1007.masks,fam1rev_1007.tsG,fam1rev_1007.R,fam1rev_1007.spikes,fam1rev_1007.spikes)

% refine ROIs
Numcells = size(fam1rev_1007.masks,3);
areaThr = 90;
maskArea = zeros(1,size(fam1rev_1007.masks,3));
for i = 1:Numcells
    maskArea(i) = bwarea(fam1rev_1007.masks(:,:,i));
end
i_area = find(maskArea>=areaThr);

% snrThr = 490;
% maskSNR = zeros(1,Numcells);
% for i = 1:Numcells
%     y = tsG(i,:);
%     maskSNR(i) = GetSn(y,[0.25,0.5],'logmexp');
% end
% i_snr = find(maskSNR<=snrThr);

satThr = 0.95; satValue = 4000; satTime = 0.3000;
Fsat = (fam1rev_1007.tsG >= satThr*satValue);
i_sat = find(mean(Fsat,ndims(fam1rev_1007.tsG))<satTime);
ind = intersect(i_area,i_sat);

masksftuned = fam1rev_1007.masks(:,:,ind);
tsGftuned = fam1rev_1007.tsG(ind,:);
tsRftuned = fam1rev_1007.tsR(ind,:);
Rftuned = fam1rev_1007.R(ind,:);
spikesftuned = fam1rev_1007.spikes(ind,:);

ind_manual = [3,4,5,7,11,25,31,38,51,52,59,63,76,78,81,84,96,99,103,104,106,108,111,...
              119,120,122,132,135,139,143,144,145,153,158,160,173,176,184,185,188,...
              190,191,195,197,204,205,207,212,213,216];
ind = 1:size(masksftuned,3);
for i = 1:length(ind_manual)
    ind = ind(ind~=ind_manual(i));
end
fam1rev_1007.masksftuned = masksftuned(:,:,ind);
fam1rev_1007.tsGftuned = tsGftuned(ind,:);
fam1rev_1007.tsRftuned = tsRftuned(ind,:);
fam1rev_1007.Rftuned = Rftuned(ind,:);
fam1rev_1007.spikesftuned = spikesftuned(ind,:);

clear Numcells areaThr maskArea i_area satThr Fsat i_sat ind masksftuned tsGftuned tsRftuned spikesftuned ind_manual 

%% Regenerate PF maps for reference image
file = '20181016_10_07_07';
fam1rev_1007.trackdata.x = x;
fam1rev_1007.trackdata.y = y;
fam1rev_1007.trackdata.r = r;
fam1rev_1007.trackdata.phi = phi;
fam1rev_1007.trackdata.speed = speed;
fam1rev_1007.trackdata.time = time;

[ fam1rev_1007.occMap, fam1rev_1007.spikeMap, fam1rev_1007.infoMap, fam1rev_1007.placeMap, fam1rev_1007.downData,...
    fam1rev_1007.activeData, fam1rev_1007.placeMap_smooth ] = ...
    generatePFmap( fam1rev_1007.spikesftuned, [], fam1rev_1007.trackdata, 30.9, 10, 1, 2, 180, 1, 5);
[ fam1rev_1007.sorted_placeMap, fam1rev_1007.normsorted_placeMap, fam1rev_1007.sortIdx ] = ...
    sortPFmap( fam1rev_1007.placeMap_smooth, fam1rev_1007.infoMap, 1);

Nbins = 180;
figure('Position',[1087 648 500 800]);
    subplot(9,5,2:5); imagesc(fam1rev_1007.occMap);
        % xticks(Nbins/6:Nbins/6:Nbins); yticks([]);
        xticks([]); yticks([]);
        title('Occupancy map'); % colorbar;
    subplot(9,5,[6,11,16,21]); imagesc(fam1rev_1007.infoMap(fam1rev_1007.sortIdx,2));
        xticks([]); yticks([]); ylabel('Cells'); 
        title('Info map'); colorbar;
    subplot(9,5,[7:10,12:15,17:20,22:25]);
        nspikeMap = fam1rev_1007.spikeMap./max(max(fam1rev_1007.spikeMap));
        imagesc(nspikeMap(fam1rev_1007.sortIdx,:)); xticks([]); yticks([]);
        title('Spike maps');
    subplot(9,5,[27:30,32:35,37:40,42:45]);    
        imagesc(fam1rev_1007.sorted_placeMap); 
        xticks(Nbins/6:Nbins/6:Nbins); yticks([]);
        degperbin = 360/Nbins; xticklabels(degperbin*(Nbins/6:Nbins/6:Nbins));
        colorbar('Location','southoutside'); caxis([0,0.01]);
        %xticks([1 30 60 90 120 150 180]); xticklabels([1 60 120 180 240 300 360]); 
        % yticklabels(sortIdx);
        title('Place field maps');
        xlabel('Position (degrees)'); ylabel('Cells');
        
%% Register at least one other image to reference image
% fam1_0944
[ fam1_0944.shift, fam1_0944.red.globalreg_meanframe ] = globalregisterImage(fam1rev_1007.red.meanregframe, fam1_0944.red.meanregframe, 1 );
fam1_0944_imG = read_file( [data_locn '20181016_09_44_06/20181016_09_44_06_2P_XYT_green_mcorr.tif'] );
fam1_0944_imR = read_file( [data_locn '20181016_09_44_06/20181016_09_44_06_2P_XYT_red_mcorr.tif'] );

[ fam1_0944_imGglobalreg, fam1_0944_imRglobalreg ] = ...
    globalregisterStack( fam1_0944_imG, fam1_0944_imR, fam1_0944.shift(1) , fam1_0944.shift(2) );
[ ~, fam1_0944.globalreg_mean_imratio ] = globalregisterImage(fam1rev_1007.mean_imratio, fam1_0944.mean_imratio, 1 );
fam1_0944.globalregmasks = registerMasks( fam1_0944.masks, fam1_0944.shift(1) , fam1_0944.shift(2) );

% fam1_0949
[ fam1_0949.shift, fam1_0949.red.globalreg_meanframe ] = globalregisterImage(fam1rev_1007.red.meanregframe, fam1_0949.red.meanregframe, 1 );
fam1_0949_imG = read_file( [data_locn '20181016_09_49_06/20181016_09_49_06_2P_XYT_green_mcorr.tif'] );
fam1_0949_imR = read_file( [data_locn '20181016_09_49_06/20181016_09_49_06_2P_XYT_red_mcorr.tif'] );

[ fam1_0949_imGglobalreg, fam1_0949_imRglobalreg ] = ...
    globalregisterStack( fam1_0949_imG, fam1_0949_imR, fam1_0949.shift(1) , fam1_0949.shift(2) );
[ ~, fam1_0949.globalreg_mean_imratio ] = globalregisterImage(fam1rev_1007.mean_imratio, fam1_0949.mean_imratio, 1 );
fam1_0949.globalregmasks = registerMasks( fam1_0949.masks, fam1_0949.shift(1) , fam1_0949.shift(2) );


%% Add any missing ROIs from registered image
fam1_0944.R = ratiometric_Ca( fam1_0944.tsG, fam1_0944.tsR, 11 );
fam1_0944.spikes = nndORoasis(fam1_0944.R, 2, 0.94, 2.4);

GUI_viewROIsSpikes_withModeChoice(fam1_0944.globalreg_mean_imratio,fam1_0944.globalregmasks,fam1_0944.tsG,fam1_0944.R,fam1_0944.spikes,fam1_0944.spikes)

ind = [1,4,24,31,41,43,57,61,63,64,69,74,86,93,94,97,98,103,104,111,112,118,119,120,135,136,139,142,148,157,162,167,171,181,198,203,208,215,225];

globalmasks(:,:,1:size(fam1rev_1007.masksftuned,3)) = fam1rev_1007.masksftuned;
globalmasks(:,:,size(fam1rev_1007.masksftuned,3)+1 : size(fam1rev_1007.masksftuned,3)+length(ind)) = fam1_0944.masks(:,:,ind);

%% Show reference image with final ROIs

figure; 
subplot(131); imagesc(fam1rev_1007.mean_imratio); colormap(jet); 
    axis square; caxis([1 6]); title('Fam1rev 10:07'); axis off;
    hold on
        for j = 1:size(fam1rev_1007.masksftuned,3)
            outline = bwboundaries(fam1rev_1007.masksftuned(:,:,j));
            trace = outline{1};
            plot(trace(:,2),trace(:,1),'w','Linewidth',2);
        end
    hold off
subplot(132); imagesc(fam1_0944.globalreg_mean_imratio); colormap(jet); 
    axis square; caxis([1 6]); title('Fam1 09:44'); axis off;
    hold on
        for j = 1:size(fam1_0944.globalregmasks,3)
            outline = bwboundaries(fam1_0944.globalregmasks(:,:,j));
            if size(outline,1) ~= 0
                trace = outline{1};
                plot(trace(:,2),trace(:,1),'w','Linewidth',2);
            end
        end
    hold off
subplot(133); imagesc(fam1rev_1007.mean_imratio); colormap(jet); 
    axis square; caxis([1 6]); title('Global ROIs'); axis off;
    hold on
        for j = 1:size(globalmasks,3)
            outline = bwboundaries(globalmasks(:,:,j));
            trace = outline{1};
            plot(trace(:,2),trace(:,1),'w','Linewidth',2);
        end
    hold off
        

%% Recalculate tsG, tsR, R, spikes for additional ROIs in reference image and regenerate PFmap
fam1rev_1007.globaltsG(1:size(fam1rev_1007.tsGftuned,1),:) = fam1rev_1007.tsGftuned;
fam1rev_1007.globaltsR(1:size(fam1rev_1007.tsRftuned,1),:) = fam1rev_1007.tsRftuned;

fam1rev_1007_imG = read_file( [data_locn '20181016_10_07_07/20181016_10_07_07_2P_XYT_green_mcorr.tif'] );
fam1rev_1007_imR = read_file( [data_locn '20181016_10_07_07/20181016_10_07_07_2P_XYT_red_mcorr.tif'] );

for i = size(fam1rev_1007.tsGftuned,1)+1:size(globalmasks,3)
    maskind = find(globalmasks(:,:,i));
    for j = 1:size(fam1rev_1007_imG,3)
        imG_reshaped = reshape( fam1rev_1007_imG(:,:,j), 512*512, 1);
        fam1rev_1007.globaltsG( i, j ) = mean( imG_reshaped(maskind) );
        imR_reshaped = reshape( fam1rev_1007_imR(:,:,j), 512*512, 1);
        fam1rev_1007.globaltsR( i, j ) = mean( imR_reshaped(maskind) );
    end
end

fam1rev_1007.globalR = ratiometric_Ca( fam1rev_1007.globaltsG, fam1rev_1007.globaltsR, 11 );
fam1rev_1007.globalspikes = nndORoasis(fam1rev_1007.globalR, 2, 0.94, 2.4);

% histogram
% [ fam1rev_1007.globaloccMap, fam1rev_1007.globalspikeMap, fam1rev_1007.globalinfoMap, fam1rev_1007.globalplaceMap, fam1rev_1007.globaldownData,...
%     fam1rev_1007.globalactiveData, fam1rev_1007.globalplaceMap_smooth ] = ...
%     generatePFmap( fam1rev_1007.globalspikes, [], fam1rev_1007.trackdata, 30.9, 10, 1, 2, 180, 1, 5);
% [ fam1rev_1007.globalsorted_placeMap, fam1rev_1007.globalnormsorted_placeMap, fam1rev_1007.globalsortIdx ] = ...
%     sortPFmap( fam1rev_1007.globalplaceMap_smooth, fam1rev_1007.globalinfoMap, 1);

% ASD
[ fam1rev_1007.globaloccMap, fam1rev_1007.globalspikeMap, fam1rev_1007.globalinfoMap, fam1rev_1007.globalplaceMap, fam1rev_1007.globaldownData,...
    fam1rev_1007.globalactiveData, fam1rev_1007.globalplaceMap_smooth ] = ...
    generatePFmap( fam1rev_1007.globalspikes, [], fam1rev_1007.trackdata, 30.9, 10, 1, 1, 180, 1, 5);
[ fam1rev_1007.globalsorted_placeMap, fam1rev_1007.globalnormsorted_placeMap, fam1rev_1007.globalsortIdx ] = ...
    sortPFmap( fam1rev_1007.globalplaceMap, fam1rev_1007.globalinfoMap, 1);

Nbins = 180;
figure('Position',[1087 648 500 800]);
    subplot(9,5,2:5); imagesc(fam1rev_1007.globaloccMap);
        % xticks(Nbins/6:Nbins/6:Nbins); yticks([]);
        xticks([]); yticks([]);
        title('Occupancy map'); % colorbar;
    subplot(9,5,[6,11,16,21]); imagesc(fam1rev_1007.globalinfoMap(fam1rev_1007.globalsortIdx,2));
        xticks([]); yticks([]); ylabel('Cells'); 
        title('Info map'); colorbar;
    subplot(9,5,[7:10,12:15,17:20,22:25]);
        nspikeMap = fam1rev_1007.globalspikeMap./max(max(fam1rev_1007.globalspikeMap));
        imagesc(nspikeMap(fam1rev_1007.globalsortIdx,:)); xticks([]); yticks([]);
        title('Spike maps');
    subplot(9,5,[27:30,32:35,37:40,42:45]);    
        imagesc(fam1rev_1007.globalsorted_placeMap); 
        xticks(Nbins/6:Nbins/6:Nbins); yticks([]);
        xticklabels(360/Nbins*(Nbins/6:Nbins/6:Nbins));
        colorbar('Location','southoutside'); caxis([0,0.01]);
        title('Place field maps');
        xlabel('Position (degrees)'); ylabel('Cells');

GUI_viewROIsSpikes_withModeChoice(fam1rev_1007.mean_imratio,globalmasks,fam1rev_1007.globaltsG,fam1rev_1007.globalR,fam1rev_1007.globalspikes,fam1rev_1007.globalspikes)


%% Recalculate tsG, tsR, R, spikes using global ROIs in other images and generate PFmaps
% fam1_0944
for i = 1:size(globalmasks,3)
    maskind = find(globalmasks(:,:,i));
    for j = 1:size(fam1_0944_imG,3)
        imG_reshaped = reshape( fam1_0944_imGglobalreg(:,:,j), 512*512, 1);
        fam1_0944.globaltsG( i, j ) = mean( imG_reshaped(maskind) );
        imR_reshaped = reshape( fam1_0944_imRglobalreg(:,:,j), 512*512, 1);
        fam1_0944.globaltsR( i, j ) = mean( imR_reshaped(maskind) );
    end
end

fam1_0944.globalR = ratiometric_Ca( fam1_0944.globaltsG, fam1_0944.globaltsR, 11 );
fam1_0944.globalspikes = nndORoasis(fam1_0944.globalR, 2, 0.94, 2.4);

fam1_0944.trackdata.x = x;
fam1_0944.trackdata.y = y;
fam1_0944.trackdata.r = r;
fam1_0944.trackdata.phi = phi;
fam1_0944.trackdata.speed = speed;
fam1_0944.trackdata.time = time;

clear filename phi r speed time TTLout w x y

% histogram
% [ fam1_0944.globaloccMap, fam1_0944.globalspikeMap, fam1_0944.globalinfoMap, fam1_0944.globalplaceMap, fam1_0944.globaldownData,...
%     fam1_0944.globalactiveData, fam1_0944.globalplaceMap_smooth ] = ...
%     generatePFmap( fam1_0944.globalspikes, [], fam1_0944.trackdata, 30.9, 10, 1, 2, 180, 1, 5);
% [ fam1_0944.globalsorted_placeMap, fam1_0944.globalnormsorted_placeMap, fam1_0944.globalsortIdx ] = ...
%     sortPFmap( fam1_0944.globalplaceMap_smooth, fam1_0944.globalinfoMap, 1);

% ASD
[ fam1_0944.globaloccMap, fam1_0944.globalspikeMap, fam1_0944.globalinfoMap, fam1_0944.globalplaceMap, fam1_0944.globaldownData,...
    fam1_0944.globalactiveData, fam1_0944.globalplaceMap_smooth ] = ...
    generatePFmap( fam1_0944.globalspikes, [], fam1_0944.trackdata, 30.9, 10, 1, 1, 180, 1, 5);
[ fam1_0944.globalsorted_placeMap, fam1_0944.globalnormsorted_placeMap, fam1_0944.globalsortIdx ] = ...
    sortPFmap( fam1_0944.globalplaceMap, fam1_0944.globalinfoMap, 1);


Nbins = 180;
figure('Position',[1087 648 500 800]);
    subplot(9,5,2:5); imagesc(fam1_0944.globaloccMap);
        % xticks(Nbins/6:Nbins/6:Nbins); yticks([]);
        xticks([]); yticks([]);
        title('Occupancy map'); % colorbar;
    subplot(9,5,[6,11,16,21]); imagesc(fam1_0944.globalinfoMap(fam1_0944.globalsortIdx,2));
        xticks([]); yticks([]); ylabel('Cells'); 
        title('Info map'); colorbar;
    subplot(9,5,[7:10,12:15,17:20,22:25]);
        nspikeMap = fam1_0944.globalspikeMap./max(max(fam1_0944.globalspikeMap));
        imagesc(nspikeMap(fam1_0944.globalsortIdx,:)); xticks([]); yticks([]);
        title('Spike maps');
    subplot(9,5,[27:30,32:35,37:40,42:45]);    
        imagesc(fam1_0944.globalsorted_placeMap); 
        xticks(Nbins/6:Nbins/6:Nbins); yticks([]);
        xticklabels(360/Nbins*(Nbins/6:Nbins/6:Nbins));
        colorbar('Location','southoutside'); caxis([0,0.01]);
        title('Place field maps');
        xlabel('Position (degrees)'); ylabel('Cells');
clear Nbins nspikeMap

% fam1_0949
for i = 1:size(globalmasks,3)
    maskind = find(globalmasks(:,:,i));
    for j = 1:size(fam1_0949_imGglobalreg,3)
        imG_reshaped = reshape( fam1_0949_imGglobalreg(:,:,j), 512*512, 1);
        fam1_0949.globaltsG( i, j ) = mean( imG_reshaped(maskind) );
        imR_reshaped = reshape( fam1_0949_imRglobalreg(:,:,j), 512*512, 1);
        fam1_0949.globaltsR( i, j ) = mean( imR_reshaped(maskind) );
    end
end
clear i j imG_reshaped imR_reshaped maskind

fam1_0949.globalR = ratiometric_Ca( fam1_0949.globaltsG, fam1_0949.globaltsR, 11 );
fam1_0949.globalspikes = nndORoasis(fam1_0949.globalR, 2, 0.94, 2.4);

fam1_0949.trackdata.x = x;
fam1_0949.trackdata.y = y;
fam1_0949.trackdata.r = r;
fam1_0949.trackdata.phi = phi;
fam1_0949.trackdata.speed = speed;
fam1_0949.trackdata.time = time;

clear filename phi r speed time TTLout w x y

% histogram
% [ fam1_0949.globaloccMap, fam1_0949.globalspikeMap, fam1_0949.globalinfoMap, fam1_0949.globalplaceMap, fam1_0949.globaldownData,...
%     fam1_0949.globalactiveData, fam1_0949.globalplaceMap_smooth ] = ...
%     generatePFmap( fam1_0949.globalspikes, [], fam1_0949.trackdata, 30.9, 10, 1, 2, 180, 1, 5);
% [ fam1_0949.globalsorted_placeMap, fam1_0949.globalnormsorted_placeMap, fam1_0949.globalsortIdx ] = ...
%     sortPFmap( fam1_0949.globalplaceMap_smooth, fam1_0949.globalinfoMap, 1);

% ASD
[ fam1_0949.globaloccMap, fam1_0949.globalspikeMap, fam1_0949.globalinfoMap, fam1_0949.globalplaceMap, fam1_0949.globaldownData,...
    fam1_0949.globalactiveData, fam1_0949.globalplaceMap_smooth ] = ...
    generatePFmap( fam1_0949.globalspikes, [], fam1_0949.trackdata, 30.9, 10, 1, 1, 180, 1, 5);
[ fam1_0949.globalsorted_placeMap, fam1_0949.globalnormsorted_placeMap, fam1_0949.globalsortIdx ] = ...
    sortPFmap( fam1_0949.globalplaceMap, fam1_0949.globalinfoMap, 1);

Nbins = 180;
figure('Position',[1087 648 500 800]);
    subplot(9,5,2:5); imagesc(fam1_0949.globaloccMap);
        % xticks(Nbins/6:Nbins/6:Nbins); yticks([]);
        xticks([]); yticks([]);
        title('Occupancy map'); % colorbar;
    subplot(9,5,[6,11,16,21]); imagesc(fam1_0949.globalinfoMap(fam1_0949.globalsortIdx,2));
        xticks([]); yticks([]); ylabel('Cells'); 
        title('Info map'); colorbar;
    subplot(9,5,[7:10,12:15,17:20,22:25]);
        nspikeMap = fam1_0949.globalspikeMap./max(max(fam1_0949.globalspikeMap));
        imagesc(nspikeMap(fam1_0949.globalsortIdx,:)); xticks([]); yticks([]);
        title('Spike maps');
    subplot(9,5,[27:30,32:35,37:40,42:45]);    
        imagesc(fam1_0949.globalsorted_placeMap); 
        xticks(Nbins/6:Nbins/6:Nbins); yticks([]);
        xticklabels(360/Nbins*(Nbins/6:Nbins/6:Nbins));
        colorbar('Location','southoutside'); caxis([0,0.01]);
        title('Place field maps');
        xlabel('Position (degrees)'); ylabel('Cells');
clear Nbins nspikeMap


%% Combine data for each environment
fam1.tsG = [fam1_0944.globaltsG, fam1_0949.globaltsG];
fam1.tsR = [fam1_0944.globaltsR, fam1_0949.globaltsR];
fam1.R = [fam1_0944.globalR, fam1_0949.globalR];
fam1.spikes = [fam1_0944.globalspikes, fam1_0949.globalspikes];

fam1.trackdata.x = [fam1_0944.trackdata.x; fam1_0949.trackdata.x];
fam1.trackdata.y = [fam1_0944.trackdata.y; fam1_0949.trackdata.y];
fam1.trackdata.r = [fam1_0944.trackdata.r; fam1_0949.trackdata.r];
fam1.trackdata.phi = [fam1_0944.trackdata.phi; fam1_0949.trackdata.phi];
fam1.trackdata.speed = [fam1_0944.trackdata.speed; fam1_0949.trackdata.speed];
fam1.trackdata.time = [fam1_0944.trackdata.time; fam1_0944.trackdata.time(end)+fam1_0949.trackdata.time];

% histogram
[ fam1.occMap, fam1.spikeMap, fam1.infoMap, fam1.placeMap, fam1.downData,...
    fam1.activeData, fam1.placeMap_smooth ] = ...
    generatePFmap( fam1.spikes, [], fam1.trackdata, 30.9, 10, 1, 2, 180, 1, 5);
[ fam1.sorted_placeMap, fam1.normsorted_placeMap, fam1.sortIdx ] = ...
    sortPFmap( fam1.placeMap_smooth, fam1.infoMap, 1);

% ASD
% [ fam1.occMap, fam1.spikeMap, fam1.infoMap, fam1.placeMap, fam1.downData,...
%     fam1.activeData, fam1.placeMap_smooth ] = ...
%     generatePFmap( fam1.spikes, [], fam1.trackdata, 30.9, 10, 1, 1, 180, 1, 5);
% [ fam1.sorted_placeMap, fam1.normsorted_placeMap, fam1.sortIdx ] = ...
%     sortPFmap( fam1.placeMap, fam1.infoMap, 1);


Nbins = 180;
figure('Position',[1087 648 500 800]);
    subplot(9,5,2:5); imagesc(fam1.occMap);
        % xticks(Nbins/6:Nbins/6:Nbins); yticks([]);
        xticks([]); yticks([]);
        title('Occupancy map'); % colorbar;
    subplot(9,5,[6,11,16,21]); imagesc(fam1.infoMap(fam1.sortIdx,2));
        xticks([]); yticks([]); ylabel('Cells'); 
        title('Info map'); colorbar;
    subplot(9,5,[7:10,12:15,17:20,22:25]);
        nspikeMap = fam1.spikeMap./max(max(fam1.spikeMap));
        imagesc(nspikeMap(fam1.sortIdx,:)); xticks([]); yticks([]);
        title('Spike maps');
    subplot(9,5,[27:30,32:35,37:40,42:45]);    
        imagesc(fam1.sorted_placeMap); 
        xticks(Nbins/6:Nbins/6:Nbins); yticks([]);
        xticklabels(360/Nbins*(Nbins/6:Nbins/6:Nbins));
        colorbar('Location','southoutside'); caxis([0,0.01]);
        title('Place field maps');
        xlabel('Position (degrees)'); ylabel('Cells');
clear Nbins nspikeMap


%% Asses remapping
figure; Nbins = 180;
subplot(131); imagesc(fam1.sorted_placeMap); 
        xticks(Nbins/6:Nbins/6:Nbins); yticks([]);
        xticklabels(360/Nbins*(Nbins/6:Nbins/6:Nbins));
        colorbar('Location','southoutside'); 
        caxis([0,0.009]);
        title('Fam1 (9:44 & 9:49)');
        xlabel('Position (degrees)'); ylabel('Cells');
subplot(132); imagesc(fam1rev_1007.globalsorted_placeMap(fam1.sortIdx,:)); 
        xticks(Nbins/6:Nbins/6:Nbins); yticks([]);
        xticklabels(360/Nbins*(Nbins/6:Nbins/6:Nbins));
        colorbar('Location','southoutside'); 
        caxis([0,0.009]);
        title('Fam1rev (10:07)'); 
        xlabel('Position (degrees)'); ylabel('Cells');
subplot(133); imagesc(fam1rev_1007.globalsorted_placeMap); 
        xticks(Nbins/6:Nbins/6:Nbins); yticks([]);
        xticklabels(360/Nbins*(Nbins/6:Nbins/6:Nbins));
        colorbar('Location','southoutside'); 
        caxis([0,0.009]);
        title('Fam1rev (10:07) resorted'); 
        xlabel('Position (degrees)'); ylabel('Cells');

figure; Nbins = 180;
subplot(131); imagesc(fam1_0944.globalsorted_placeMap); 
        xticks(Nbins/6:Nbins/6:Nbins); yticks([]);
        xticklabels(360/Nbins*(Nbins/6:Nbins/6:Nbins));
        colorbar('Location','southoutside'); caxis([0,0.009]);
        title('Fam1 09:44');
        xlabel('Position (degrees)'); ylabel('Cells');
subplot(132); imagesc(fam1rev_1007.globalsorted_placeMap(fam1_0944.globalsortIdx,:)); 
        xticks(Nbins/6:Nbins/6:Nbins); yticks([]);
        xticklabels(360/Nbins*(Nbins/6:Nbins/6:Nbins));
        colorbar('Location','southoutside'); caxis([0,0.009]);
        title('Fam1rev 10:07'); 
        xlabel('Position (degrees)'); ylabel('Cells');
subplot(133); imagesc(fam1rev_1007.globalsorted_placeMap); 
        xticks(Nbins/6:Nbins/6:Nbins); yticks([]);
        xticklabels(360/Nbins*(Nbins/6:Nbins/6:Nbins));
        colorbar('Location','southoutside'); caxis([0,0.009]);
        title('Fam1rev 10:07 resorted'); 
        xlabel('Position (degrees)'); ylabel('Cells');

figure; Nbins = 180;
subplot(131); imagesc(fam1_0949.globalsorted_placeMap); 
        xticks(Nbins/6:Nbins/6:Nbins); yticks([]);
        xticklabels(360/Nbins*(Nbins/6:Nbins/6:Nbins));
        colorbar('Location','southoutside'); caxis([0,0.009]);
        title('Fam1 09:49');
        xlabel('Position (degrees)'); ylabel('Cells');
subplot(132); imagesc(fam1rev_1007.globalsorted_placeMap(fam1_0949.globalsortIdx,:)); 
        xticks(Nbins/6:Nbins/6:Nbins); yticks([]);
        xticklabels(360/Nbins*(Nbins/6:Nbins/6:Nbins));
        colorbar('Location','southoutside'); caxis([0,0.009]);
        title('Fam1rev 10:07'); 
        xlabel('Position (degrees)'); ylabel('Cells');
subplot(133); imagesc(fam1rev_1007.globalsorted_placeMap); 
        xticks(Nbins/6:Nbins/6:Nbins); yticks([]);
        xticklabels(360/Nbins*(Nbins/6:Nbins/6:Nbins));
        colorbar('Location','southoutside'); caxis([0,0.009]);
        title('Fam1rev 10:07 resorted'); 
        xlabel('Position (degrees)'); ylabel('Cells');

%% overlay 2 images
% A = fam1_0944.mean_imratio;
% B = fam1_0944.globalreg_mean_imratio;
% C = imfuse(A./max(max(A)),B./max(max(B)),'falsecolor','Scaling','joint','ColorChannels',[1 2 0]);
% figure; imshow(C);





