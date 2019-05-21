% data_locn = '/Volumes/thefarm2/live/CrazyEights/AD_2PCa/Data/20181016/Processed/';
data_locn = '/Volumes/RDS/project/thefarm2/live/CrazyEights/AD_2PCa/Data/20190404/Processed/';

open_1543.green = green;
open_1543.red = red;
clear green red

open_1543.tsG = cell_tsG;
open_1543.tsR = cell_tsR;
open_1543.masks = masks;
open_1543.mean_imratio = mean_imratio;
open_1543.R = R;
open_1543.spikes = spikes;
clear cell_tsG cell_tsR masks mean_imratio R spikes params


%% Compare red images to see which ones are most similar
figure;
subplot(241);imagesc(open_1449.red.meanregframe); title('Open 14:49'); caxis([1500 6500]); axis off; axis square;
subplot(242);imagesc(open_1501.red.meanregframe); title('Open 15:01'); caxis([1500 6500]); axis off; axis square;
subplot(243);imagesc(open_1506.red.meanregframe); title('Open 15:06'); caxis([1500 6500]); axis off; axis square;
subplot(244);imagesc(open_1512.red.meanregframe); title('Open 15:12'); caxis([1500 6500]); axis off; axis square;
subplot(245);imagesc(open_1521.red.meanregframe); title('Open 15:21'); caxis([1500 6500]); axis off; axis square;
subplot(246);imagesc(open_1528.red.meanregframe/1.5); title('Open 15:28'); caxis([1500 6500]); axis off; axis square;
subplot(247);imagesc(open_1535.red.meanregframe); title('Open 15:35'); caxis([1500 6500]); axis off; axis square;
subplot(248);imagesc(open_1543.red.meanregframe); title('Open 15:43'); caxis([1500 6500]); axis off; axis square;


%% Compare trajectories
open_1543.trackdata.x = x;
open_1543.trackdata.y = y;
open_1543.trackdata.r = r;
open_1543.trackdata.phi = phi;
open_1543.trackdata.speed = speed;
open_1543.trackdata.time = time;
clear x y r phi speed TTLout w alpha time

figure; 
subplot(251); plot(open_1449.trackdata.x,open_1449.trackdata.y); axis off; title('Open 14:49');
subplot(252); plot(open_1501.trackdata.x,open_1501.trackdata.y); axis off; title('Open 15:01');
subplot(253); plot(open_1506.trackdata.x,open_1506.trackdata.y); axis off; title('Open 15:06');
subplot(254); plot(open_1512.trackdata.x,open_1512.trackdata.y); axis off; title('Open 15:12');
subplot(255); 
    plot(open_1449.trackdata.x,open_1449.trackdata.y); hold on
    plot(open_1501.trackdata.x,open_1501.trackdata.y);
    plot(open_1506.trackdata.x,open_1506.trackdata.y);
    plot(open_1512.trackdata.x,open_1512.trackdata.y); 
    %plot(open_1521.trackdata.x,open_1521.trackdata.y); 
    hold off
    axis off; title('14:49, 15:01, 15:06, 15:12'); 
subplot(256); plot(open_1521.trackdata.x,open_1521.trackdata.y); axis off; title('Open 15:21');
subplot(257); plot(open_1528.trackdata.x,open_1528.trackdata.y); axis off; title('Open 15:28');
subplot(258); plot(open_1535.trackdata.x,open_1535.trackdata.y); axis off; title('Open 15:35');
subplot(259); plot(open_1543.trackdata.x,open_1543.trackdata.y); axis off; title('Open 15:43');
subplot(2,5,10); 
    plot(open_1521.trackdata.x,open_1521.trackdata.y); hold on
    plot(open_1528.trackdata.x,open_1528.trackdata.y);
    plot(open_1535.trackdata.x,open_1535.trackdata.y);
    plot(open_1512.trackdata.x,open_1512.trackdata.y); hold off
    axis off; title('15:21, 15:28, 15:35, 15:12');

        
%% Register images to the first one
open_1449_imG = read_file( [data_locn '20190404_14_49_51/20190404_14_49_51_2P_XYT_green_mcorr.tif'] );
open_1449_imR = read_file( [data_locn '20190404_14_49_51/20190404_14_49_51_2P_XYT_red_mcorr.tif'] );

% open_1501
[ open_1501.shift, open_1501.red.globalreg_meanframe ] = globalregisterImage(open_1449.red.meanregframe, open_1501.red.meanregframe, 1 );
open_1501_imG = read_file( [data_locn '20190404_15_01_29/20190404_15_01_29_2P_XYT_green_mcorr.tif'] );
open_1501_imR = read_file( [data_locn '20190404_15_01_29/20190404_15_01_29_2P_XYT_red_mcorr.tif'] );

[ open_1501_imGglobalreg, open_1501_imRglobalreg ] = ...
    globalregisterStack( open_1501_imG, open_1501_imR, open_1501.shift(1) , open_1501.shift(2) );

% open_1506
[ open_1506.shift, open_1506.red.globalreg_meanframe ] = globalregisterImage(open_1449.red.meanregframe, open_1506.red.meanregframe, 1 );
open_1506_imG = read_file( [data_locn '20190404_15_06_46/20190404_15_06_46_2P_XYT_green_mcorr.tif'] );
open_1506_imR = read_file( [data_locn '20190404_15_06_46/20190404_15_06_46_2P_XYT_red_mcorr.tif'] );

[ open_1506_imGglobalreg, open_1506_imRglobalreg ] = ...
    globalregisterStack( open_1506_imG, open_1506_imR, open_1506.shift(1) , open_1506.shift(2) );

% open_1512
[ open_1512.shift, open_1512.red.globalreg_meanframe ] = globalregisterImage(open_1449.red.meanregframe, open_1512.red.meanregframe, 1 );
open_1512_imG = read_file( [data_locn '20190404_15_12_38/20190404_15_12_38_2P_XYT_green_mcorr.tif'] );
open_1512_imR = read_file( [data_locn '20190404_15_12_38/20190404_15_12_38_2P_XYT_red_mcorr.tif'] );

[ open_1512_imGglobalreg, open_1512_imRglobalreg ] = ...
    globalregisterStack( open_1512_imG, open_1512_imR, open_1512.shift(1) , open_1512.shift(2) );

% open_1521
[ open_1521.shift, open_1521.red.globalreg_meanframe ] = globalregisterImage(open_1449.red.meanregframe, open_1521.red.meanregframe, 1 );
open_1521_imG = read_file( [data_locn '20190404_15_21_48/20190404_15_21_48_2P_XYT_green_mcorr.tif'] );
open_1521_imR = read_file( [data_locn '20190404_15_21_48/20190404_15_21_48_2P_XYT_red_mcorr.tif'] );

[ open_1521_imGglobalreg, open_1521_imRglobalreg ] = ...
    globalregisterStack( open_1521_imG, open_1521_imR, open_1521.shift(1) , open_1521.shift(2) );


%% Calculate tsG, tsR, R, spikes 
% open_1449
Nmasks = size(globalmasks,2);

for i = 1:Nmasks
    maskind = globalmasks{i};
    for j = 1:size(open_1449_imG,3)
        imG_reshaped = reshape( open_1449_imG(:,:,j), 512*512, 1);
        open_1449.globaltsG( i, j ) = mean( imG_reshaped(maskind) );
        imR_reshaped = reshape( open_1449_imR(:,:,j), 512*512, 1);
        open_1449.globaltsR( i, j ) = mean( imR_reshaped(maskind) );
    end
end
open_1449.globalR = ratiometric_Ca( open_1449.globaltsG, open_1449.globaltsR, 11 );
open_1449.globalspikes = nndORoasis(open_1449.globalR, 2, 0.94, 2.4);

% open_1501
for i = 1:Nmasks
    maskind = globalmasks{i};
    for j = 1:size(open_1501_imG,3)
        imG_reshaped = reshape( open_1501_imGglobalreg(:,:,j), 512*512, 1);
        open_1501.globaltsG( i, j ) = mean( imG_reshaped(maskind) );
        imR_reshaped = reshape( open_1501_imRglobalreg(:,:,j), 512*512, 1);
        open_1501.globaltsR( i, j ) = mean( imR_reshaped(maskind) );
    end
end
open_1501.globalR = ratiometric_Ca( open_1501.globaltsG, open_1501.globaltsR, 11 );
open_1501.globalspikes = nndORoasis(open_1501.globalR, 2, 0.94, 2.4);

% open_1506
for i = 1:Nmasks
    maskind = globalmasks{i};
    for j = 1:size(open_1506_imGglobalreg,3)
        imG_reshaped = reshape( open_1506_imGglobalreg(:,:,j), 512*512, 1);
        open_1506.globaltsG( i, j ) = mean( imG_reshaped(maskind) );
        imR_reshaped = reshape( open_1506_imRglobalreg(:,:,j), 512*512, 1);
        open_1506.globaltsR( i, j ) = mean( imR_reshaped(maskind) );
    end
end
open_1506.globalR = ratiometric_Ca( open_1506.globaltsG, open_1506.globaltsR, 11 );
open_1506.globalspikes = nndORoasis(open_1506.globalR, 2, 0.94, 2.4);

% open_1512
for i = 1:Nmasks
    maskind = globalmasks{i};
    for j = 1:size(open_1512_imG,3)
        imG_reshaped = reshape( open_1512_imGglobalreg(:,:,j), 512*512, 1);
        open_1512.globaltsG( i, j ) = mean( imG_reshaped(maskind) );
        imR_reshaped = reshape( open_1512_imRglobalreg(:,:,j), 512*512, 1);
        open_1512.globaltsR( i, j ) = mean( imR_reshaped(maskind) );
    end
end
open_1512.globalR = ratiometric_Ca( open_1512.globaltsG, open_1512.globaltsR, 11 );
open_1512.globalspikes = nndORoasis(open_1512.globalR, 2, 0.94, 2.4);

% open_1521
for i = 1:Nmasks
    maskind = globalmasks{i};
    for j = 1:size(open_1521_imG,3)
        imG_reshaped = reshape( open_1521_imGglobalreg(:,:,j), 512*512, 1);
        open_1521.globaltsG( i, j ) = mean( imG_reshaped(maskind) );
        imR_reshaped = reshape( open_1521_imRglobalreg(:,:,j), 512*512, 1);
        open_1521.globaltsR( i, j ) = mean( imR_reshaped(maskind) );
    end
end
open_1521.globalR = ratiometric_Ca( open_1521.globaltsG, open_1521.globaltsR, 11 );
open_1521.globalspikes = nndORoasis(open_1521.globalR, 2, 0.94, 2.4);


%% Combine data for each environment
open.tsG    = [open_1449.globaltsG,    open_1501.globaltsG,    open_1506.globaltsG,    open_1512.globaltsG,    open_1521.globaltsG   ];
open.tsR    = [open_1449.globaltsR,    open_1501.globaltsR,    open_1506.globaltsR,    open_1512.globaltsR,    open_1521.globaltsR   ];
open.R      = [open_1449.globalR,      open_1501.globalR,      open_1506.globalR,      open_1512.globalR,      open_1521.globalR     ];
open.spikes = [open_1449.globalspikes, open_1501.globalspikes, open_1506.globalspikes, open_1512.globalspikes, open_1521.globalspikes];

open.trackdata.x     = [open_1449.trackdata.x,     open_1501.trackdata.x,         open_1506.trackdata.x,           open_1512.trackdata.x,           open_1521.trackdata.x];
open.trackdata.y     = [open_1449.trackdata.y,     open_1501.trackdata.y,         open_1506.trackdata.y,           open_1512.trackdata.y,           open_1521.trackdata.y];
open.trackdata.r     = [open_1449.trackdata.r,     open_1501.trackdata.r,         open_1506.trackdata.r,           open_1512.trackdata.r,           open_1521.trackdata.r];
open.trackdata.phi   = [open_1449.trackdata.phi,   open_1501.trackdata.phi,      open_1506.trackdata.phi,         open_1512.trackdata.phi,         open_1521.trackdata.phi];
open.trackdata.speed = [open_1449.trackdata.speed, open_1501.trackdata.speed,     open_1506.trackdata.speed,       open_1512.trackdata.speed,       open_1521.trackdata.speed];
tmax = open_1449.trackdata.time(end);
open.trackdata.time  = [open_1449.trackdata.time,  tmax+open_1501.trackdata.time, 2*tmax+open_1506.trackdata.time, 3*tmax+open_1512.trackdata.time, 4*tmax+open_1521.trackdata.time];

% downsample data
tracktime   = open.trackdata.time;
x           = open.trackdata.x;
y           = open.trackdata.y;
r           = open.trackdata.r;
phi         = open.trackdata.phi;
speed       = open.trackdata.speed;
spikes      = open.spikes;

t0 = tracktime(1);
Nt = size(spikes,2);
dt = 1/30.91;
t = (t0:dt:Nt*dt)';
Vthr = 10;

% Downsample tracking to Ca trace
downphi     = interp1(tracktime,phi,t,'linear');
downx       = interp1(tracktime,x,t,'linear');
downy       = interp1(tracktime,y,t,'linear');
downspeed   = interp1(tracktime,speed,t,'linear'); % mm/s
downr       = interp1(tracktime,r,t,'linear'); % mm/s

% Consider only samples when the mouse is active
activex     = downx(downspeed > Vthr);
activey     = downy(downspeed > Vthr);
activephi   = downphi(downspeed > Vthr);
activespeed = speed(downspeed > Vthr);
activer     = r(downspeed > Vthr);
activespikes = spikes(:,downspeed > Vthr);
activet     = t(downspeed > Vthr);

% save
open.downData.x = downx;
open.downData.y = downy;
open.downData.phi = downphi;
open.downData.speed = downspeed;
open.downData.t = t;
open.downData.r = downr;

open.activeData.x = activex;
open.activeData.y = activey;
open.activeData.phi = activephi;
open.activeData.speed = activespeed;
open.activeData.t = activet;
open.activeData.spikes = activespikes;
open.activeData.r = activer;

clear activephi activer activespeed activespikes activet activex activey ans
clear downphi downr downspeed downx downy dt Nt phi r speed spikes t t0 tracktime Vthr x y


% histogram
[ fam1.occMap, fam1.spikeMap, fam1.infoMap, fam1.placeMap, fam1.downData,...
    fam1.activeData, fam1.placeMap_smooth ] = ...
    generatePFmap( fam1.spikes, [], fam1.trackdata, 30.9, 10, 1, 2, 180, 1, 5);
[ fam1.sorted_placeMap, fam1.normsorted_placeMap, fam1.sortIdx ] = ...
    sortPFmap( fam1.placeMap_smooth, fam1.infoMap, 1);

% ASD
[ fam1.occMap, fam1.spikeMap, fam1.infoMap, fam1.placeMap, fam1.downData,...
    fam1.activeData, fam1.placeMap_smooth ] = ...
    generatePFmap( fam1.spikes, [], fam1.trackdata, 30.9, 10, 1, 1, 180, 1, 5);
[ fam1.sorted_placeMap, fam1.normsorted_placeMap, fam1.sortIdx ] = ...
    sortPFmap( fam1.placeMap, fam1.infoMap, 1);


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


%% Assess remapping
figure; Nbins = 180;
subplot(131); imagesc(fam1.sorted_placeMap); 
        xticks(Nbins/6:Nbins/6:Nbins); yticks([]);
        xticklabels(360/Nbins*(Nbins/6:Nbins/6:Nbins));
        colorbar('Location','southoutside'); caxis([0,0.01]);
        title('Fam1 (9:44 & 9:49)');
        xlabel('Position (degrees)'); ylabel('Cells');
subplot(132); imagesc(fam1rev_1007.globalsorted_placeMap(fam1.sortIdx,:)); 
        xticks(Nbins/6:Nbins/6:Nbins); yticks([]);
        xticklabels(360/Nbins*(Nbins/6:Nbins/6:Nbins));
        colorbar('Location','southoutside'); caxis([0,0.01]);
        title('Fam1rev (10:07)'); 
        xlabel('Position (degrees)'); ylabel('Cells');
subplot(133); imagesc(fam1rev_1007.globalsorted_placeMap); 
        xticks(Nbins/6:Nbins/6:Nbins); yticks([]);
        xticklabels(360/Nbins*(Nbins/6:Nbins/6:Nbins));
        colorbar('Location','southoutside'); caxis([0,0.01]);
        title('Fam1rev (10:07) resorted'); 
        xlabel('Position (degrees)'); ylabel('Cells');

figure; Nbins = 180;
subplot(131); imagesc(fam1_0944.globalsorted_placeMap); 
        xticks(Nbins/6:Nbins/6:Nbins); yticks([]);
        xticklabels(360/Nbins*(Nbins/6:Nbins/6:Nbins));
        colorbar('Location','southoutside'); caxis([0,0.01]);
        title('Fam1 09:44');
        xlabel('Position (degrees)'); ylabel('Cells');
subplot(132); imagesc(fam1rev_1007.globalsorted_placeMap(fam1_0944.globalsortIdx,:)); 
        xticks(Nbins/6:Nbins/6:Nbins); yticks([]);
        xticklabels(360/Nbins*(Nbins/6:Nbins/6:Nbins));
        colorbar('Location','southoutside'); caxis([0,0.01]);
        title('Fam1rev 10:07'); 
        xlabel('Position (degrees)'); ylabel('Cells');
subplot(133); imagesc(fam1rev_1007.globalsorted_placeMap); 
        xticks(Nbins/6:Nbins/6:Nbins); yticks([]);
        xticklabels(360/Nbins*(Nbins/6:Nbins/6:Nbins));
        colorbar('Location','southoutside'); caxis([0,0.01]);
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
        title('Fam1rev 09:49'); 
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





