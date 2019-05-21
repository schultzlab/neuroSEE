data_locn = '/Volumes/thefarm2/live/CrazyEights/AD_2PCa/';
if ~exist(data_locn,'dir')
    data_locn = '/Volumes/RDS/project/thefarm2/live/CrazyEights/AD_2PCa/';
end
if ~exist(data_locn,'dir')
    data_locn = '/rds/general/user/mgo/projects/thefarm2/live/CrazyEights/AD_2PCa/';
end

filedir = [data_locn 'Data/20190406/Processed/'];

addpath(genpath('../behaviour'));
addpath(genpath('../intervideo_processing'));
addpath(genpath('../motion_correction'));
addpath(genpath('../PF_mapping'));
addpath(genpath('../ROI_segmentation'));
addpath(genpath('../spike_extraction'));
addpath(genpath('../utilities'));
addpath(genpath('../pipelines'));


% load('openfield_20190406_m82.mat');
% Nmasks = size(globalmasks,2);


% Load data for each recording
open_1816.green = green;
open_1816.red = red;
clear green red

% open_1816.tsG = cell_tsG;
% open_1816.tsR = cell_tsR;
% open_1816.masks = masks;
% open_1816.mean_imratio = mean_imratio;
% open_1816.R = R;
% open_1816.spikes = spikes;
% clear cell_tsG cell_tsR masks mean_imratio R spikes params

open_1816.trackdata.x = x;
open_1816.trackdata.y = y;
open_1816.trackdata.r = r;
open_1816.trackdata.phi = phi;
open_1816.trackdata.speed = speed;
open_1816.trackdata.time = time;
clear x y r phi speed TTLout w alpha time params


%% Compare red images to see which ones are most similar
% figure;
% subplot(241);imagesc(open_2014.red.meanregframe); title('Open 20:14'); caxis([3000 35000]); axis off; axis square;
% subplot(242);imagesc(open_2027.red.meanregframe); title('Open 20:27'); caxis([3000 35000]); axis off; axis square;
% subplot(243);imagesc(open_2033.red.meanregframe); title('Open 20:33'); caxis([3000 35000]); axis off; axis square;
% subplot(244);imagesc(open_2038.red.meanregframe); title('Open 20:38'); caxis([3000 35000]); axis off; axis square;
% subplot(245);imagesc(open_2041.red.meanregframe); title('Open 20:41'); caxis([3000 35000]); axis off; axis square;
% subplot(246);imagesc(open_2045.red.meanregframe); title('Open 20:45'); caxis([3000 35000]); axis off; axis square;
% subplot(247);imagesc(open_2051.red.meanregframe); title('Open 20:51'); caxis([3000 35000]); axis off; axis square;
% subplot(248);imagesc(open_2056.red.meanregframe); title('Open 20:56'); caxis([3000 35000]); axis off; axis square;


%% Compare trajectories
% figure; 
% subplot(251); plot(open_2014.trackdata.x,open_2014.trackdata.y); axis off; title('Open 20:14');
% subplot(252); plot(open_2027.trackdata.x,open_2027.trackdata.y); axis off; title('Open 20:27');
% subplot(253); plot(open_2033.trackdata.x,open_2033.trackdata.y); axis off; title('Open 20:33');
% subplot(254); plot(open_2038.trackdata.x,open_2038.trackdata.y); axis off; title('Open 20:38');
% subplot(255); 
%     plot(open_2014.trackdata.x,open_2014.trackdata.y); hold on
%     plot(open_2033.trackdata.x,open_2033.trackdata.y);
%     plot(open_2027.trackdata.x,open_2027.trackdata.y);
%     plot(open_2038.trackdata.x,open_2038.trackdata.y);
%     hold off
%     axis off; title('20:14, 20:27, 20:33, 20:38'); 
% subplot(256); plot(open_2041.trackdata.x,open_2041.trackdata.y); axis off; title('Open 20:41');
% subplot(257); plot(open_2045.trackdata.x,open_2045.trackdata.y); axis off; title('Open 20:45');
% subplot(258); plot(open_2051.trackdata.x,open_2051.trackdata.y); axis off; title('Open 20:51');
% subplot(259); plot(open_2056.trackdata.x,open_2056.trackdata.y); axis off; title('Open 20:56');
% subplot(2,5,10); 
%     plot(open_2041.trackdata.x,open_2041.trackdata.y); hold on
%     plot(open_2045.trackdata.x,open_2045.trackdata.y);
%     plot(open_2051.trackdata.x,open_2051.trackdata.y);
%     plot(open_2056.trackdata.x,open_2056.trackdata.y);
%         plot(open_2014.trackdata.x,open_2014.trackdata.y); 
%         plot(open_2033.trackdata.x,open_2033.trackdata.y);
%         plot(open_2027.trackdata.x,open_2027.trackdata.y);
%         plot(open_2038.trackdata.x,open_2038.trackdata.y);
%     hold off
%     axis off; title('all'); %title('20:41, 20:45, 20:51, 20:56');

        
%% Register images to the first one and calculate tsG, tsR, R, spikes
% reference image
% cprintf('Processing open_2014')
% open_2014_imG = read_file( [filedir '20190406_20_14_42/20190406_20_14_42_2P_XYT_green_mcorr.tif'] );
% open_2014_imR = read_file( [filedir '20190406_20_14_42/20190406_20_14_42_2P_XYT_red_mcorr.tif'] );
% 
% for i = 1:Nmasks
%     maskind = globalmasks{i};
%     for j = 1:size(open_2014_imG,3)
%         imG_reshaped = reshape( open_2014_imG(:,:,j), 512*512, 1);
%         open_2014.globaltsG( i, j ) = mean( imG_reshaped(maskind) );
%         imR_reshaped = reshape( open_2014_imR(:,:,j), 512*512, 1);
%         open_2014.globaltsR( i, j ) = mean( imR_reshaped(maskind) );
%     end
% end
% clear open_2014_imG open_2014_imR
% 
% open_2014.globalR = ratiometric_Ca( open_2014.globaltsG, open_2014.globaltsR, 11 );
% open_2014.globalspikes = nndORoasis(open_2014.globalR, 2, 0.94, 2.4);
% 
% % open_2027
% cprintf('Processing open_2027')
% [ open_2027.shift, open_2027.red.globalreg_meanframe ] = globalregisterImage(open_2014.red.meanregframe, open_2027.red.meanregframe, 1 );
% open_2027_imG = read_file( [filedir '20190406_20_27_07/20190406_20_27_07_2P_XYT_green_mcorr.tif'] );
% open_2027_imR = read_file( [filedir '20190406_20_27_07/20190406_20_27_07_2P_XYT_red_mcorr.tif'] );
% 
% [ open_2027_imGglobalreg, open_2027_imRglobalreg ] = ...
%     globalregisterStack( open_2027_imG, open_2027_imR, open_2027.shift(1) , open_2027.shift(2) );
% 
% clear open_2027_imG open_2027_imR
% 
% for i = 1:Nmasks
%     maskind = globalmasks{i};
%     for j = 1:size(open_2027_imGglobalreg,3)
%         imG_reshaped = reshape( open_2027_imGglobalreg(:,:,j), 512*512, 1);
%         open_2027.globaltsG( i, j ) = mean( imG_reshaped(maskind) );
%         imR_reshaped = reshape( open_2027_imRglobalreg(:,:,j), 512*512, 1);
%         open_2027.globaltsR( i, j ) = mean( imR_reshaped(maskind) );
%     end
% end
% clear open_2027_imGglobalreg open_2027_imRglobalreg
% 
% open_2027.globalR = ratiometric_Ca( open_2027.globaltsG, open_2027.globaltsR, 11 );
% open_2027.globalspikes = nndORoasis(open_2027.globalR, 2, 0.94, 2.4);
% 
% % open_2033
% cprintf('Processing open_2033')
% [ open_2033.shift, open_2033.red.globalreg_meanframe ] = globalregisterImage(open_2014.red.meanregframe, open_2033.red.meanregframe, 1 );
% open_2033_imG = read_file( [filedir '20190406_20_33_01/20190406_20_33_01_2P_XYT_green_mcorr.tif'] );
% open_2033_imR = read_file( [filedir '20190406_20_33_01/20190406_20_33_01_2P_XYT_red_mcorr.tif'] );
% 
% [ open_2033_imGglobalreg, open_2033_imRglobalreg ] = ...
%     globalregisterStack( open_2033_imG, open_2033_imR, open_2033.shift(1) , open_2033.shift(2) );
% 
% clear open_2033_imG open_2033_imR
% 
% for i = 1:Nmasks
%     maskind = globalmasks{i};
%     for j = 1:size(open_2033_imGglobalreg,3)
%         imG_reshaped = reshape( open_2033_imGglobalreg(:,:,j), 512*512, 1);
%         open_2033.globaltsG( i, j ) = mean( imG_reshaped(maskind) );
%         imR_reshaped = reshape( open_2033_imRglobalreg(:,:,j), 512*512, 1);
%         open_2033.globaltsR( i, j ) = mean( imR_reshaped(maskind) );
%     end
% end
% open_2033.globalR = ratiometric_Ca( open_2033.globaltsG, open_2033.globaltsR, 11 );
% open_2033.globalspikes = nndORoasis(open_2033.globalR, 2, 0.94, 2.4);
% 
% clear open_2033_imGglobalreg open_2033_imRglobalreg
% 
% % open_2038
% cprintf('Processing open_2038')
% [ open_2038.shift, open_2038.red.globalreg_meanframe ] = globalregisterImage(open_2014.red.meanregframe, open_2038.red.meanregframe, 1 );
% open_2038_imG = read_file( [filedir '20190406_20_38_41/20190406_20_38_41_2P_XYT_green_mcorr.tif'] );
% open_2038_imR = read_file( [filedir '20190406_20_38_41/20190406_20_38_41_2P_XYT_red_mcorr.tif'] );
% 
% [ open_2038_imGglobalreg, open_2038_imRglobalreg ] = ...
%     globalregisterStack( open_2038_imG, open_2038_imR, open_2038.shift(1) , open_2038.shift(2) );
% 
% clear open_2038_imG open_2038_imR
% 
% for i = 1:Nmasks
%     maskind = globalmasks{i};
%     for j = 1:size(open_2038_imGglobalreg,3)
%         imG_reshaped = reshape( open_2038_imGglobalreg(:,:,j), 512*512, 1);
%         open_2038.globaltsG( i, j ) = mean( imG_reshaped(maskind) );
%         imR_reshaped = reshape( open_2038_imRglobalreg(:,:,j), 512*512, 1);
%         open_2038.globaltsR( i, j ) = mean( imR_reshaped(maskind) );
%     end
% end
% open_2038.globalR = ratiometric_Ca( open_2038.globaltsG, open_2038.globaltsR, 11 );
% open_2038.globalspikes = nndORoasis(open_2038.globalR, 2, 0.94, 2.4);
% 
% clear open_2038_imGglobalreg open_2038_imRglobalreg
% 
% % open_2041
% cprintf('Processing open_2041')
% [ open_2041.shift, open_2041.red.globalreg_meanframe ] = globalregisterImage(open_2014.red.meanregframe, open_2041.red.meanregframe, 1 );
% open_2041_imG = read_file( [filedir '20190406_20_41_07/20190406_20_41_07_2P_XYT_green_mcorr.tif'] );
% open_2041_imR = read_file( [filedir '20190406_20_41_07/20190406_20_41_07_2P_XYT_red_mcorr.tif'] );
% 
% [ open_2041_imGglobalreg, open_2041_imRglobalreg ] = ...
%     globalregisterStack( open_2041_imG, open_2041_imR, open_2041.shift(1) , open_2041.shift(2) );
% 
% clear open_2041_imG open_2041_imR
% 
% for i = 1:Nmasks
%     maskind = globalmasks{i};
%     for j = 1:size(open_2041_imGglobalreg,3)
%         imG_reshaped = reshape( open_2041_imGglobalreg(:,:,j), 512*512, 1);
%         open_2041.globaltsG( i, j ) = mean( imG_reshaped(maskind) );
%         imR_reshaped = reshape( open_2041_imRglobalreg(:,:,j), 512*512, 1);
%         open_2041.globaltsR( i, j ) = mean( imR_reshaped(maskind) );
%     end
% end
% open_2041.globalR = ratiometric_Ca( open_2041.globaltsG, open_2041.globaltsR, 11 );
% open_2041.globalspikes = nndORoasis(open_2041.globalR, 2, 0.94, 2.4);
% 
% clear open_2041_imGglobalreg open_2041_imRglobalreg
% 
% % open_2045
% cprintf('Processing open_2045')
% [ open_2045.shift, open_2045.red.globalreg_meanframe ] = globalregisterImage(open_2014.red.meanregframe, open_2045.red.meanregframe, 1 );
% open_2045_imG = read_file( [filedir '20190406_20_45_35/20190406_20_45_35_2P_XYT_green_mcorr.tif'] );
% open_2045_imR = read_file( [filedir '20190406_20_45_35/20190406_20_45_35_2P_XYT_red_mcorr.tif'] );
% 
% [ open_2045_imGglobalreg, open_2045_imRglobalreg ] = ...
%     globalregisterStack( open_2045_imG, open_2045_imR, open_2045.shift(1) , open_2045.shift(2) );
% 
% clear open_2045_imG open_2045_imR
% 
% for i = 1:Nmasks
%     maskind = globalmasks{i};
%     for j = 1:size(open_2045_imGglobalreg,3)
%         imG_reshaped = reshape( open_2045_imGglobalreg(:,:,j), 512*512, 1);
%         open_2045.globaltsG( i, j ) = mean( imG_reshaped(maskind) );
%         imR_reshaped = reshape( open_2045_imRglobalreg(:,:,j), 512*512, 1);
%         open_2045.globaltsR( i, j ) = mean( imR_reshaped(maskind) );
%     end
% end
% open_2045.globalR = ratiometric_Ca( open_2045.globaltsG, open_2045.globaltsR, 11 );
% open_2045.globalspikes = nndORoasis(open_2045.globalR, 2, 0.94, 2.4);
% 
% clear open_2045_imGglobalreg open_2045_imRglobalreg
% 
% % open_2051
% cprintf('Processing open_2051')
% [ open_2051.shift, open_2051.red.globalreg_meanframe ] = globalregisterImage(open_2014.red.meanregframe, open_2051.red.meanregframe, 1 );
% open_2051_imG = read_file( [filedir '20190406_20_51_29/20190406_20_51_29_2P_XYT_green_mcorr.tif'] );
% open_2051_imR = read_file( [filedir '20190406_20_51_29/20190406_20_51_29_2P_XYT_red_mcorr.tif'] );
% 
% [ open_2051_imGglobalreg, open_2051_imRglobalreg ] = ...
%     globalregisterStack( open_2051_imG, open_2051_imR, open_2051.shift(1) , open_2051.shift(2) );
% 
% clear open_2051_imG open_2051_imR
% 
% for i = 1:Nmasks
%     maskind = globalmasks{i};
%     for j = 1:size(open_2051_imGglobalreg,3)
%         imG_reshaped = reshape( open_2051_imGglobalreg(:,:,j), 512*512, 1);
%         open_2051.globaltsG( i, j ) = mean( imG_reshaped(maskind) );
%         imR_reshaped = reshape( open_2051_imRglobalreg(:,:,j), 512*512, 1);
%         open_2051.globaltsR( i, j ) = mean( imR_reshaped(maskind) );
%     end
% end
% open_2051.globalR = ratiometric_Ca( open_2051.globaltsG, open_2051.globaltsR, 11 );
% open_2051.globalspikes = nndORoasis(open_2051.globalR, 2, 0.94, 2.4);
% 
% clear open_2051_imGglobalreg open_2051_imRglobalreg
% 
% % open_2056
% cprintf('Processing open_2056')
% [ open_2056.shift, open_2056.red.globalreg_meanframe ] = globalregisterImage(open_2014.red.meanregframe, open_2056.red.meanregframe, 1 );
% open_2056_imG = read_file( [filedir '20190406_20_56_30/20190406_20_56_30_2P_XYT_green_mcorr.tif'] );
% open_2056_imR = read_file( [filedir '20190406_20_56_30/20190406_20_56_30_2P_XYT_red_mcorr.tif'] );
% 
% [ open_2056_imGglobalreg, open_2056_imRglobalreg ] = ...
%     globalregisterStack( open_2056_imG, open_2056_imR, open_2056.shift(1) , open_2056.shift(2) );
% 
% clear open_2056_imG open_2056_imR
% 
% for i = 1:Nmasks
%     maskind = globalmasks{i};
%     for j = 1:size(open_2056_imGglobalreg,3)
%         imG_reshaped = reshape( open_2056_imGglobalreg(:,:,j), 512*512, 1);
%         open_2056.globaltsG( i, j ) = mean( imG_reshaped(maskind) );
%         imR_reshaped = reshape( open_2056_imRglobalreg(:,:,j), 512*512, 1);
%         open_2056.globaltsR( i, j ) = mean( imR_reshaped(maskind) );
%     end
% end
% open_2056.globalR = ratiometric_Ca( open_2056.globaltsG, open_2056.globaltsR, 11 );
% open_2056.globalspikes = nndORoasis(open_2056.globalR, 2, 0.94, 2.4);
% 
% clear open_2056_imGglobalreg open_2056_imRglobalreg
% clear i j maskind imG_reshaped imR_reshaped
% 
% 
% %% Combine data for each environment
% cprintf('Combining data for open field')
% open.tsG    = [open_2014.globaltsG,    open_2027.globaltsG,    open_2033.globaltsG,    open_2038.globaltsG,    open_2041.globaltsG,    open_2045.globaltsG,    open_2051.globaltsG,    open_2056.globaltsG   ];
% open.tsR    = [open_2014.globaltsR,    open_2027.globaltsR,    open_2033.globaltsR,    open_2038.globaltsR,    open_2041.globaltsR,    open_2045.globaltsR,    open_2051.globaltsR,    open_2056.globaltsR   ];
% open.R      = [open_2014.globalR,      open_2027.globalR,      open_2033.globalR,      open_2038.globalR,      open_2041.globalR,      open_2045.globalR,      open_2051.globalR,      open_2056.globalR     ];
% open.spikes = [open_2014.globalspikes, open_2027.globalspikes, open_2033.globalspikes, open_2038.globalspikes, open_2041.globalspikes, open_2045.globalspikes, open_2051.globalspikes, open_2056.globalspikes];
% 
% open.trackdata.x     = [open_2014.trackdata.x,     open_2027.trackdata.x,     open_2033.trackdata.x,     open_2038.trackdata.x,     open_2041.trackdata.x,     open_2045.trackdata.x,     open_2051.trackdata.x,     open_2056.trackdata.x    ];
% open.trackdata.y     = [open_2014.trackdata.y,     open_2027.trackdata.y,     open_2033.trackdata.y,     open_2038.trackdata.y,     open_2041.trackdata.y,     open_2045.trackdata.y,     open_2051.trackdata.y,     open_2056.trackdata.y    ];
% open.trackdata.r     = [open_2014.trackdata.r,     open_2027.trackdata.r,     open_2033.trackdata.r,     open_2038.trackdata.r,     open_2041.trackdata.r,     open_2045.trackdata.r,     open_2051.trackdata.r,     open_2056.trackdata.r    ];
% open.trackdata.phi   = [open_2014.trackdata.phi,   open_2027.trackdata.phi,   open_2033.trackdata.phi,   open_2038.trackdata.phi,   open_2041.trackdata.phi,   open_2045.trackdata.phi,   open_2051.trackdata.phi,   open_2056.trackdata.phi  ];
% open.trackdata.speed = [open_2014.trackdata.speed, open_2027.trackdata.speed, open_2033.trackdata.speed, open_2038.trackdata.speed, open_2041.trackdata.speed, open_2045.trackdata.speed, open_2051.trackdata.speed, open_2056.trackdata.speed];
% tmax = open_2014.trackdata.time(end);
% open.trackdata.time  = [open_2014.trackdata.time,  tmax+open_2027.trackdata.time, 2*tmax+open_2033.trackdata.time, 3*tmax+open_2038.trackdata.time, 4*tmax+open_2041.trackdata.time, 5*tmax+open_2045.trackdata.time, 6*tmax+open_2051.trackdata.time, 7*tmax+open_2056.trackdata.time];
% 
% % downsample data
% tracktime   = open.trackdata.time;
% x           = open.trackdata.x;
% y           = open.trackdata.y;
% r           = open.trackdata.r;
% phi         = open.trackdata.phi;
% speed       = open.trackdata.speed;
% spikes      = open.spikes;
% 
% t0 = tracktime(1);
% Nt = size(spikes,2);
% dt = 1/30.91;
% t = (t0:dt:Nt*dt)';
% Vthr = 10;
% 
% % Downsample tracking to Ca trace
% downphi     = interp1(tracktime,phi,t,'linear');
% downx       = interp1(tracktime,x,t,'linear');
% downy       = interp1(tracktime,y,t,'linear');
% downspeed   = interp1(tracktime,speed,t,'linear'); % mm/s
% downr       = interp1(tracktime,r,t,'linear'); % mm/s
% 
% % Consider only samples when the mouse is active
% activex     = downx(downspeed > Vthr);
% activey     = downy(downspeed > Vthr);
% activephi   = downphi(downspeed > Vthr);
% activespeed = speed(downspeed > Vthr);
% activer     = r(downspeed > Vthr);
% activespikes = spikes(:,downspeed > Vthr);
% activet     = t(downspeed > Vthr);
% 
% % save
% open.downData.x = downx;
% open.downData.y = downy;
% open.downData.phi = downphi;
% open.downData.speed = downspeed;
% open.downData.t = t;
% open.downData.r = downr;
% 
% open.activeData.x = activex;
% open.activeData.y = activey;
% open.activeData.phi = activephi;
% open.activeData.speed = activespeed;
% open.activeData.t = activet;
% open.activeData.spikes = activespikes;
% open.activeData.r = activer;
% 
% clear activephi activer activespeed activespikes activet activex activey ans
% clear downphi downr downspeed downx downy dt Nt phi r speed spikes t t0 tracktime Vthr x y
% save 'openfield_20190406_m82.mat'

%% Generate PF maps with ASD
cprintf('Generating PF maps')
spikes = open.activeData.spikes;
xp1 = open.activeData.x;
xp2 = open.activeData.y;

% addpath(genpath([pwd,'/fastASD'])); % add ASD folder

npart       = 1;                 % number of epochs to divide trial in
ASD         = 1;                 % perform ASD estimation
track_plot  = 0;                 % display tracking and firing for sample cells
fsamp       = 30;                % in Hz
n1=100; n2=100; nks=[n1,n2];     % env discretisation for ASD estimation
h1 = 20; h2 = 20; hs = [h1,h2];
%%%%%%%%%%%%
% spikes = activeData.spikes; % extract spiking matrix
%%%%%%%%%%%%
% initialise  variables to store results
kasd    = zeros(npart, n1, n2, size(spikes,1));
hpf     = zeros(npart, h1, h2, size(spikes,1));
hpfs    = zeros(npart, h1, h2, size(spikes,1));
infoH   = zeros(npart, size(spikes,1), 2);
infoASD = zeros(npart, size(spikes,1), 2);

% process tracking-environment data
%%%%%%%%%%%%%
% xp1 = activeData.x;
% xp2 = activeData.y;
%%%%%%%%%%%%%
xp1 = (xp1-min(xp1))/(max(xp1)-min(xp1)); % normalised 0 mean
xp2 = (xp2-min(xp2))/(max(xp2)-min(xp2));
xp  = [xp1,xp2];

% discretize x and y position for histogram estimation
x1 = linspace(0,1.0001,h1+1);
x2 = linspace(0,1.0001,h2+1);
xh = floor((xp(:,1)-x1(1))/(x1(2)-x1(1)))+1;
yh = floor((xp(:,2)-x2(1))/(x2(2)-x2(1)))+1;
occMap = full(sparse(xh,yh,1,h1,h2));
mode = 0; % the mask is obtained by imfill only
envMask_h = getEnvEdgePrior(occMap,mode); % hist

% discretize x and y position for ASD estimation
x1 = linspace(0,1.0001,n1+1);
x2 = linspace(0,1.0001,n2+1);
xi = floor((xp(:,1)-x1(1))/(x1(2)-x1(1)))+1;
yi = floor((xp(:,2)-x2(1))/(x2(2)-x2(1)))+1;
xind = sub2ind([n1,n2],xi,yi); % flatten bin tracking (for ASD)
occMap = full(sparse(xi,yi,1,n1,n2));
mode = 2; % the mask is obtained by dilation and imfill
envMask_asd = getEnvEdgePrior(occMap,mode); % ASD

% find which neurons are spiking
for ii = 1:size(spikes,1)
    a(ii) = sum(spikes(ii,:));
end
spk_idx = find(a); % store indices

for id = 1:length(spk_idx)
    z = spikes(spk_idx(id),:);
    if track_plot % plot sample tracking & spikes
        figure; hold on;
        plot(xp1,xp2); plot(xp1(z>0),xp2(z>0),'r.','markersize',10);
    end
    % separate exploration in smaller intervals
    pp_idx = round(linspace(1,length(z),npart+1));
    for pp = 1:npart
            xpp = xind(pp_idx(pp):pp_idx(pp+1));
            zpp = z(pp_idx(pp):pp_idx(pp+1));
        % ASD estimation
        if ASD
            [asd,~] = runASD_2d(xpp',zpp',nks,envMask_asd);
             if min(asd)<0; asd = asd-min(asd); end
             kasd(pp,:,:,spk_idx(id)) = asd;
        end
        % histogram estimation
        occMap = full(sparse(xh,yh,1,h1,h2));
        counts = full(sparse(xh,yh,z,h1,h2));
        hhh = counts./occMap;
        hhh(isnan(hhh)) = 0;
        hpf(pp,:,:,spk_idx(id)) = hhh;
        hhh = imgaussfilt(hhh,h1/10); hhh(~envMask_h) = 0;
        hpfs(pp,:,:,spk_idx(id)) = hhh;
        % info estimation
        [infoASD(pp,spk_idx(id),1), infoASD(pp,spk_idx(id),2)] =...
            infoMeasures(asd',ones(n1,n2),0);
        [infoH(pp,spk_idx(id),1), infoH(pp,spk_idx(id),2)] = infoMeasures(hhh,occMap,0);
    end
end

%% plot PF estimate results

%n_info = 1/10; % ratio of units to consider as not informative
info_type = 2; % 1 is info/sec, 2 is info/spk\

%% ASD
% if ASD
%     for pp = 1:npart
%         for ii = 0:length(spk_idx)/16
%             figure
%             ha = tight_subplot(4,4,[.03 .005],[.01 .05],[.01 .01]);
%             for jj=1:16
%                 if (ii*16+jj)<=length(spk_idx)
%                     axes(ha(jj));
%                     imagesc(squeeze(kasd(pp,:,:,spk_idx(ii*16+jj)))');
%                     axis off; colorbar;
%                     title(['cell ',num2str(spk_idx(ii*16+jj))],'fontsize',15)
%                 end
%             end
%             ax = subtitle('ASD'); set(ax, 'fontsize',20);
%         end 
%     end
% end

%% HIST (not smoothed)
% for pp = 1:npart
%     for ii = 0:length(spk_idx)/16
%         figure
%         ha = tight_subplot(4,4,[.01 .005],[.01 .07],[.01 .01]);
%         for jj=1:16
%             if (ii*16+jj)<=length(spk_idx)
%                 axes(ha(jj));
%                 imagesc(squeeze(hpf(pp,:,:,spk_idx(ii*16+jj)))');
%                 axis off; colorbar;
%                 title(['cell ',num2str(spk_idx(ii*16+jj))],'fontsize',15)
%             end
%         end
%         ax = subtitle('RAW HIST'); set(ax, 'fontsize',20);
%     end 
% end

%% HIST (smoothed)
% for pp = 1:npart
%     for ii = 0:length(spk_idx)/16
%         figure
%         ha = tight_subplot(4,4,[.01 .005],[.01 .07],[.01 .01]);
%         for jj=1:16
%             if (ii*16+jj)<=length(spk_idx)
%                 axes(ha(jj));
%                 imagesc(squeeze(hpfs(pp,:,:,spk_idx(ii*16+jj)))');
%                 axis off; colorbar;
%                 title(['cell ',num2str(spk_idx(ii*16+jj))],'fontsize',15)
%             end
%         end
%         ax = subtitle('SMOOTH HIST'); set(ax, 'fontsize',20);
%     end 
% end

%% TRACKING
% for pp = 1:npart
%     for ii = 0:length(spk_idx)/16
%         figure
%         ha = tight_subplot(4,4,[.01 .005],[.01 .07],[.01 .01]);
%         for jj=1:16
%             if (ii*16+jj)<=length(spk_idx)
%                 axes(ha(jj));
%                 z = spikes(spk_idx(ii*16+jj),:);
%                 hold on; axis off;
%                 plot(xp1,-xp2); plot(xp1(z>0),-xp2(z>0),'r.','markersize',10);
%                 title(['cell ',num2str(spk_idx(ii*16+jj))],'fontsize',15)
%             end
%         end
%         ax = subtitle('TRACKING'); set(ax, 'fontsize',20);
%     end 
% end

%%  ALL RESULTS TOGETHER
nPlot = 4;
for pp = 1:npart
    for ii=0:length(spk_idx)/nPlot
        figure
        ha = tight_subplot(nPlot,4,[.01 .005],[.01 .07],[.01 .01]);
        for jj=0:3
            if (ii*nPlot+jj) <= length(spk_idx)
                axes(ha(jj*nPlot+1));
                z = spikes(spk_idx(ii*nPlot+jj+1),:);
                hold on; axis off;
                plot(xp1,-xp2); plot(xp1(z>0),-xp2(z>0),'r.','markersize',10);
                axes(ha(jj*nPlot+2));
                imagesc(squeeze(hpf(pp,:,:,spk_idx(ii*nPlot+jj+1)))');
                axis off; colorbar; %caxis([0 0.06]);
                title(['cell ',num2str(spk_idx(ii*nPlot+jj+1))],'fontsize',15)
                axes(ha(jj*nPlot+3)); 
                imagesc(squeeze(hpfs(pp,:,:,spk_idx(ii*nPlot+jj+1)))');
                axis off; colorbar; %caxis([0 0.005]);
                axes(ha(jj*nPlot+4));
                imagesc(squeeze(kasd(pp,:,:,spk_idx(ii*nPlot+jj+1)))');
                axis off; colorbar; %caxis([0 0.003]);
            end
        end
    end 
end

% % Plot selected cells
nPlot = 3;
% star_idx = [6,7,20,26,36,47,59,55,71,73,78,94,99,104,110,120,145,144]; % varying colorbars
star_idx = [151, 128, 101, 86, 70, 60]; % constant colorbar
for pp = 1:npart
    for ii=0:length(star_idx)/nPlot
        figure;
        ha = tight_subplot(6,3,[.01 .01],[.01 .07],[.01 .01]);
        for jj=0:5
            if (ii*nPlot+jj) <= length(star_idx)
                axes(ha(jj*nPlot+1));
                z1 = spikes(star_idx(ii*nPlot+jj+1),:);
                hold on; axis off;
                plot(xp1,-xp2); plot(xp1(z1>0),-xp2(z1>0),'r.','markersize',10);
                axes(ha(jj*nPlot+2)); colorbar;
                z2 = hpf(pp,:,:,star_idx(ii*nPlot+jj+1));
                imagesc(squeeze(z2)');
                axis off; colorbar; %caxis([0 0.01]);
                %title(['cell ',num2str(star_idx(ii*nPlot+jj+1))],'fontsize',15)
                axes(ha(jj*nPlot+3));
                imagesc(squeeze(hpfs(pp,:,:,star_idx(ii*nPlot+jj+1)))');
                axis off; colorbar; % caxis([0 0.007]);
%                 axes(ha(jj*nPlot+4));
%                 imagesc(squeeze(kasd(pp,:,:,star_idx(ii*nPlot+jj+1)))');
%                 axis off; colorbar; %caxis([0 0.005]);
            end
        end
    end 
end

figure; 
id = 151;
subplot(311); plot(open.tsG(id,:));
subplot(312); plot(open.R(id,:));
subplot(313); plot(open.spikes(id,:));


% save 'openfield_20190406_m82_ASD.mat'
