clear; close all;
tic

%% Basic setup
addpath(genpath('../behaviour'));
addpath(genpath('../intervideo_processing'));
addpath(genpath('../motion_correction'));
addpath(genpath('../PF_mapping'));
addpath(genpath('../pipelines'));
addpath(genpath('../ROI_segmentation'));
addpath(genpath('../spike_extraction'));
addpath(genpath('../utilities'));

% Data location
data_locn = '/Volumes/thefarm2/live/CrazyEights/AD_2PCa/';
if ~exist(data_locn,'dir')
    data_locn = '/Volumes/RDS/project/thefarm2/live/CrazyEights/AD_2PCa/';
end
if ~exist(data_locn,'dir')
    data_locn = '/rds/general/user/mgo/projects/thefarm2/live/CrazyEights/AD_2PCa/';
end

filedir = [data_locn 'Data/20190406/Processed/'];
refIm = '20190406_20_27_07';
% Files to be registered
files = ['20190406_20_14_42';...
         '20190406_20_33_01';...
         '20190406_20_38_41';...
         '20190406_20_41_07';...
         '20190406_20_45_35';...
         '20190406_20_51_29';...
         '20190406_20_56_30'];
     
% Read the masks
filedir2 = [data_locn 'Summaries/m82/'];
ROIfile = [filedir2 'RoiSet.zip'];
[cvsROIs] = ReadImageJROI(ROIfile);
[rawmasks] = ROIs2Regions(cvsROIs, [512 512]);
Nmasks = size(cvsROIs,2);

for n = 1:Nmasks
    rawmasks_ind{n} = rawmasks.PixelIdxList{n};
    mask1 = zeros(512,512);
    mask1(rawmasks_ind{n}) = 1;
    mask = mask1';
    masks{n} = find(mask);
end

clear ROIfile cvsROIs rawmasks
    
%% Reference image
str = sprintf('Reading reference image: %s\n', refIm);
cprintf(str)
load([filedir refIm '/' refIm '_2P_mcorr_output.mat'])
ref.green = green;
ref.red = red;
clear green red params

ref_imG = read_file( [filedir refIm '/' refIm '_2P_XYT_green_mcorr.tif'] );
ref_imR = read_file( [filedir refIm '/' refIm '_2P_XYT_red_mcorr.tif'] );

for n = 1:Nmasks
    maskind = masks{n};
    for j = 1:size(ref_imG,3)
        imG_reshaped = reshape( ref_imG(:,:,j), 512*512, 1);
        ref.tsG( n, j ) = mean( imG_reshaped(maskind) );
        imR_reshaped = reshape( ref_imR(:,:,j), 512*512, 1);
        ref.tsR( n, j ) = mean( imR_reshaped(maskind) );
    end
end
ref.R = ratiometric_Ca( ref.tsG, ref.tsR, 11 );
ref.spikes = nndORoasis(ref.R, 2, 0.94, 2.4);

clear ref_imG ref_imR imG_reshaped imR_reshaped maskind i j

%% Other images
for n = 1:size(files,1)
    str = sprintf('Processing %s ...\n', files(n,:));
    cprintf(str)

    open_imG = read_file( [filedir files(n,:) '/' files(n,:) '_2P_XYT_green_mcorr_globalreg.tif'] );
    open_imR = read_file( [filedir files(n,:) '/' files(n,:) '_2P_XYT_red_mcorr_globalreg.tif'] );

    for i = 1:Nmasks
        maskind = masks{i};
        for j = 1:size(open_imG,3)
            imG_reshaped = reshape( open_imG(:,:,j), 512*512, 1);
            open(n).tsG( i, j ) = mean( imG_reshaped(maskind) );
            imR_reshaped = reshape( open_imR(:,:,j), 512*512, 1);
            open(n).tsR( i, j ) = mean( imR_reshaped(maskind) );
        end
    end
    clear open_imG open_imR imG_reshaped imR_reshaped maskind i j

    open(n).R = ratiometric_Ca( open(n).tsG, open(n).tsR, 11 );
    open(n).spikes = nndORoasis(open(n).R, 2, 0.94, 2.4);
end

clear Nmasks data_locn filedir filedir2 files refIm str n 
%% Save output
open_2027 = ref;
open_2014 = open(1);
open_2033 = open(2);
open_2038 = open(3);
open_2041 = open(4);
open_2045 = open(5);
open_2051 = open(6);
open_2056 = open(7);

save('m82_openfield_20190406_timeseries.mat')

toc

%% Refine spike extraction
% load('m82_openfield_20190406_segment.mat')

%% Combine data for each environment
% cprintf('Combining data for open field')
% open.tsG    = [open_2014.tsG,    open_2027.globaltsG,    open_2033.globaltsG,    open_2038.globaltsG,    open_2041.globaltsG,    open_2045.globaltsG,    open_2051.globaltsG,    open_2056.globaltsG   ];
% open.tsR    = [open_2014.tsR,    open_2027.globaltsR,    open_2033.globaltsR,    open_2038.globaltsR,    open_2041.globaltsR,    open_2045.globaltsR,    open_2051.globaltsR,    open_2056.globaltsR   ];
% open.R      = [open_2014.R,      open_2027.globalR,      open_2033.globalR,      open_2038.globalR,      open_2041.globalR,      open_2045.globalR,      open_2051.globalR,      open_2056.globalR     ];
% open.spikes = [open_2014.spikes, open_2027.globalspikes, open_2033.globalspikes, open_2038.globalspikes, open_2041.globalspikes, open_2045.globalspikes, open_2051.globalspikes, open_2056.globalspikes];
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
% 
% %% Generate PF maps with ASD
% cprintf('Generating PF maps')
% spikes = open.activeData.spikes;
% xp1 = open.activeData.x;
% xp2 = open.activeData.y;
% 
% % addpath(genpath([pwd,'/fastASD'])); % add ASD folder
% 
% npart       = 1;                 % number of epochs to divide trial in
% ASD         = 1;                 % perform ASD estimation
% track_plot  = 0;                 % display tracking and firing for sample cells
% fsamp       = 30;                % in Hz
% n1=100; n2=100; nks=[n1,n2];     % env discretisation for ASD estimation
% h1 = 20; h2 = 20; hs = [h1,h2];
% %%%%%%%%%%%%
% % spikes = activeData.spikes; % extract spiking matrix
% %%%%%%%%%%%%
% % initialise  variables to store results
% kasd    = zeros(npart, n1, n2, size(spikes,1));
% hpf     = zeros(npart, h1, h2, size(spikes,1));
% hpfs    = zeros(npart, h1, h2, size(spikes,1));
% infoH   = zeros(npart, size(spikes,1), 2);
% infoASD = zeros(npart, size(spikes,1), 2);
% 
% % process tracking-environment data
% %%%%%%%%%%%%%
% % xp1 = activeData.x;
% % xp2 = activeData.y;
% %%%%%%%%%%%%%
% xp1 = (xp1-min(xp1))/(max(xp1)-min(xp1)); % normalised 0 mean
% xp2 = (xp2-min(xp2))/(max(xp2)-min(xp2));
% xp  = [xp1,xp2];
% 
% % discretize x and y position for histogram estimation
% x1 = linspace(0,1.0001,h1+1);
% x2 = linspace(0,1.0001,h2+1);
% xh = floor((xp(:,1)-x1(1))/(x1(2)-x1(1)))+1;
% yh = floor((xp(:,2)-x2(1))/(x2(2)-x2(1)))+1;
% occMap = full(sparse(xh,yh,1,h1,h2));
% mode = 0; % the mask is obtained by imfill only
% envMask_h = getEnvEdgePrior(occMap,mode); % hist
% 
% % discretize x and y position for ASD estimation
% x1 = linspace(0,1.0001,n1+1);
% x2 = linspace(0,1.0001,n2+1);
% xi = floor((xp(:,1)-x1(1))/(x1(2)-x1(1)))+1;
% yi = floor((xp(:,2)-x2(1))/(x2(2)-x2(1)))+1;
% xind = sub2ind([n1,n2],xi,yi); % flatten bin tracking (for ASD)
% occMap = full(sparse(xi,yi,1,n1,n2));
% mode = 2; % the mask is obtained by dilation and imfill
% envMask_asd = getEnvEdgePrior(occMap,mode); % ASD
% 
% % find which neurons are spiking
% for ii = 1:size(spikes,1)
%     a(ii) = sum(spikes(ii,:));
% end
% spk_idx = find(a); % store indices
% 
% for id = 1:length(spk_idx)
%     z = spikes(spk_idx(id),:);
%     if track_plot % plot sample tracking & spikes
%         figure; hold on;
%         plot(xp1,xp2); plot(xp1(z>0),xp2(z>0),'r.','markersize',10);
%     end
%     % separate exploration in smaller intervals
%     pp_idx = round(linspace(1,length(z),npart+1));
%     for pp = 1:npart
%             xpp = xind(pp_idx(pp):pp_idx(pp+1));
%             zpp = z(pp_idx(pp):pp_idx(pp+1));
%         % ASD estimation
%         if ASD
%             [asd,~] = runASD_2d(xpp',zpp',nks,envMask_asd);
%              if min(asd)<0; asd = asd-min(asd); end
%              kasd(pp,:,:,spk_idx(id)) = asd;
%         end
%         % histogram estimation
%         occMap = full(sparse(xh,yh,1,h1,h2));
%         counts = full(sparse(xh,yh,z,h1,h2));
%         hhh = counts./occMap;
%         hhh(isnan(hhh)) = 0;
%         hpf(pp,:,:,spk_idx(id)) = hhh;
%         hhh = imgaussfilt(hhh,h1/10); hhh(~envMask_h) = 0;
%         hpfs(pp,:,:,spk_idx(id)) = hhh;
%         % info estimation
%         [infoASD(pp,spk_idx(id),1), infoASD(pp,spk_idx(id),2)] =...
%             infoMeasures(asd',ones(n1,n2),0);
%         [infoH(pp,spk_idx(id),1), infoH(pp,spk_idx(id),2)] = infoMeasures(hhh,occMap,0);
%     end
% end
% 
% %% plot PF estimate results
% 
% %n_info = 1/10; % ratio of units to consider as not informative
% info_type = 2; % 1 is info/sec, 2 is info/spk\
% 
% %% ASD
% % if ASD
% %     for pp = 1:npart
% %         for ii = 0:length(spk_idx)/16
% %             figure
% %             ha = tight_subplot(4,4,[.03 .005],[.01 .05],[.01 .01]);
% %             for jj=1:16
% %                 if (ii*16+jj)<=length(spk_idx)
% %                     axes(ha(jj));
% %                     imagesc(squeeze(kasd(pp,:,:,spk_idx(ii*16+jj)))');
% %                     axis off; colorbar;
% %                     title(['cell ',num2str(spk_idx(ii*16+jj))],'fontsize',15)
% %                 end
% %             end
% %             ax = subtitle('ASD'); set(ax, 'fontsize',20);
% %         end 
% %     end
% % end
% 
% %% HIST (not smoothed)
% % for pp = 1:npart
% %     for ii = 0:length(spk_idx)/16
% %         figure
% %         ha = tight_subplot(4,4,[.01 .005],[.01 .07],[.01 .01]);
% %         for jj=1:16
% %             if (ii*16+jj)<=length(spk_idx)
% %                 axes(ha(jj));
% %                 imagesc(squeeze(hpf(pp,:,:,spk_idx(ii*16+jj)))');
% %                 axis off; colorbar;
% %                 title(['cell ',num2str(spk_idx(ii*16+jj))],'fontsize',15)
% %             end
% %         end
% %         ax = subtitle('RAW HIST'); set(ax, 'fontsize',20);
% %     end 
% % end
% 
% %% HIST (smoothed)
% % for pp = 1:npart
% %     for ii = 0:length(spk_idx)/16
% %         figure
% %         ha = tight_subplot(4,4,[.01 .005],[.01 .07],[.01 .01]);
% %         for jj=1:16
% %             if (ii*16+jj)<=length(spk_idx)
% %                 axes(ha(jj));
% %                 imagesc(squeeze(hpfs(pp,:,:,spk_idx(ii*16+jj)))');
% %                 axis off; colorbar;
% %                 title(['cell ',num2str(spk_idx(ii*16+jj))],'fontsize',15)
% %             end
% %         end
% %         ax = subtitle('SMOOTH HIST'); set(ax, 'fontsize',20);
% %     end 
% % end
% 
% %% TRACKING
% % for pp = 1:npart
% %     for ii = 0:length(spk_idx)/16
% %         figure
% %         ha = tight_subplot(4,4,[.01 .005],[.01 .07],[.01 .01]);
% %         for jj=1:16
% %             if (ii*16+jj)<=length(spk_idx)
% %                 axes(ha(jj));
% %                 z = spikes(spk_idx(ii*16+jj),:);
% %                 hold on; axis off;
% %                 plot(xp1,-xp2); plot(xp1(z>0),-xp2(z>0),'r.','markersize',10);
% %                 title(['cell ',num2str(spk_idx(ii*16+jj))],'fontsize',15)
% %             end
% %         end
% %         ax = subtitle('TRACKING'); set(ax, 'fontsize',20);
% %     end 
% % end
% 
% %%  ALL RESULTS TOGETHER
% nPlot = 4;
% for pp = 1:npart
%     for ii=0:length(spk_idx)/nPlot
%         figure
%         ha = tight_subplot(nPlot,4,[.01 .005],[.01 .07],[.01 .01]);
%         for jj=0:3
%             if (ii*nPlot+jj) <= length(spk_idx)
%                 axes(ha(jj*nPlot+1));
%                 z = spikes(spk_idx(ii*nPlot+jj+1),:);
%                 hold on; axis off;
%                 plot(xp1,-xp2); plot(xp1(z>0),-xp2(z>0),'r.','markersize',10);
%                 axes(ha(jj*nPlot+2));
%                 imagesc(squeeze(hpf(pp,:,:,spk_idx(ii*nPlot+jj+1)))');
%                 axis off; colorbar; %caxis([0 0.06]);
%                 title(['cell ',num2str(spk_idx(ii*nPlot+jj+1))],'fontsize',15)
%                 axes(ha(jj*nPlot+3)); 
%                 imagesc(squeeze(hpfs(pp,:,:,spk_idx(ii*nPlot+jj+1)))');
%                 axis off; colorbar; %caxis([0 0.005]);
%                 axes(ha(jj*nPlot+4));
%                 imagesc(squeeze(kasd(pp,:,:,spk_idx(ii*nPlot+jj+1)))');
%                 axis off; colorbar; %caxis([0 0.003]);
%             end
%         end
%     end 
% end
% 
% % % Plot selected cells
% nPlot = 3;
% % star_idx = [6,7,20,26,36,47,59,55,71,73,78,94,99,104,110,120,145,144]; % varying colorbars
% star_idx = [151, 128, 101, 86, 70, 60]; % constant colorbar
% for pp = 1:npart
%     for ii=0:length(star_idx)/nPlot
%         figure;
%         ha = tight_subplot(6,3,[.01 .01],[.01 .07],[.01 .01]);
%         for jj=0:5
%             if (ii*nPlot+jj) <= length(star_idx)
%                 axes(ha(jj*nPlot+1));
%                 z1 = spikes(star_idx(ii*nPlot+jj+1),:);
%                 hold on; axis off;
%                 plot(xp1,-xp2); plot(xp1(z1>0),-xp2(z1>0),'r.','markersize',10);
%                 axes(ha(jj*nPlot+2)); colorbar;
%                 z2 = hpf(pp,:,:,star_idx(ii*nPlot+jj+1));
%                 imagesc(squeeze(z2)');
%                 axis off; colorbar; %caxis([0 0.01]);
%                 %title(['cell ',num2str(star_idx(ii*nPlot+jj+1))],'fontsize',15)
%                 axes(ha(jj*nPlot+3));
%                 imagesc(squeeze(hpfs(pp,:,:,star_idx(ii*nPlot+jj+1)))');
%                 axis off; colorbar; % caxis([0 0.007]);
% %                 axes(ha(jj*nPlot+4));
% %                 imagesc(squeeze(kasd(pp,:,:,star_idx(ii*nPlot+jj+1)))');
% %                 axis off; colorbar; %caxis([0 0.005]);
%             end
%         end
%     end 
% end
% 
% figure; 
% id = 151;
% subplot(311); plot(open.tsG(id,:));
% subplot(312); plot(open.R(id,:));
% subplot(313); plot(open.spikes(id,:));
% 
% 
% % save 'openfield_20190406_m82_ASD.mat'
