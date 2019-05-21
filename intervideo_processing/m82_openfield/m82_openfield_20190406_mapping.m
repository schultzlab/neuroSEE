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

% load('m82_openfield_20190406_timeseries.mat')
% 
% % Data location
% data_locn = '/Volumes/thefarm2/live/CrazyEights/AD_2PCa/';
% if ~exist(data_locn,'dir')
%     data_locn = '/Volumes/RDS/project/thefarm2/live/CrazyEights/AD_2PCa/';
% end
% if ~exist(data_locn,'dir')
%     data_locn = '/rds/general/user/mgo/projects/thefarm2/live/CrazyEights/AD_2PCa/';
% end
% 
% filedir = [data_locn 'Data/20190406/Processed/'];
% 
%     
% %% Load tsG, tsR, R, spikes
% matfiles = ['m82_open_2014_tsRspikes_refined.mat';...
%             'm82_open_2027_tsRspikes_refined.mat';...
%             'm82_open_2033_tsRspikes_refined.mat';...
%             'm82_open_2038_tsRspikes_refined.mat';...
%             'm82_open_2041_tsRspikes_refined.mat';...
%             'm82_open_2045_tsRspikes_refined.mat';...
%             'm82_open_2051_tsRspikes_refined.mat';...
%             'm82_open_2056_tsRspikes_refined.mat'];
% 
% for n = 1:size(matfiles,1)
%     openrefined(n) = load(matfiles(n,:));
% end
% 
% clear n temp 
% 
% %% Combine tsG, tsR, R, spikes data for different videos
% cprintf('Combining spike data for open field...\n')
% 
% openspliced.tsG = [];
% openspliced.tsR = [];
% openspliced.R = [];
% openspliced.spikes = [];
% 
% for n = 1:size(matfiles,1)
%     openspliced.tsG     = [openspliced.tsG, openrefined(n).cell_tsG];
%     openspliced.tsR     = [openspliced.tsR, openrefined(n).cell_tsR];
%     openspliced.R       = [openspliced.R,   openrefined(n).R];
%     openspliced.spikes  = [openspliced.spikes, openrefined(n).spikes];
% end
% 
% %% Load tracking data and combine data for different videos
% Cafiles = ['20190406_20_14_42';...
%             '20190406_20_27_07';...
%             '20190406_20_33_01';...
%             '20190406_20_38_41';...
%             '20190406_20_41_07';...
%             '20190406_20_45_35';...
%             '20190406_20_51_29';...
%             '20190406_20_56_30'];
% trackfiles = ['20190406_20_14_46';...
%          '20190406_20_27_19';...
%          '20190406_20_33_13';...
%          '20190406_20_38_44';...
%          '20190406_20_41_10';...
%          '20190406_20_45_38';...
%          '20190406_20_51_32';...
%          '20190406_20_56_33'];
% 
% for n = 1:size(trackfiles,1)
%     opentrack(n) = load([filedir Cafiles(n,:) '/Track-' trackfiles(n,1:4) '-' trackfiles(n,5:6) '-' trackfiles(n,7:8) '-' trackfiles(n,10:11) '-' trackfiles(n,13:14) '-' trackfiles(n,16:17) '.mat']); 
% end
% clear matfiles Cafiles data_locn filedir n
% 
% % Combine track data
% cprintf('Combining track data for open field...\n')
% openspliced.time = [];
% openspliced.x = [];
% openspliced.y = [];
% openspliced.r = [];
% openspliced.phi = [];
% openspliced.speed = [];
% tmax = opentrack(1).time(end);
% 
% for n = 1:size(trackfiles,1)
%     openspliced.x      = [openspliced.x,      opentrack(n).x];
%     openspliced.y      = [openspliced.y,      opentrack(n).y];
%     openspliced.r      = [openspliced.r,      opentrack(n).r];
%     openspliced.phi    = [openspliced.phi,    opentrack(n).phi];
%     openspliced.speed  = [openspliced.speed,  opentrack(n).speed];
%     openspliced.time   = [openspliced.time,   (n-1)*tmax+opentrack(n).time];
% end
% clear trackfiles tmax
% 
% % downsample data
% tracktime   = openspliced.time;
% x           = openspliced.x;
% y           = openspliced.y;
% r           = openspliced.r;
% phi         = openspliced.phi;
% speed       = openspliced.speed;
% spikes      = openspliced.spikes;
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
% openspliced.downData.x = downx;
% openspliced.downData.y = downy;
% openspliced.downData.phi = downphi;
% openspliced.downData.speed = downspeed;
% openspliced.downData.t = t;
% openspliced.downData.r = downr;
% 
% openspliced.activeData.x = activex;
% openspliced.activeData.y = activey;
% openspliced.activeData.phi = activephi;
% openspliced.activeData.speed = activespeed;
% openspliced.activeData.t = activet;
% openspliced.activeData.spikes = activespikes;
% openspliced.activeData.r = activer;
% 
% clear activephi activer activespeed activespikes activet activex activey ans
% clear downphi downr downspeed downx downy dt Nt phi r speed spikes t t0 tracktime Vthr x y
% save 'm82_openfield_20190406_splicedSpikesTrack.mat'

load('m82_openfield_20190406_splicedSpikesTrack.mat'
%% Generate PF maps with ASD
cprintf('Generating PF maps...\n')
spikes = openspliced.activeData.spikes;
xp1 = openspliced.activeData.x;
xp2 = openspliced.activeData.y;

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
for pp = 2:npart
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

save 'm82_openfield_20190406_ASD.mat'

toc
