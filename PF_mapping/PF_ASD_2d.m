% close all

% load data ...
%load('20190220_14_12_22_openfield.mat')

spikes = open.activeData.spikes;
xp1 = open.activeData.x;
xp2 = open.activeData.y;

addpath(genpath([pwd,'/fastASD'])); % add ASD folder

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
        hhh = imgaussfilt(hhh,h1/5); hhh(~envMask_h) = 0;
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
                axis off; colorbar; caxis([0 0.06]);
                title(['cell ',num2str(spk_idx(ii*nPlot+jj+1))],'fontsize',15)
                axes(ha(jj*nPlot+3)); 
                imagesc(squeeze(hpfs(pp,:,:,spk_idx(ii*nPlot+jj+1)))');
                axis off; colorbar; caxis([0 0.005]);
                axes(ha(jj*nPlot+4));
                imagesc(squeeze(kasd(pp,:,:,spk_idx(ii*nPlot+jj+1)))');
                axis off; colorbar; caxis([0 0.003]);
            end
        end
    end 
end

% Plot selected cells
nPlot = 4;
star_idx = [187,195,168,155,123,91,81,64,56,15]; % constant colorbar
for pp = 1:npart
    for ii=0:length(star_idx)/nPlot
        figure;
        ha = tight_subplot(nPlot,4,[.01 .005],[.01 .07],[.01 .01]);
        for jj=0:3
            if (ii*nPlot+jj) <= length(star_idx)
                axes(ha(jj*nPlot+1));
                z = spikes(star_idx(ii*nPlot+jj+1),:);
                hold on; axis off;
                plot(xp1,-xp2); plot(xp1(z>0),-xp2(z>0),'r.','markersize',10);
                axes(ha(jj*nPlot+2));
                imagesc(squeeze(hpf(pp,:,:,star_idx(ii*nPlot+jj+1)))');
                axis off; colorbar; caxis([0 0.08]);
                title(['cell ',num2str(star_idx(ii*nPlot+jj+1))],'fontsize',15)
                axes(ha(jj*nPlot+3));
                imagesc(squeeze(hpfs(pp,:,:,star_idx(ii*nPlot+jj+1)))');
                axis off; colorbar; caxis([0 0.007]);
                axes(ha(jj*nPlot+4));
                imagesc(squeeze(kasd(pp,:,:,star_idx(ii*nPlot+jj+1)))');
                axis off; colorbar; caxis([0 0.005]);
            end
        end
    end 
end

% figure; 
% id = 195;
% subplot(311); plot(open.tsG(id,:));
% subplot(312); plot(open.R(id,:));
% subplot(313); plot(open.spikes(id,:));


% save 'openfield_20190404_ASD.mat'