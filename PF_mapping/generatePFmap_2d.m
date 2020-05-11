% Written by Ann Go, adapted from Giuseppe's PF_ASD_2d.m

function [hist, asd, activeData, PFdata] = generatePFmap_2d(spikes, downTrackdata, params, doasd)

if nargin < 4, doasd = true; end
hist.Nbins = params.PFmap.Nbins;
Nepochs = params.PFmap.Nepochs;
Vthr = params.PFmap.Vthr;
histsmoothWin = params.PFmap.histsmoothWin;
prctile_thr = params.PFmap.prctile_thr;
Ncells = size(spikes,1);

% Input data
x = downTrackdata.x;
y = downTrackdata.y;
phi = downTrackdata.phi;
r = downTrackdata.r;
speed = downTrackdata.speed;
t = downTrackdata.time;
ind = find(abs(diff(t))>200);
    if numel(ind)>1
        dt = mean(diff(t(1:ind(1))));
    else
        dt = mean(diff(t));
    end

% Consider only samples when the mouse is active
activex     = x(speed > Vthr);
activey     = y(speed > Vthr);
activephi   = phi(speed > Vthr);
activespk   = spikes(:,speed > Vthr);
activet     = t(speed > Vthr);
activespeed = speed(speed > Vthr);
activer     = r(speed > Vthr);
clear x y phi r speed t

xp1 = activex;
xp2 = activey;

h1 = hist.Nbins(1); h2 = hist.Nbins(2);     % env discretisation
n1 = 64; n2 = 64; asd.Nbins = [n1,n2];      

% process tracking-environment data
xp1 = (xp1-min(xp1))/(max(xp1)-min(xp1)); % normalised 0 mean
xp2 = (xp2-min(xp2))/(max(xp2)-min(xp2));
xp  = [xp1,xp2];

% discretize x and y position for histogram estimation
x1 = linspace(0,1.0001,h1+1);
x2 = linspace(0,1.0001,h2+1);
xh = floor((xp(:,1)-x1(1))/(x1(2)-x1(1)))+1;
yh = floor((xp(:,2)-x2(1))/(x2(2)-x2(1)))+1;
hist.bin_pos = sub2ind([h1,h2],xh,yh); % flatten bin tracking 
hist.occMap = full(sparse(xh,yh,1,h1,h2));
mode = 0; % the mask is obtained by imfill only
envMask_h = getEnvEdgePrior(hist.occMap,mode); % hist

if doasd
    % discretize x and y position for ASD estimation
    x1 = linspace(0,1.0001,n1+1);
    x2 = linspace(0,1.0001,n2+1);
    xa = floor((xp(:,1)-x1(1))/(x1(2)-x1(1)))+1;
    ya = floor((xp(:,2)-x2(1))/(x2(2)-x2(1)))+1;
    asd.bin_pos = sub2ind([n1,n2],xa,ya); % flatten bin tracking (for ASD)
    asd.occMap = full(sparse(xa,ya,1,n1,n2));
    mode = 2; % the mask is obtained by dilation and imfill
    envMask_asd = getEnvEdgePrior(asd.occMap,mode); % ASD
end

% find which neurons are spiking
a = zeros(size(spikes,1),1);
for ii = 1:Ncells
    a(ii) = sum(spikes(ii,:));
end
spkIdx = find(a); % store indices
Nspk = length(spkIdx);

% initialise  variables to store results
spkPeak = zeros(Nspk,1);
spkMean = zeros(Nspk,1);

hist.spkMap = zeros(h1, h2, Nspk);
hist.rMap = zeros(h1, h2, Nspk);
hist.rMap_sm = zeros(h1, h2, Nspk);
hist.infoMap = zeros(Nspk, 2);

if doasd
    asd.spkMap = zeros(n1, n2, Nspk);
    asd.rMap = zeros(n1, n2, Nspk);
    asd.infoMap = zeros(Nspk, 2);
end

for id = 1:Nspk
    z = activespk(spkIdx(id),:);
    spkPeak(id) = max(z);
    spkMean(id) = mean(z);
    
    % Histogram estimation
    hist.spkMap(:,:,id) = full(sparse(xh,yh,z,h1,h2));
    hhh = hist.spkMap(:,:,id)./(hist.occMap(:,:)*dt);
    hhh(isnan(hhh)) = 0; hhh(isinf(hhh)) = 0; 
    hist.rMap(:,:,id) = hhh;
    hhh = imgaussfilt(hhh,2); 
    hhh(~envMask_h) = 0;
    hist.rMap_sm(:,:,id) = hhh;
    hist.normrMap_sm(:,:,id) = hist.rMap_sm(:,:,id)./max(max(hist.rMap_sm(:,:,id)));

    % ASD estimation
    if doasd
        [aaa,~] = runASD_2d(asd.bin_pos, z', asd.Nbins, envMask_asd);
        if min(aaa)<0; aaa = aaa-min(aaa); end
        asd.rMap(:,:,id) = aaa;
        asd.normrMap(:,:,id) = asd.rMap(:,:,id)./max(max(asd.rMap(:,:,id)));
    end

    % info estimation
    [hist.infoMap(id,1), hist.infoMap(id,2)] = infoMeasures(hhh,hist.occMap,0);
    if doasd
        [asd.infoMap(id,1), asd.infoMap(id,2)] = infoMeasures(aaa',ones(n1,n2),0);
    end
end

% Find location preference and field size
[ hist.centroid, hist.fieldSize ] = prefLoc_fieldSize_2d( hist.rMap_sm );
if doasd
    [ asd.centroid, asd.fieldSize ] = prefLoc_fieldSize_2d( asd.rMap );
end


%% PLACE CELLS
% Identify place cells. The cells are sorted in descending order of info content
[hist.SIsec.pcIdx, hist.SIspk.pcIdx, hist.SIsec.nonpcIdx, hist.SIspk.nonpcIdx] ...
    = identifyPCs_2d( hist.bin_pos, activespk, hist.infoMap, hist.Nbins, prctile_thr, 1000, 'hist' );

if doasd
    [asd.SIsec.pcIdx, asd.SIspk.pcIdx, asd.SIsec.nonpcIdx, asd.SIspk.nonpcIdx] ...
    = identifyPCs_2d( asd.bin_pos, activespk, asd.infoMap, asd.Nbins, prctile_thr, 200, 'asd' );
end


%% Finalise place field maps, recalculate if Nepochs > 1
if Nepochs >1
    % Calculate PF maps for each epoch
    % Initialise matrices
    spkPeak = zeros(Nspk, Nepochs);
    spkMean = zeros(Nspk, Nepochs);

    bin_pos_e = zeros(h1*h2, Nepochs);
    hist.spkMap = zeros(h1, h2, Nspk, Nepochs);
    hist.occMap = zeros(h1, h2, Nepochs);
    hist.rMap = zeros(h1, h2, Nspk, Nepochs);
    hist.rMap_sm = zeros(h1, h2, Nspk, Nepochs);
    hist.infoMap = zeros(Nspk, 2, Nepochs);
    hist.centroid = zeros(Nspk, Nepochs);
    hist.fieldSize = zeros(Nspk, 2, Nepochs);

    if doasd
        asd.spkMap = zeros(n1, n2, Nspk, Nepochs);
        asd.rMap = zeros(n1, n2, Nspk, Nepochs);
        asd.infoMap = zeros(Nspk, 2, Nepochs);
        asd.centroid = zeros(Nspk, Nepochs);
        asd.fieldSize = zeros(Nspk, 2, Nepochs);
    end

    % Calculate PF maps
    e_bound = round( linspace(1,size(activespk,2),Nepochs+1) );
    for id = 1:Nspk
        z = activespk(id,:);
        
        % separate exploration in smaller intervals
        for e = 1:Nepochs
            activer(:,e) = activer(e_bound(e):e_bound(e+1));

            % histogram estimation
            bin_pos_e(:,e) = hist.bin_pos(e_bound(e):e_bound(e+1));
            hist.occMap(:,:,e) = full(sparse(xh(e_bound(e):e_bound(e+1)),yh(e_bound(e):e_bound(e+1)),1,h1,h2));
            hist.spkMap(:,:,id,e) = full(sparse(xh(e_bound(e):e_bound(e+1)), yh(e_bound(e):e_bound(e+1)), z(e_bound(e):e_bound(e+1)), h1, h2));
            hhh = hist.spkMap(:,:,id,e)./hist.occMap(:,:,e);
            hhh(isnan(hhh)) = 0; hhh(isinf(hhh)) = 0;
            hist.rMap(:,:,id,e) = hhh;
            hhh = imgaussfilt(hhh,2); 
            hhh(~envMask_h) = 0;
            hist.rMap_sm(:,:,id,e) = hhh;
            hist.normrMap_sm(:,:,id,e) = hist.rMap_sm(:,:,id,e)./max(max(hist.rMap_sm(:,:,id,e)));
            
            % ASD estimation
            if doasd
                [aaa,~] = runASD_2d(asd.bin_pos(e_bound(e):e_bound(e+1)), z(e_bound(e):e_bound(e+1))', asd.Nbins, envMask_asd);
                if min(aaa)<0; aaa = aaa-min(aaa); end
                asd.rMap(:,:,id,e) = aaa;
                asd.normrMap(:,:,id,e) = asd.rMap(:,:,id,e)./max(max(asd.rMap(:,:,id,e)));
            end

            % info estimation
            [hist.infoMap(id,1,e), hist.infoMap(id,2,e)] = infoMeasures(hhh,hist.occMap(:,:,e),0);
            if doasd
                [asd.infoMap(id,1,e), asd.infoMap(id,2,e)] = infoMeasures(aaa',ones(n1,n2),0);
            end
    
            % Find location preference and field size
            [ hist.centroid(id,e), hist.fieldSize(id,e) ] = prefLoc_fieldSize_2d( hist.rMap_sm(:,:,id,e) );
            if doasd
                [ asd.centroid(id,e), asd.fieldSize(id,e) ] = prefLoc_fieldSize_2d( asd.rMap(:,:,id,e) );
            end
        end
        
        hist.bin_pos = bin_pos_e;
    end
end


%% Final data
% histogram estimation
if ~isempty(hist.SIsec.pcIdx)
    hist.SIsec.spkMean = spkMean(hist.SIsec.pcIdx);
    hist.SIsec.spkPeak = spkPeak(hist.SIsec.pcIdx);
    
    hist.SIsec.spkMap = hist.spkMap(:,:,hist.SIsec.pcIdx);
    hist.SIsec.pfMap = hist.rMap(:,:,hist.SIsec.pcIdx);
    hist.SIsec.pfMap_sm = hist.rMap_sm(:,:,hist.SIsec.pcIdx);
    hist.SIsec.normpfMap_sm = hist.normrMap_sm(:,:,hist.SIsec.pcIdx);
    hist.SIsec.infoMap = hist.infoMap(hist.SIsec.pcIdx,1,:);
    hist.SIsec.centroid = hist.centroid(hist.SIsec.pcIdx,:);
    hist.SIsec.fieldSize = hist.fieldSize(hist.SIsec.pcIdx,:,:);
end

if ~isempty(hist.SIspk.pcIdx)
    hist.SIspk.spkMean = spkMean(hist.SIspk.pcIdx);
    hist.SIspk.spkPeak = spkPeak(hist.SIspk.pcIdx);
    
    hist.SIspk.spkMap = hist.spkMap(:,:,hist.SIspk.pcIdx);
    hist.SIspk.pfMap = hist.rMap(:,:,hist.SIspk.pcIdx);
    hist.SIspk.pfMap_sm = hist.rMap_sm(:,:,hist.SIspk.pcIdx);
    hist.SIspk.normpfMap_sm = hist.normrMap_sm(:,:,hist.SIspk.pcIdx);
    hist.SIspk.infoMap = hist.infoMap(hist.SIspk.pcIdx,2,:);
    hist.SIspk.centroid = hist.centroid(hist.SIspk.pcIdx,:);
    hist.SIspk.fieldSize = hist.fieldSize(hist.SIspk.pcIdx,:,:);
end

if doasd
    %asd
    if ~isempty(asd.SIsec.pcIdx)
        asd.SIsec.spkMean = spkMean(asd.SIsec.pcIdx);
        asd.SIsec.spkPeak = spkPeak(asd.SIsec.pcIdx);

        asd.SIsec.pfMap = asd.rMap(:,:,asd.SIsec.pcIdx);
        asd.SIsec.normpfMap = asd.normrMap(:,:,asd.SIsec.pcIdx);
        asd.SIsec.infoMap = asd.infoMap(asd.SIsec.pcIdx,1,:);
        asd.SIsec.centroid = asd.centroid(asd.SIsec.pcIdx,:);
        asd.SIsec.fieldSize = asd.fieldSize(asd.SIsec.pcIdx,:,:);
    end

    if ~isempty(asd.SIspk.pcIdx)
        asd.SIspk.spkMean = spkMean(asd.SIspk.pcIdx);
        asd.SIspk.spkPeak = spkPeak(asd.SIspk.pcIdx);

        asd.SIspk.pfMap = asd.rMap(:,:,asd.SIspk.pcIdx);
        asd.SIspk.normpfMap = asd.normrMap(:,:,asd.SIspk.pcIdx);
        asd.SIspk.infoMap = asd.infoMap(asd.SIspk.pcIdx,2,:);
        asd.SIspk.centroid = asd.centroid(asd.SIspk.pcIdx,:);
        asd.SIspk.fieldSize = asd.fieldSize(asd.SIspk.pcIdx,:,:);
    end
end

%% Outputs
activeData.x = activex;
activeData.y = activey;
activeData.r = activer;
activeData.phi = activephi;
activeData.speed = activespeed;
activeData.t = activet;
activeData.spikes = activespk;

PFdata.spkMean = spkMean;
PFdata.spkPeak = spkPeak;

if ~doasd, asd = []; end

end