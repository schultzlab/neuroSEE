% Written by Ann Go, adapted from Giuseppe's PF_ASD_2d.m

function [occMap, spkMap, spkIdx, hist, asd, downData, activeData] = generatePFmap_2d(spikes, imtime, trackData, params)

fr = params.fr;
Nbins = params.PFmap.Nbins;
Nepochs = params.PFmap.Nepochs;     % number of epochs to divide trial in
Vthr = params.PFmap.Vthr;
histsmoothFac = params.PFmap.histsmoothFac;

%% Pre-process tracking data
tracktime = trackData.time;
xcont = trackData.x;
ycont = trackData.y;
r = trackData.r;
phi = trackData.phi;
speed = trackData.speed;

t0 = tracktime(1);                  % initial time in tracking data
nspikes = spikes; %bsxfun( @rdivide, bsxfun(@minus, spikes, min(spikes,[],2)), range(spikes,2) ); % normalisation
Nt = size(spikes,2);                % number of timestamps for spikes

% Convert -180:180 to 0:360
if min(phi)<0
   phi(phi<0) = phi(phi<0)+360;
end

% If no timestamps were recorded for Ca images, generate timestamps
% using known image frame rate
if isempty(imtime)
   dt = 1/fr;
   t = (t0:dt:Nt*dt)';
end

% Downsample tracking to Ca trace
downphi = interp1(tracktime,phi,t,'linear');
downx   = interp1(tracktime,xcont,t,'linear');
downy   = interp1(tracktime,ycont,t,'linear');
downspeed = interp1(tracktime,speed,t,'linear'); % mm/s
downr   = interp1(tracktime,r,t,'linear'); % mm/s

% Consider only samples when the mouse is active
activex    = downx(downspeed > Vthr);
activey    = downy(downspeed > Vthr);
activephi  = downphi(downspeed > Vthr);
activespk  = nspikes(:,downspeed > Vthr);
activet     = t(downspeed > Vthr);
activespeed = speed(downspeed > Vthr);
activer    = r(downspeed > Vthr);

xp1 = activex;
xp2 = activey;

n1=100; n2=100; nks=[n1,n2];     % env discretisation for ASD estimation
h1 = Nbins(1); h2 = Nbins(2); % hs = [h1,h2];

% process tracking-environment data
xp1 = (xp1-min(xp1))/(max(xp1)-min(xp1)); % normalised 0 mean
xp2 = (xp2-min(xp2))/(max(xp2)-min(xp2));
xp  = [xp1,xp2];

% discretize x and y position for histogram estimation
x1 = linspace(0,1.0001,h1+1);
x2 = linspace(0,1.0001,h2+1);
xh = floor((xp(:,1)-x1(1))/(x1(2)-x1(1)))+1;
yh = floor((xp(:,2)-x2(1))/(x2(2)-x2(1)))+1;
occMap_hist = full(sparse(xh,yh,1,h1,h2));
mode = 0; % the mask is obtained by imfill only
envMask_h = getEnvEdgePrior(occMap_hist,mode); % hist

% discretize x and y position for ASD estimation
x1 = linspace(0,1.0001,n1+1);
x2 = linspace(0,1.0001,n2+1);
xi = floor((xp(:,1)-x1(1))/(x1(2)-x1(1)))+1;
yi = floor((xp(:,2)-x2(1))/(x2(2)-x2(1)))+1;
xind = sub2ind([n1,n2],xi,yi); % flatten bin tracking (for ASD)
occMap_asd = full(sparse(xi,yi,1,n1,n2));
mode = 2; % the mask is obtained by dilation and imfill
envMask_asd = getEnvEdgePrior(occMap_asd,mode); % ASD

% find which neurons are spiking
for ii = 1:size(spikes,1)
    a(ii) = sum(spikes(ii,:));
end
spkIdx = find(a); % store indices
Nspk = length(spkIdx);

% initialise  variables to store results
occMap          = zeros(h1, h2, Nepochs);
spkMap          = zeros(h1, h2, Nspk, Nepochs);
hist.pfMap      = zeros(h1, h2, Nspk, Nepochs);
hist.pfMap_sm   = zeros(h1, h2, Nspk, Nepochs);
hist.infoMap    = zeros(Nspk, 2, Nepochs);
asd.pfMap       = zeros(n1, n2, Nspk, Nepochs);
asd.infoMap     = zeros(Nspk, 2, Nepochs);

for id = 1:Nspk
    z = activespk(spkIdx(id),:);
    % separate exploration in smaller intervals
    e_bound = round(linspace(1,length(z),Nepochs+1));
    for e = 1:Nepochs
        % ASD estimation
        x_e = xind(e_bound(e):e_bound(e+1));
        z_e = z(e_bound(e):e_bound(e+1));
        [aaa,~] = runASD_2d(x_e',z_e',nks,envMask_asd);
        if min(aaa)<0; aaa = aaa-min(aaa); end
        asd.pfMap(:,:,id,e) = aaa;
        
        % histogram estimation
        x_e = xh(e_bound(e):e_bound(e+1));
        y_e = yh(e_bound(e):e_bound(e+1));
        occMap(:,:,e) = full(sparse(x_e,y_e,1,h1,h2));
        spkMap(:,:,id,e) = full(sparse(x_e,y_e,z_e,h1,h2));
        hhh = spkMap(:,:,id,e)./occMap(:,:,e);
        hhh(isnan(hhh)) = 0;
        hist.pfMap(:,:,id,e) = hhh;
        hhh = imgaussfilt(hhh,h1/histsmoothFac); hhh(~envMask_h) = 0;
        hist.pfMap_sm(:,:,id,e) = hhh;
        
        % info estimation
        [asd.infoMap(id,1,e), asd.infoMap(id,2,e)] =...
            infoMeasures(aaa',ones(n1,n2),0);
        [hist.infoMap(id,1,e), hist.infoMap(id,2,e)] = infoMeasures(hhh,occMap(:,e),0);
    end
end

% Outputs
downData.x = downx;
downData.y = downy;
downData.r = downr;
downData.phi = downphi;
downData.speed = downspeed;
downData.t = t;

activeData.x = activex;
activeData.y = activey;
activeData.r = activer;
activeData.phi = activephi;
activeData.speed = activespeed;
activeData.t = activet;
activeData.spikes = activespk;


end