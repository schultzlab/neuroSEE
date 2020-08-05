function [PFdata, hist, asd, xh, yh, xa, ya] = calcPFdata_2d(activex, activey, activespk, Nbins_hist, Nbins_asd, gaussfiltSigma, fr, doasd)
if nargin<8, doasd = false; end
if nargin<7, params = neuroSEE_setparams; fr = params.PFmap.fr; end
if nargin<6, params = neuroSEE_setparams; gaussfiltSigma = params.PFmap.gaussfiltSigma; end

% process tracking environment data
xp = activex;
yp = activey;
xp = (xp-min(xp))/(max(xp)-min(xp));   
yp = (yp-min(yp))/(max(yp)-min(yp));
h1 = Nbins_hist(1); h2 = Nbins_hist(2);     % env discretisation
hist.Nbins = [h1,h2];

% discretize x and y position for histogram estimation
xx = linspace(0,1.0001,h1+1);
yy = linspace(0,1.0001,h2+1);
xh = floor((xp-xx(1))/(xx(2)-xx(1)))+1;
yh = floor((yp-yy(1))/(yy(2)-yy(1)))+1;
hist.bin_pos = sub2ind([h1,h2],xh,yh);      % flatten bin tracking 
hist.occMap = full(sparse(yh,xh,1,h1,h2));  % the reverse order of y&x lets the image pixel indexing 
                                            % match the matrix row and column indexing
mode = 0; % the mask is obtained by imfill only
envMask_h = getEnvEdgePrior(hist.occMap,mode); 

if doasd
    n1 = Nbins_asd(1); n2 = Nbins_asd(2); asd.Nbins = [n1,n2];      
    % discretize x and y position for ASD estimation
    xx = linspace(0,1.0001,n1+1);
    yy = linspace(0,1.0001,n2+1);
    xa = floor((xp-xx(1))/(xx(2)-xx(1)))+1;
    ya = floor((yp-yy(1))/(yy(2)-yy(1)))+1;
    asd.bin_pos = sub2ind([n1,n2],xa,ya);   % flatten bin tracking (for ASD)
    asd.occMap = full(sparse(ya,xa,1,n1,n2));
    mode = 2; % the mask is obtained by dilation and imfill
    envMask_asd = getEnvEdgePrior(asd.occMap,mode); % ASD
end

% initialise  variables to store results
Ncells = size(activespk,1);
spk_eventrate = zeros(Ncells,1);            % average event ("spike" rate)
spk_rate = zeros(Ncells,1);                 % average event amplitude rate
bin_activet = zeros(Ncells, h1*h2);         % fraction of dwell time in bin for which cell was active   
hist.spkMap = zeros(h1, h2, Ncells);
hist.rMap = zeros(h1, h2, Ncells);
hist.rMap_sm = zeros(h1, h2, Ncells);
hist.infoMap = zeros(Ncells, 2);
L = size(activespk,2);
T = L / fr;

if doasd
    asd.spkMap = zeros(n1, n2, Ncells);
    asd.rMap = zeros(n1, n2, Ncells);
    asd.infoMap = zeros(Ncells, 2);
end

%% rate map calculations
dt = 1/fr;
hist.occMap = hist.occMap*dt;
hist.occMap(hist.occMap < 1) = 0;   % bin locations hardly visited might lead to artificially high event rates
for c = 1:Ncells
    z = activespk(c,:);
    spk_idx = z>0;
    spk_ampl = z(spk_idx);
    spk_eventrate(c) = length(spk_ampl)/T;
    spk_rate(c) = sum(spk_ampl)/T;
    for bin = 1:h1*h2
       bin_idx = find(hist.bin_pos == bin);
       bin_activet(c,bin) = length(find(z(bin_idx)>0))/length(bin_idx);
    end
    
    % Histogram estimation
    hist.spkMap(:,:,c) = full(sparse(yh,xh,z,h1,h2));   % the reverse order of y&x lets the image pixel indexing 
                                                        % match the matrix row and column indexing
    hhh = hist.spkMap(:,:,c)./hist.occMap;
    hhh(isnan(hhh)) = 0; hhh(isinf(hhh)) = 0; 
    hhh(~envMask_h) = -0.01;
    hist.rMap(:,:,c) = hhh;
    hhh(~envMask_h) = 0;
    hhh = imgaussfilt(hhh,gaussfiltSigma); 
    hhh(~envMask_h) = 0;
    [hist.infoMap(c,1), hist.infoMap(c,2)] = infoMeasures(hhh,hist.occMap,0);
    hhh(~envMask_h) = -0.01;     % this is so pf map can have a circular boundary
    hist.rMap_sm(:,:,c) = hhh;
    hist.normrMap_sm(:,:,c) = hist.rMap_sm(:,:,c)./max(max(hist.rMap_sm(:,:,c)));
    
    

    % ASD estimation
    if doasd
        [aaa,~] = runASD_2d(asd.bin_pos, z', asd.Nbins, envMask_asd);
        if min(aaa)<0; aaa = aaa-min(aaa); end
        asd.rMap(:,:,c) = aaa;
        asd.normrMap(:,:,c) = asd.rMap(:,:,c)./max(max(asd.rMap(:,:,c)));
        [asd.infoMap(c,1), asd.infoMap(c,2)] = infoMeasures(aaa',ones(n1,n2),0);
    end

    % info estimation
    
    if doasd
        
    end
end

%% Find location preference and field size
[ hist.centroid, hist.fieldSize, hist.pfBins ] = prefLoc_fieldSize_2d( hist.rMap_sm );
if doasd
    [ asd.centroid, asd.fieldSize, asd.pfBins ] = prefLoc_fieldSize_2d( asd.rMap );
end

%% calculate fraction of dwell time within place field in which cell was active
for c = 1:Ncells
   % get list of spikes for this cell
   z = activespk(c,:);
   pfBins_activet = 0;
   pfBins_t = 0;
   for k = 1:length(hist.pfBins{c})
       bin = hist.pfBins{c}(k);
       bin_idx = find(hist.bin_pos == bin);
       pfBins_activet = pfBins_activet + length(find(z(bin_idx)>0));
       pfBins_t = pfBins_t + length(bin_idx);
   end
   hist.pf_activet(c) = pfBins_activet/pfBins_t;
   
   if doasd
       for k = 1:length(asd.pfBins{c})
           bin = asd.pfBins{c}(k);
           bin_idx = find(asd.bin_pos == bin);
           pfBins_activet = pfBins_activet + length(find(z(bin_idx)>0));
           pfBins_t = pfBins_t + length(bin_idx);
       end
       asd.pf_activet(c) = pfBins_activet/pfBins_t;
   end
end

% output
PFdata.spk_eventrate = spk_eventrate;
PFdata.spk_rate = spk_rate;
PFdata.bin_activet = bin_activet;
if ~doasd, asd = []; xa = []; ya = []; end
