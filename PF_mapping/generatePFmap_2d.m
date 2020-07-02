% Written by Ann Go, adapted from Giuseppe's PF_ASD_2d.m

function [hist, asd, PFdata] = generatePFmap_2d(spikes, downTrackdata, params)
if nargin<3, params = neuroSEE_setparams; end

doasd = params.methods.doasd; 
fr = params.PFmap.fr;
Nepochs = params.PFmap.Nepochs;
Vthr = params.PFmap.Vthr;
gaussfiltSigma = params.PFmap.gaussfiltSigma;
prctile_thr = params.PFmap.prctile_thr;
pfactivet_thr = params.PFmap.pfactivet_thr;
Nrand = params.PFmap.Nrand;
hist.Nbins = params.PFmap.Nbins;
asd.Nbins = [64,64];
Ncells = size(spikes,1);

% Input data
x = downTrackdata.x;
y = downTrackdata.y;
phi = downTrackdata.phi;
r = downTrackdata.r;
speed = downTrackdata.speed;
t = downTrackdata.time;
dt = 1/fr;

% Consider only samples when the mouse is active
activex     = x(speed > Vthr);
activey     = y(speed > Vthr);
activephi   = phi(speed > Vthr);
activespk   = spikes(:,speed > Vthr);
activet     = t(speed > Vthr);
activespeed = speed(speed > Vthr);
activer     = r(speed > Vthr);
clear x y phi r speed t

%% ALL CELLS: calculate PF data for entire duration (one epoch)
[PFdata, hist, asd, xh, yh, xa, ya] = calcPFdata_2d(activex, activey, activespk, hist.Nbins, asd.Nbins, gaussfiltSigma, fr, doasd);

%% PLACE CELLS
% Identify place cells. The cells are sorted in descending order of info content
[hist.SIsec.pcIdx, hist.SIspk.pcIdx, hist.SIsec.nonpcIdx, hist.SIspk.nonpcIdx] ...
    = identifyPCs_2d( activespk, xh, yh, hist.infoMap, hist.pf_activet, hist.Nbins, prctile_thr, pfactivet_thr, Nrand, 'hist', 2, gaussfiltSigma );

if doasd
    [asd.SIsec.pcIdx, asd.SIspk.pcIdx, asd.SIsec.nonpcIdx, asd.SIspk.nonpcIdx] ...
    = identifyPCs_2d( activespk, xa, ya, asd.infoMap, asd.pf_activet, asd.Nbins, prctile_thr, pfactivet_thr, Nrand, 'asd', gaussfiltSigma );
end


%% Finalise place field maps, recalculate if Nepochs > 1
if Nepochs >1
    % Calculate PF maps for each epoch
    % Initialise matrices

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
    for c = 1:Nspk
        z = activespk(c,:);
        
        % separate exploration in smaller intervals
        for e = 1:Nepochs
            activer(:,e) = activer(e_bound(e):e_bound(e+1));

            % histogram estimation
            bin_pos_e(:,e) = hist.bin_pos(e_bound(e):e_bound(e+1));
            hist.occMap(:,:,e) = full(sparse(xh(e_bound(e):e_bound(e+1)),yh(e_bound(e):e_bound(e+1)),1,h1,h2));
            hist.spkMap(:,:,c,e) = full(sparse(xh(e_bound(e):e_bound(e+1)), yh(e_bound(e):e_bound(e+1)), z(e_bound(e):e_bound(e+1)), h1, h2));
            hhh = hist.spkMap(:,:,c,e)./hist.occMap(:,:,e);
            hhh(isnan(hhh)) = 0; hhh(isinf(hhh)) = 0;
            hist.rMap(:,:,c,e) = hhh;
            hhh = imgaussfilt(hhh,gaussfiltSigma); 
            hhh(~envMask_h) = 0;
            hist.rMap_sm(:,:,c,e) = hhh;
            hist.normrMap_sm(:,:,c,e) = hist.rMap_sm(:,:,c,e)./max(max(hist.rMap_sm(:,:,c,e)));
            
            % ASD estimation
            if doasd
                [aaa,~] = runASD_2d(asd.bin_pos(e_bound(e):e_bound(e+1)), z(e_bound(e):e_bound(e+1))', asd.Nbins, envMask_asd);
                if min(aaa)<0; aaa = aaa-min(aaa); end
                asd.rMap(:,:,c,e) = aaa;
                asd.normrMap(:,:,c,e) = asd.rMap(:,:,c,e)./max(max(asd.rMap(:,:,c,e)));
            end

            % info estimation
            [hist.infoMap(c,1,e), hist.infoMap(c,2,e)] = infoMeasures(hhh,hist.occMap(:,:,e),0);
            if doasd
                [asd.infoMap(c,1,e), asd.infoMap(c,2,e)] = infoMeasures(aaa',ones(n1,n2),0);
            end
    
            % Find location preference and field size
            [ hist.centroid(c,e), hist.fieldSize(c,e) ] = prefLoc_fieldSize_2d( hist.rMap_sm(:,:,c,e) );
            if doasd
                [ asd.centroid(c,e), asd.fieldSize(c,e) ] = prefLoc_fieldSize_2d( asd.rMap(:,:,c,e) );
            end
        end
        
        hist.bin_pos = bin_pos_e;
    end
end


%% Final data
% histogram estimation
if ~isempty(hist.SIsec.pcIdx)
    hist.SIsec.pfMap = hist.rMap(:,:,hist.SIsec.pcIdx);
    hist.SIsec.pfMap_sm = hist.rMap_sm(:,:,hist.SIsec.pcIdx);
    hist.SIsec.normpfMap_sm = hist.normrMap_sm(:,:,hist.SIsec.pcIdx);
end

if ~isempty(hist.SIspk.pcIdx)
    hist.SIspk.pfMap = hist.rMap(:,:,hist.SIspk.pcIdx);
    hist.SIspk.pfMap_sm = hist.rMap_sm(:,:,hist.SIspk.pcIdx);
    hist.SIspk.normpfMap_sm = hist.normrMap_sm(:,:,hist.SIspk.pcIdx);
end

if doasd
    %asd
    if ~isempty(asd.SIsec.pcIdx)
        asd.SIsec.pfMap = asd.rMap(:,:,asd.SIsec.pcIdx);
        asd.SIsec.normpfMap = asd.normrMap(:,:,asd.SIsec.pcIdx);
    end

    if ~isempty(asd.SIspk.pcIdx)
        asd.SIspk.pfMap = asd.rMap(:,:,asd.SIspk.pcIdx);
        asd.SIspk.normpfMap = asd.normrMap(:,:,asd.SIspk.pcIdx);
    end
else
    asd = [];
end

end