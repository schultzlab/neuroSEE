% Written by Ann Go (some parts adapted from Giuseppe's PF_ASD_1d.m)
%
% This function maps place fields
%
% INPUTS:
%   spikes      : spike estimates obtained with oasisAR2
%   imtime      : imaging timestamps
%   trackData   : cell of tracking data with fields x, y, r, phi, w,
%                  speed, time, alpha, TTLout
%   params.
%    fr                  : imaging frame rate [default: 30.9 Hz]
%    PFmap.Nbins         : number of location bins
%    PFmap.Nepochs       : number of epochs for each 4 min video [default: 1]
%    PFmap.Vthr          : speed threshold (mm/s) [default: 20]
%    PFmap.histsmoothFac : Gaussian smoothing window for histogram
%                           estimation [default: 10]

% OUTPUTS:
%   occMap      : occupancy map
%   spikeMap   : spike map (Ncells rows x Nbins columns)
%   infoMap     : information map 
%   place field maps
%     placeMap    : place field map obtained with histogram estimation 
%     placeMap_smooth : smoothed version of placeMap 
%     placeMap_asd    : place field map obtained with ASD
%     placeMap_pertrial
%     normplaceMap_pertrial
%   downData    : tracking data downsampled to imaging frequency, fields are
%                 x, y, r, phi, w, speed, time, alpha, TTLout
%   activaData  : downsampled tracking data for when animal was moving, fields are
%                 x, y, r, phi, w, speed, time, alpha, TTLout
% 
%   varargout: [ occMap, spikeMap, infoMap, infoMap_asd, placeMap, placeMap_smooth, placeMap_asd, ...
%                placeMap_pertrial, normplaceMap_pertrial, PCidx, downData, activeData ]
% 

function varargout = generatePFmap_1D( spikes, imtime, trackData, params )
    fr = params.fr;
    Nbins = params.PFmap.Nbins;
    Nepochs = params.PFmap.Nepochs;
    Vthr = params.PFmap.Vthr;
    histsmoothFac = params.PFmap.histsmoothFac;
    
    % identify place cells
    
    
    % tracking data
    tracktime = trackData.time;
    x = trackData.x;
    y = trackData.y;
    r = trackData.r;
    phi = trackData.phi;
    speed = trackData.speed;
    
    t0 = tracktime(1);                  % initial time in tracking data
    nspikes = spikes; %bsxfun( @rdivide, bsxfun(@minus, spikes, min(spikes,[],2)), range(spikes,2) ); % normalisation
    Ncells = size(spikes,1);            % number of cells
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
    downx = interp1(tracktime,x,t,'linear');
    downy = interp1(tracktime,y,t,'linear');
    downspeed = interp1(tracktime,speed,t,'linear'); % mm/s
    downr = interp1(tracktime,r,t,'linear'); % mm/s

    % Consider only samples when the mouse is active
    activex    = downx(downspeed > Vthr);
    activey    = downy(downspeed > Vthr);
    activephi  = downphi(downspeed > Vthr);
    activespikes = nspikes(:,downspeed > Vthr);
    activet     = t(downspeed > Vthr);
    activespeed = speed(downspeed > Vthr);
    activer = r(downspeed > Vthr);
    
    % Bin phi data
    [bin_phi,~] = discretize(activephi,Nbins);

    % Initialise matrices
    occMap = zeros(Nepochs, Nbins);             % occupancy map
    spikeMap = zeros(Ncells, Nbins, Nepochs);   % spike map
    placeMap = zeros(Nepochs, Nbins, Ncells);   % place field map
    infoMap = zeros(Nepochs, Ncells, 2);        % Skagg's info

    for id = 1:Ncells
        z = activespikes(id,:);

        % separate exploration in smaller intervals
        e_bound = round( linspace(1,length(z),Nepochs+1) );
        for e = 1:Nepochs
            bin_phi_e = bin_phi(e_bound(e):e_bound(e+1));
            spike_e = z(e_bound(e):e_bound(e+1));

            % Occupancy and spike rate maps
            occMap(e,:) = histcounts(bin_phi_e,Nbins);
            for i = 1:Nbins
                spikeMap(id,i,e) = sum(spike_e(bin_phi_e == i));
            end

            % histogram estimation
            placeMap(e,:,id) = spikeMap(id,:,e)./occMap(e,:);
            placeMap(isnan(placeMap)) = 0;
            [infoMap(e,id,1), infoMap(e,id,2)] = infoMeasures(placeMap(e,:,id),occMap(e,:),0);
            
            % ASD estimation
            [placeMap_asd(e,:,id),~] = runASD_1d(bin_phi_e,(spike_e)',Nbins);
            [infoMap_asd(e,id,1), infoMap_asd(e,id,2)] =...
              infoMeasures(squeeze(placeMap_asd(e,:,id))',ones(Nbins,1),0);
        end
    end
    
    placeMap = permute(placeMap,[3 2 1]);
    infoMap = permute(infoMap,[2 3 1]);
    placeMap_smooth = smoothdata(placeMap,2,'gaussian',Nbins/histsmoothFac);
    
    placeMap_asd = permute(placeMap_asd,[3 2 1]);
    infoMap_asd = permute(infoMap_asd,[2 3 1]);
   
    % Outputs
    downData.x = downx;
    downData.y = downy;
    downData.phi = downphi;
    downData.speed = downspeed;
    downData.t = t;
    downData.r = downr;
    
    activeData.x = activex;
    activeData.y = activey;
    activeData.phi = activephi;
    activeData.speed = activespeed;
    activeData.t = activet;
    activeData.spikes = activespikes;
    activeData.r = activer;
    
    varargout{1} = occMap;
    varargout{2} = spikeMap;
    varargout{3} = infoMap;
    varargout{4} = infoMap_asd;
    varargout{5} = placeMap;
    varargout{6} = placeMap_smooth;
    varargout{7} = placeMap_asd;
    varargout{8} = placeMap_pertrial;
    varargout{9} = normplaceMap_pertrial;
    varargout{10} = PCidx;
    varargout{11} = downData;
    varargout{12} = activeData;
end



