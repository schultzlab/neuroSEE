% Written by Ann Go (adapted from Giuseppe's PF_ASD_1d.m)
%
% This function maps place fields
%
%   INPUTS:
%       spikes      : non-negative derivative spike estimates
%       imtime      : imaging timestamps
%       trackdata   : cell of tracking data with fields x, y, phi, speed, time
%       params.
%           mode_dim    : 1 for 1D, 2 for 2D 
%           mode_method : 1 for ASD, 2 for histogram estimation
%           Nbins       : number of location bins
%           Nepochs     : number of epochs for each 4 min video
%           histsmoothFac : Gaussian smoothing window for histogram estimation
%           Vthr        : speed threshold (mm/s) Note: David Dupret uses 20 mm/s, Neurotar uses 8 mm/s

%   OUTPUTS:
%       occMap      : occupancy map, 1D: Ncells rows x Nbins columns
%                                    2D: Nbins rows, Nbins, columns, Ncells stack
%       spikesMap   : spike map (same size as occMap)
%       infoMap     : information map (same size as occMap)
%       placeMap    : place field map (same size as occMap)
%       downData    : tracking data downsampled to imaging frequency, fields are
%                       x, y, phi, speed, time
%       activaData  : downsampled tracking data for when animal was moving, fields are
%                       x, y, phi, speed, time, spikes
%       placeMap_smooth : smoothed version of placeMap (for when histogram estimation is used)
% 
%   varargout: [ occMap, spikeMap, infoMap, placeMap, downData, activeData, placeMap_smooth ]

function varargout = generatePFmap( spikes, imtime, trackdata, imrate, Vthr, mode_dim, mode_method, Nbins, Nepochs, histsmoothFac )

    % Read input variables
    tracktime = trackdata.time;
    x = trackdata.x;
    y = trackdata.y;
    r = trackdata.r;
    phi = trackdata.phi;
    speed = trackdata.speed;
    
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
       dt = 1/imrate;
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

    if mode_dim == 1 % 1D
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
                
                % ASD estimation
                if mode_method == 1
                    [placeMap(e,:,id),~] = runASD_1d(bin_phi_e,(spike_e)',Nbins);
                    [infoMap(e,id,1), infoMap(e,id,2)] =...
                    infoMeasures(squeeze(placeMap(e,:,id))',ones(Nbins,1),0);
                else % histogram estimation
                    placeMap(e,:,id) = spikeMap(id,:,e)./occMap(e,:);
                    placeMap(isnan(placeMap)) = 0;
                    [infoMap(e,id,1), infoMap(e,id,2)] = infoMeasures(placeMap(e,:,id),occMap(e,:),0);
                end
            end
        end
        placeMap = permute(placeMap,[3 2 1]);
        infoMap = permute(infoMap,[2 3 1]);
        placeMap_smooth = smoothdata(placeMap,2,'gaussian',Nbins/histsmoothFac);
        if mode_method == 1
            placeMap_smooth = [];
        end

    elseif mode_dim == 2 % 2D


    end
   
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
    varargout{4} = placeMap;
    varargout{5} = downData;
    varargout{6} = activeData;
    varargout{7} = placeMap_smooth;
end



