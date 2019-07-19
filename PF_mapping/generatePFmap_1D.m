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
%   spikeMap    : spike map (Ncells rows x Nbins columns)
%   infoMap     : information map 
%   place field maps
%     pfMap     : place field map obtained with histogram estimation 
%     pfMap_sm  : smoothed version of placeMap 
%     pfMap_asd     : place field map obtained with ASD
%     normpfMap     : place field map obtained with histogram estimation 
%     normpfMap_sm  : smoothed version of placeMap 
%     normpfMap_asd : place field map obtained with ASD
%     pfMap_pertrial
%     normpfMap_pertrial
%   downData    : tracking data downsampled to imaging frequency, fields are
%                 x, y, r, phi, speed, t
%   activeData  : downsampled tracking data for when animal was moving, fields are
%                 x, y, r, phi, speed, t, spikes, spikes_pc 
% 
%   varargout: [ occMap, spikeMap, infoMap, infoMap_asd, pfMap, pfMap_sm, pfMap_asd, ...
%                normpfMap, normpfMap_smooth, normpfMap_asd, pfMap_pertrial, normpfMap_pertrial, ...
%                pcIdx, downData, activeData ]
% 

function varargout = generatePFmap_1d( spikes, imtime, trackData, params )
    
fr = params.fr;
Nbins = params.PFmap.Nbins;
Nepochs = params.PFmap.Nepochs;
Vthr = params.PFmap.Vthr;
histsmoothFac = params.PFmap.histsmoothFac;
Ncells = size(spikes,1);

%% Pre-process tracking data
tracktime = trackData.time;
x = trackData.x;
y = trackData.y;
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


%% Identify place cells by first calculating PF maps for entire session
% (i.e. Nepochs = 1)

% Initialise matrices
spikeMap = zeros(Ncells, Nbins);         % spike map
pfMap = zeros(Ncells, Nbins);            % place field map
pfMap_asd = zeros(Ncells, Nbins);        % place field map for asd
infoMap = zeros(Ncells, 2);              % mutual info
infoMap_asd = zeros(Ncells, 2);          % mutual info for asd

% Calculate PF maps
occMap = histcounts(bin_phi,Nbins);
for id = 1:Ncells
    z = activespikes(id,:);

    % Spike rate maps
    for i = 1:Nbins
        spikeMap(id,i) = sum(z(bin_phi == i));
    end

    % histogram estimation
    pfMap(id,:) = spikeMap(id,:)./occMap;
    [infoMap(id,1), infoMap(id,2)] = infoMeasures(pfMap(id,:),occMap,0);

    % ASD estimation
    [pfMap_asd(id,:),~] = runASD_1d(bin_phi,z',Nbins);
    [infoMap_asd(id,:), infoMap_asd(id,:)] =...
        infoMeasures(squeeze(pfMap_asd(id,:)),ones(Nbins,1),0);
end

% Identify place cells
[pcIdx,pcIdx_asd] = filter_nonPC(bin_phi, activespikes, infoMap, infoMap_asd, Nbins);
activespikes_pc = activespikes(pcIdx,:);
activespikes_pcasd = activespikes(pcIdxasd,:);
Npcs = length(pcIdx);
Npcs_asd = length(pcIdx_asd);


%% Finalise place field maps, recalculate if Nepochs > 1
if Nepochs == 1
    spikeMap 
    


    % Initialise matrices
    occMap = zeros(Nepochs, Nbins);                 % occupancy map
    spikeMap = zeros(Npcs, Nbins, Nepochs);         % spike map
    pfMap = zeros(Nepochs, Nbins, Npcs);            % place field map
    pfMap_sm = zeros(Nepochs, Nbins, Npcs);         % smoothened place field map
    pfMap_asd = zeros(Nepochs, Nbins, Npcs);        % place field map for asd
    normpfMap = zeros(Nepochs, Nbins, Npcs);        % normalised place field map
    normpfMap_sm = zeros(Nepochs, Nbins, Npcs);     % normalised smoothened place field map
    normpfMap_asd = zeros(Nepochs, Nbins, Npcs);    % normalised place field map for asd
    infoMap = zeros(Nepochs, Npcs, 2);              % mutual info
    infoMap_asd = zeros(Nepochs, Npcs, 2);          % mutual info for asd

    % Calculate PF maps
    for id = 1:Npcs
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
            pfMap(e,:,id) = spikeMap(id,:,e)./occMap(e,:);
            pfMap(isnan(pfMap)) = 0;
            pfMap_sm(e,:,id) = smoothdata(pfMap(e,:,id),'gaussian',Nbins/histsmoothFac);

            normpfMap(e,:,id) = pfMap(e,:,id)./max(pfMap(e,:,id));
            normpfMap_sm(e,:,id) = pfMap_sm(e,:,id)./max(pfMap_sm(e,:,id));
            [infoMap(e,id,1), infoMap(e,id,2)] = infoMeasures(pfMap(e,:,id),occMap(e,:),0);

            % ASD estimation
            [pfMap_asd(e,:,id),~] = runASD_1d(bin_phi_e,(spike_e)',Nbins);
            normpfMap_asd(e,:,id) = pfMap_asd(e,:,id)./max(pfMap_asd(e,:,id));
            [infoMap_asd(e,id,1), infoMap_asd(e,id,2)] =...
                infoMeasures(squeeze(pfMap_asd(e,:,id))',ones(Nbins,1),0);

        end
    end


% Permute so we get Npcs rows, Nbins cols, Nepochs depth
pfMap = permute(pfMap,[3 2 1]);
infoMap = permute(infoMap,[2 3 1]);

pfMap_asd = permute(pfMap_asd,[3 2 1]);
infoMap_asd = permute(infoMap_asd,[2 3 1]);

% Calculate PF maps per trial
dthr = 190;
phi_bound = find( abs(diff(activephi)) > dthr );
Ntrials = length(phi_bound)-1;

spikeMap_pertrial = zeros(length(phi_bound)-1,Nbins,Npcs);
normspikeMap_pertrial = zeros(length(phi_bound)-1,Nbins,Npcs);

for i = 1:Npcs
    z = activespikes_pc(i,:);
    for j = 1:Ntrials
        phi = bin_phi(phi_bound(j):phi_bound(j+1)-1);
        spike = z(phi_bound(j):phi_bound(j+1)-1);

        for k = 1:Nbins
            spikeMap_pertrial(j,k,i) = sum(spike(phi == k));
        end

        normspikeMap_pertrial(j,:,i) = spikeMap_pertrial(j,:,i)./max(spikeMap_pertrial(j,:,i));
        normspikeMap_pertrial(isnan(normspikeMap_pertrial)) = 0;
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
activeData.spikes = activespikes;
activeData.spikes_pc = activespikes_pc; 

varargout{1} = occMap;
varargout{2} = spikeMap;
varargout{3} = infoMap;
varargout{4} = infoMap_asd;
varargout{5} = pfMap;
varargout{6} = pfMap_sm;
varargout{7} = pfMap_asd;
varargout{8} = normpfMap;
varargout{9} = normpfMap_sm;
varargout{10} = normpfMap_asd;
varargout{11} = spikeMap_pertrial;
varargout{12} = normspikeMap_pertrial;
varargout{13} = pcIdx;
varargout{14} = downData;
varargout{15} = activeData;

end


