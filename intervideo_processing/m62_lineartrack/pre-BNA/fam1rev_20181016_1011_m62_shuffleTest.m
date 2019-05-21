clear all; close all;

%% Basic setup

load('fam1rev_20181016_1011_m62.mat');

addpath(genpath('../behaviour'));
addpath(genpath('../intervideo_processing'));
addpath(genpath('../motion_correction'));
addpath(genpath('../PF_mapping'));
addpath(genpath('../ROI_segmentation'));
addpath(genpath('../spike_extraction'));
addpath(genpath('../utilities'));
addpath(genpath('../pipelines'));

%% Declare variables
tracktime = time;
rawphi = phi;

Vthr = 10;
Nbins = 180;
Nepochs = 1;
mode_method = 2;
histsmoothFac = 5;

%% Bootstrap shuffle test
% randind = randperm(numel(rawphi));
% phi = rawphi(randind);

% phi = rawphi;
rawspikes = spikes;
for i = 1:size(rawspikes,1)
    spikes(i,:) = randsample(rawspikes(i,:),size(rawspikes,2));
end

%% Generate PF maps
t0 = tracktime(1);                  % initial time in tracking data
nspikes = spikes; %bsxfun( @rdivide, bsxfun(@minus, spikes, min(spikes,[],2)), range(spikes,2) ); % normalisation
Ncells = size(spikes,1);            % number of cells
Nt = size(spikes,2);                % number of timestamps for spikes

% Convert -180:180 to 0:360
if min(phi)<0
   phi(phi<0) = phi(phi<0)+360;
end

imrate = 30.91;
dt = 1/imrate;
t = (t0:dt:Nt*dt)';


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

%% Sort PF maps

n_info = 1/10; % ratio of units to consider as not informative
info_type = 2; % 1 is info/sec, 2 is info/spk
Ncells = size(placeMap,1);
Nbins = size(placeMap,2);
hdl = figure;

% Initialise arrays
sorted_placeMap = zeros(Ncells,Nbins,Nepochs);
normsorted_placeMap = zeros(Ncells,Nbins,Nepochs);
allIdx = zeros(Ncells,Nepochs);

for e = 1:Nepochs
    mat = squeeze(placeMap_smooth(:,:,e));                                 % Ncells x Nbins
    % sort by information
    [~,idx] = sort( squeeze(infoMap(:,:,e)) );                      % Ncells x 2
    idx_info = idx( round(n_info*Ncells)+1:Ncells, info_type );     % info idx (<Ncells x 1)
    idx_no = idx( 1:round(n_info*Ncells), info_type );              % no info idx (<Ncells x 1)
    pf_info = mat(idx_info,:);                                      % subset of placeMap with info
    maxloc = zeros(length(idx_info),1);
    for i = 1:length(idx_info)
        [~,maxloc(i)] = max(pf_info(i,:));
    end
    [~,sortIdx(:,e)] = sort(maxloc);                                % sort according to phi of max info
    rawsorted_placeMap = [ pf_info(sortIdx(:,e),:); mat(flip(idx_no),:) ]; % cat 2 sorted subsamples

    % normalise for visualisation
    for i = 1:Ncells
        sorted_placeMap(i,:,e) = rawsorted_placeMap(i,:);
        normsorted_placeMap(i,:,e) = rawsorted_placeMap(i,:)/max(rawsorted_placeMap(i,:));
    end

    allIdx(:,e) = [idx_info(sortIdx(:,e)); flip(idx_no)];

    % plot
    if Nepochs > 1
        figure(hdl); 
        subplot(1,Nepochs,e); imagesc(sorted_placeMap(allIdx,:,e)); 
        % colorbar; % caxis([0,0.005]);
        xticks([1 30 60 90 120 150 180]); xticklabels([1 60 120 180 240 300 360]); 
        yticks(1:2:Ncells); 
        yticklabels(allIdx(1:2:Ncells,e)); 
        str = sprintf('Epoch %g', e);
        title(str); 
        xlabel('Position (degrees)'); ylabel('Cell number');
    else
        figure(hdl); subplot(121);
        imagesc(normsorted_placeMap(:,:)); 
        % colorbar; % caxis([0,0.005]);
        % xticks([1 30 60 90 120 150 180]); xticklabels([1 60 120 180 240 300 360]); 
        % yticks(1:2:Ncells); 
        %yticklabels(allIdx); %(1:2:Ncells,e)); 
        title('Normalised place field maps');
        xlabel('Position (degrees)'); ylabel('Cell number');
        subplot(122);
        imagesc(sorted_placeMap(:,:)); 
        colorbar; caxis([0,0.0145]);
        % xticks([1 30 60 90 120 150 180]); xticklabels([1 60 120 180 240 300 360]); 
        % yticks(1:2:Ncells); 
        %yticklabels(allIdx); %(1:2:Ncells,e)); 
        title('Place field maps');
        xlabel('Position (degrees)'); ylabel('Cell number');
    end
end

save('fam1rev_20181016_1011_m62_bootstrapshuffle_hist.mat')