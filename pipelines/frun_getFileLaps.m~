function [ laps_per_file, phi_trials ] = frun_getFileLaps( list, force )
if nargin<2, force = false; end
Nt = 7420;
params = neuroSEE_setparams;
fr = params.PFmap.fr;

% Load module folders and define data directory
[data_locn, ~, err] = load_neuroSEEmodules;
if ~isempty(err)
    beep
    cprintf('Errors',err);    
    return
end

% list
listfile = [data_locn 'Digital_Logbook/lists_imaging/' list];
files = extractFilenamesFromTxtfile( listfile );

laps_per_file = zeros(size(files,1),1);
phi_trials = cell(size(files,1),1);
for n = 1:size(files,1)
    % image file
    file = files(n,:);

    % find tracking file then load it
    trackfile = findMatchingTrackingFile( data_locn, file );
    dir_proc = [data_locn 'Data/' file(1:8) '/Processed/' file '/'];
    if ~exist([dir_proc 'behaviour/' file '_downTrackdata.mat'],'file')
        Trackdata = load_trackfile(data_locn, file, trackfile, force);
        downTrackdata = downsample_trackData( Trackdata, Nt, fr, [] );
    else
        M = load([data_locn 'Data/' file(1:8) '/Processed/' file '/behaviour/' file '_downTrackdata.mat']);
        downTrackdata = M.downTrackdata;
    end

    % find number of laps per file
    params = neuroSEE_setparams;
    Vthr = params.PFmap.Vthr;
    Nbins = params.PFmap.Nbins_1D;
    activephi = downTrackdata.phi(downTrackdata.speed > Vthr);
    bin_phi = discretize(activephi,Nbins);

    % find the delineations for diff trials (laps)
    idx_tr = [];
    for li = 1:numel(bin_phi)
        if bin_phi(li) <= bin_phi(1)+1 && bin_phi(li) >= bin_phi(1)-1
            idx_tr = [idx_tr; li];
        end
    end
    % Eliminate false delineations - if indices are less than Nbins
    % apart, they are probably due to animal staying in the same place
    % or moving back and forth (if animal moves Nbins (102 cm) in Nbins
    % time points (30.9 Hz recording), this corresponds to >600 mm/s!
    % Animal doesn't run this fast, so Nbins is a safe threshold.
    temp = idx_tr;
    for l = 1:numel(idx_tr)-1
        if (temp(l+1) - temp(l)) <= Nbins 
            idx_tr(l+1) = 0;
        end
    end
    idx_tr = idx_tr( idx_tr > 0 );
    idx1 = []; idx2 = []; idx3 = []; idx4 = [];
    temp = idx_tr;
    for l = 1:numel(idx_tr)-1
        bp_tr = bin_phi(temp(l):temp(l+1)); 
        for li = 1:numel(bp_tr)
            % This trial must contain a bin close to the 360th degree position
            if bp_tr(li) <= Nbins && bp_tr(li) >= Nbins-4
                idx1 = [idx1; li];
            end
            % This trial must contain a bin close to the 270th degree position
            if bp_tr(li) <= (3*Nbins/4)+2 && bp_tr(li) >= (3*Nbins/4)-2
                idx2 = [idx2; li];
            end
            % This trial must contain a bin close to the 180th degree position
            if bp_tr(li) <= (Nbins/2)+2 && bp_tr(li) >= (Nbins/2)-2
                idx3 = [idx3; li];
            end
            % This trial must contain a bin close to the 90th degree position
            if bp_tr(li) <= (Nbins/4)+2 && bp_tr(li) >= (Nbins/4)-2
                idx4 = [idx4; li];
            end
        end
        if any([isempty(idx1) isempty(idx2) isempty(idx3) isempty(idx4)])
            idx_tr(l+1) = 0;
        end
    end
    idx_tr = idx_tr( idx_tr > 0 );
    laps_per_file(n) = numel(idx_tr)-1;
    if laps_per_file(n) > 0
        idx_tr(1) = 0;
        for k = 1:numel(idx_tr)-1
            phi_trials{n}{k} = activephi(idx_tr(k)+1:idx_tr(k+1));
        end
    else
        phi_trials{n} = [];
    end

end