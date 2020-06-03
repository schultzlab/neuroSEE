function [ laps_per_file, phi_trials ] = frun_getFileLaps( list, force )
    if nargin<2, force = false; end

% Load module folders and define data directory
[data_locn, ~, err] = load_neuroSEEmodules(false);
if ~isempty(err)
    beep
    cprintf('Errors',err);    
    return
end

% list
listfile = [data_locn 'Digital_Logbook/lists_imaging/' list];
files = extractFilenamesFromTxtfile( listfile );

laps_per_file = zeros(size(files,1),1);
for n = 1:size(files,1)
    % image file
    file = files(n,:);

    % find tracking file then load it
    trackfile = findMatchingTrackingFile( data_locn, file );
    dir_proc = [data_locn 'Data/' file(1:8) '/Processed/' file '/'];
    if ~exist([dir_proc 'behaviour/' file '_downTrackdata.mat'],'file')
        Trackdata = load_trackfile(data_locn, file, trackfile, force);
        ds = 0;
    else
        M = load([data_locn 'Data/' file(1:8) '/Processed/' file '/behaviour/' file '_downTrackdata.mat']);
        downTrackdata = M.downTrackdata;
        ds = 1;
    end

    % find number of laps per file
    Vthr = 20; Nbins = 50;
    if ds
        data = downTrackdata;
    else
        data = Trackdata;
    end
    activephi = data.phi(data.speed > Vthr);
    bin_phi = discretize(activephi,Nbins);

    % find the delineations per trial (i.e. loop)
    idx_tr = find( bin_phi == bin_phi(1) );
    for k = numel(idx_tr):-1:2
        if (idx_tr(k) - idx_tr(k-1)) <= Nbins 
            idx_tr(k) = 0;
        end
    end
    idx_tr = idx_tr( idx_tr > 0 );
    laps_per_file(n) = numel(idx_tr)-1;
    if laps_per_file(n) > 0
        for k = 1:numel(idx_tr)-1
            phi_trials{n}{k} = activephi(idx_tr(k)+1:idx_tr(k+1));
        end
    else
        phi_trials{n} = [];
    end

end