% Written by Ann Go
% Script which takes as input a list of Y-maze experiment files and outputs
% the following per file:
% 1) # of trials
% 2) # of correct trials
% 3) # of left turns
% 4) # of right turns
% 5) # of correct left turns
% 6) # of correct right turns
% 7) the time indices for different trials
% 8) direction and accuracy of each trial

% Flag force to load tracking data from raw file instead of mat file.

function [Ntrials, Ncorrecttrials, Nleftturns, Nrightturns, Ncorrectleft, Ncorrectright, Ttrials, trials_info] = frun_getFileYTrials( list, force )
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

% initialise arrays
Ntrials = zeros(size(files,1),1);
Ncorrecttrials = zeros(size(files,1),1);
Nleftturns = zeros(size(files,1),1);
Nrightturns = zeros(size(files,1),1);
Ncorrectleft = zeros(size(files,1),1);
Ncorrectright = zeros(size(files,1),1);
Ttrials = cell(size(files,1),1);
trials_info = cell(size(files,1),1);

for n = 1:size(files,1)
    % image file
    file = files(n,:);

    % find tracking file, load it, downsample and save it
    trackfile = findMatchingTrackingFile( data_locn, file, force );
    dir_proc = [data_locn 'Data/' file(1:8) '/Processed/' file '/'];
    if or(force, ~exist([dir_proc 'behaviour/' file '_downTrackdata.mat'],'file')) 
        Trackdata = load_trackfile(data_locn, file, trackfile, force);
        downTrackdata = downsample_trackData( Trackdata, Nt, fr, [] );
    else
        M = load([data_locn 'Data/' file(1:8) '/Processed/' file '/behaviour/' file '_downTrackdata.mat']);
        downTrackdata = M.downTrackdata;
        if ~isfield(downTrackdata,'zone')
            trackfile = findMatchingTrackingFile( data_locn, file, true );
            Trackdata = load_trackfile(data_locn, file, trackfile, true);
            downTrackdata = downsample_trackData( Trackdata, Nt, fr, [] );
        end
    end

    % find the delineations for diff trials: find the times the animal
    % returns to zone 9
    zone = downTrackdata.zone;
    tzone9 = find(zone==9);
    idx_tr = find(diff(tzone9)>20)+1;
    Ntrials(n) = numel(idx_tr);
    
    % if the maximum phi for the trial is > 200, the animal turned right
    clear turns;
    for j = 1:numel(idx_tr)-1
        if j == 1
            if mean(downTrackdata.phi(tzone9(1):tzone9(idx_tr(1))-1)) > 200
                turns(1) = 'r';
            else
                turns(1) = 'l';
            end
            Ttrials{n}.trials{1} = tzone9(1):tzone9(idx_tr(1))-1;
        else
            if mean(downTrackdata.phi(tzone9(idx_tr(j-1)):tzone9(idx_tr(j))-1)) > 200
                turns(j) = 'r';
            else
                turns(j) = 'l';
            end
            Ttrials{n}.trials{j} = tzone9(idx_tr(j-1)):tzone9(idx_tr(j))-1;
        end
    end
    if mean(downTrackdata.phi(tzone9(idx_tr(j)):tzone9(idx_tr(end))-1)) > 200
        turns(Ntrials(n)) = 'r';
    else
        turns(Ntrials(n)) = 'l';
    end
    Ttrials{n}.trials{Ntrials(n)} = tzone9(idx_tr(j)):tzone9(idx_tr(end))-1;
    
    % number of left and right turns
    Nleftturns(n) = numel(strfind(turns,'l'));
    Nrightturns(n) = numel(strfind(turns,'r'));
    
    % find number of correct left and right turns
    % consider the first turn as correct
    clear accuracy;
    Ncorrecttrials(n) = 1;
    accuracy(1) = 1;
    if strcmp(turns(1),'l')
        Ncorrectleft(n) = 1; Ncorrectright(n) = 0;
    else
        Ncorrectright(n) = 1; Ncorrectleft(n) = 0;
    end
    % succeeding turns are correct if they are different from previous turn
    for j = 2:Ntrials(n)
        if ~strcmp(turns(j),turns(j-1))
            Ncorrecttrials(n) = Ncorrecttrials(n) + 1;
            accuracy(j) = 1; % correct turn
            if strcmp(turns(j),'l')
                Ncorrectleft(n) = Ncorrectleft(n) + 1;
            else
                Ncorrectright(n) = Ncorrectright(n) + 1;
            end
        else
            accuracy(j) = 0; % wrong turn
        end
    end

    trials_info{n}.direction = turns;
    trials_info{n}.accuracy = accuracy;
end

