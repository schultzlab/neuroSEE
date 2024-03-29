% Written by Ann Go
% This function calculates activity (% of training session) and mean speed
% for the list of training sessions in list_BT.
% Flag 'force' to override previous processing.

function [ activity, meanspeed ] = getActivityMeanspeed( list_BT, Vthr, force, fsave, figclose )
if nargin<2, Vthr = 13; end %(mm/s) - speed threshold for qualifying presence of activity
if nargin<3, force = false; end
if nargin<4, fsave = true; end
if nargin<5, figclose = true; end

% Load module folders and define data directory
[data_locn, ~, err] = load_neuroSEEmodules;
if ~isempty(err)
    beep
    cprintf('Errors',err);    
    return
end

% list
[mouseid,~] = find_mouseIDexpname( list_BT );
listfile = [data_locn 'Digital_Logbook/lists_BT/' list_BT];
fileID = fopen(listfile);
lines = textscan(fileID,'%s','delimiter','\n');
fclose(fileID);
lines = lines{1}; 

% check if data have been processed
sname_mat = [data_locn 'Analysis/' mouseid '/training_data/' mouseid '_trainingdata.mat'];
if exist(sname_mat,'file') && ~force
    M = load(sname_mat);
    activity = M.activity;
    meanspeed = M.meanspeed;
    fprintf('%s: Processed training data found and loaded.\n', mouseid);
    return
end

% convert filenames to timestamps
format = 'Track_yyyy-mm-dd-HH-MM-SS'; 
Nsessions = 1;
for n = 1:size(lines,1)
    % find timestamp of tracking file
    file = lines{n};
    if length(file) > 6
        timestamps(n,:)  = extractTimeFromFilename( file, format );

        % load tracking file
        trackfolder = [data_locn 'Data/' datestr(timestamps(n,:),'yyyymmdd') '/Neurotar/' file '/'];
        if exist(trackfolder,'dir')
            matfile = subdir(fullfile(trackfolder,['*.','mat']));
            if isempty(matfile) 
                csvfile = subdir(fullfile(trackfolder,['*.','csv'])); % Files from 2018
                if numel(csvfile) == 1
                    fname_track = csvfile.name;
                else 
                    tdmsfile = subdir(fullfile(trackfolder,['*.','tdms'])); % Files from 2019
                    if numel(tdmsfile) == 1
                        fname_track = tdmsfile.name;
                    else
                        fprintf('%s: Tracking file does not exist in folder\n', file);
                        fname_track = '';
                        trackdata_temp{n} = [];
                    end
                end
                if ~isempty(fname_track)
                    try
                        trackdata_temp{n} = load_trackfile( data_locn, file, fname_track, force, true );    
                    catch
                        trackdata_temp{n} = [];
                        fprintf('%s: Error reading tracking file\n', file);
                    end
                end
            else
                fname_track = matfile.name;
                trackdata_temp{n} = load_trackfile( data_locn, file, fname_track, force, true );
            end            
        else
            trackdata_temp{n} = [];
            fprintf('%s: Tracking folder does not exist.\n', file);
        end

        if n>1
            % concatenate files for which timestamps are less than an hour apart
            if all( timestamps(n-1,:) ~= 0 ) && ( etime(timestamps(n,:), timestamps(n-1,:)) < 60*60 )
                if ~isempty(trackdata_temp{n})
                    if size(trackdata_temp{n}.time,2) > size(trackdata_temp{n}.time,1)
                        trackdata{Nsessions}.time   = [trackdata{Nsessions}.time';    trackdata{Nsessions}.time(end) + trackdata_temp{n}.time'];
                        trackdata{Nsessions}.r      = [trackdata{Nsessions}.r';       trackdata_temp{n}.r'];
                        trackdata{Nsessions}.phi    = [trackdata{Nsessions}.phi';     trackdata_temp{n}.phi'];
                        trackdata{Nsessions}.alpha  = [trackdata{Nsessions}.alpha';   trackdata_temp{n}.alpha'];
                        trackdata{Nsessions}.x      = [trackdata{Nsessions}.x';       trackdata_temp{n}.x'];
                        trackdata{Nsessions}.y      = [trackdata{Nsessions}.y';       trackdata_temp{n}.y'];
                        trackdata{Nsessions}.w      = [trackdata{Nsessions}.w';       trackdata_temp{n}.w'];
                        trackdata{Nsessions}.TTLout = [trackdata{Nsessions}.TTLout';  trackdata_temp{n}.TTLout'];
                        trackdata{Nsessions}.speed  = [trackdata{Nsessions}.speed';   trackdata_temp{n}.speed'];
                    else
                        trackdata{Nsessions}.time   = [trackdata{Nsessions}.time;    trackdata{Nsessions}.time(end) + trackdata_temp{n}.time];
                        trackdata{Nsessions}.r      = [trackdata{Nsessions}.r;       trackdata_temp{n}.r];
                        trackdata{Nsessions}.phi    = [trackdata{Nsessions}.phi;     trackdata_temp{n}.phi];
                        trackdata{Nsessions}.alpha  = [trackdata{Nsessions}.alpha;   trackdata_temp{n}.alpha];
                        trackdata{Nsessions}.x      = [trackdata{Nsessions}.x;       trackdata_temp{n}.x];
                        trackdata{Nsessions}.y      = [trackdata{Nsessions}.y;       trackdata_temp{n}.y];
                        trackdata{Nsessions}.w      = [trackdata{Nsessions}.w;       trackdata_temp{n}.w];
                        trackdata{Nsessions}.TTLout = [trackdata{Nsessions}.TTLout;  trackdata_temp{n}.TTLout];
                        trackdata{Nsessions}.speed  = [trackdata{Nsessions}.speed;   trackdata_temp{n}.speed];
                    end
                else
                    Nsessions = Nsessions + 1;
                    trackdata{Nsessions} = [];
                end
            else
                Nsessions = Nsessions + 1;
                trackdata{Nsessions} = trackdata_temp{n};
            end
        else
            trackdata{Nsessions} = trackdata_temp{n};
        end
    else
        Nsessions = Nsessions + 1;
        activity(Nsessions) = str2double(file)/100;
        meanspeed(Nsessions) = str2double(file);
        trackdata{Nsessions} = [];
        timestamps(n,:) = zeros(1,6);
    end
end

% calculate activity and meanspeed
for s = 1:Nsessions
    if ~isempty(trackdata{s})
        ind = find(trackdata{s}.time>45*60);
        if ~isempty(ind)
            k = ind(1);
            speed = trackdata{s}.speed(1:k);
            speed(isinf(speed)) = [];
            activity(s) = numel(find(speed>Vthr))/numel(speed);
            meanspeed(s) = mean(speed);
        else
            speed = trackdata{s}.speed;
            speed(isinf(speed)) = [];
            activity(s) = numel(find(speed>Vthr))/numel(speed);
            meanspeed(s) = mean(speed);
        end
    end
end
                
activity = activity';
meanspeed = meanspeed';

% save activity and meanspeed data
if fsave
    sdir = [data_locn 'Analysis/' mouseid '/training_data/'];
        if ~exist(sdir,'dir'), mkdir(sdir); end
    save([sdir mouseid '_trainingdata.mat'], 'trackdata', 'activity', 'meanspeed');
    save([sdir mouseid '_trainingdata_summary.mat'], 'activity', 'meanspeed');
end

fh1 = figure;
if any(isnan(activity))
    dotitle = 1;
else
    dotitle = 0;
end
% activity(isnan(activity)) = -0.1;
plot(1:Nsessions,100*activity,'-bo'); 
axis([1 Nsessions -10 100]);
ylabel('Activity (%)'); xlabel('Session no.'); 
if dotitle, title('(Negative values correspond to missing files. Check logbook for notes.)'); end
if fsave
    saveas(fh1, [sdir mouseid '_activity'], 'fig');
    saveas(fh1, [sdir mouseid '_activity'], 'png');
end
if figclose, close(fh1); end

fh2 = figure;
if any(isnan(meanspeed))
    dotitle = 1;
else
    dotitle = 0;
end
% meanspeed(isnan(meanspeed)) = -5;
plot(1:Nsessions,meanspeed,'-bo'); 
axis([1 Nsessions -5 100]);
ylabel('Mean speed (%)'); xlabel('Session no.'); 
if dotitle, title('(Negative values correspond to missing files.)'); end
if fsave
    saveas(fh2, [sdir mouseid '_meanspeed'], 'fig');
    saveas(fh2, [sdir mouseid '_meanspeed'], 'png');
end
if figclose, close(fh2); end


end