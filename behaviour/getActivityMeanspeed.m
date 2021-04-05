% Written by Ann Go
% This function calculated activity (% of training session) and mean speed
% for the list of training sessions in list_BT.
% Flag 'force' to override previous processing.

function [ activity, meanspeed ] = getActivityMeanspeed( list_BT, Vthr, force )
if nargin<2, Vthr = 10; end %(mm/s) - speed threshold for qualifying presence of activity
if nargin<3, force = 0; end

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
files = extractFilenamesFromTxtfile( listfile );

% convert filenames to timestamps
format = 'Track_yyyy-mm-dd-HH-MM-SS'; 
Nsessions = 1;
for n = 1:size(files,1)
    % find timestamp of tracking file
    file = files(n,:);
    timestamps(n,:)  = extractTimeFromFilename( file, format );
    
    % identify tracking file type
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
                end
            end
            % load tracking file
            if ~isempty(fname_track)
                try
                    trackdata_temp{n} = load_trackfile( data_locn, [], fname_track, force );    
                catch
                    trackdata_temp{n} = [];
                    fprintf('%s: Tracking file does not exist in folder\n', file);
                end
            else
                trackdata_temp{n} = [];
                fprintf('%s: Tracking file does not exist in folder\n', file);
            end
        else
            fname_track = matfile.name;
            trackdata_temp{n} = load_trackfile( data_locn, [], fname_track, force );
        end            
    else
        trackdata_temp{n} = [];
        fprintf('%s: Tracking file does not exist in folder.', file);
    end
    
    if n>1
        % concatenate files for which timestamps are less than an hour apart
        if etime(timestamps(n,:), timestamps(n-1,:)) < 60*60
            trackdata{Nsessions}.time   = [trackdata{Nsessions}.time    trackdata{Nsessions}.time(end) + trackdata_temp{n}.time];
            trackdata{Nsessions}.r      = [trackdata{Nsessions}.r       trackdata_temp{n}.r];
            trackdata{Nsessions}.phi    = [trackdata{Nsessions}.phi     trackdata_temp{n}.phi];
            trackdata{Nsessions}.alpha  = [trackdata{Nsessions}.alpha   trackdata_temp{n}.alpha];
            trackdata{Nsessions}.x      = [trackdata{Nsessions}.x       trackdata_temp{n}.x];
            trackdata{Nsessions}.y      = [trackdata{Nsessions}.y       trackdata_temp{n}.y];
            trackdata{Nsessions}.w      = [trackdata{Nsessions}.w       trackdata_temp{n}.w];
            trackdata{Nsessions}.TTLout = [trackdata{Nsessions}.TTLout  trackdata_temp{n}.TTLout];
            trackdata{Nsessions}.speed  = [trackdata{Nsessions}.speed   trackdata_temp{n}.speed];
        else
            Nsessions = Nsessions + 1;
            trackdata{Nsessions} = trackdata_temp{n};
        end
    else
        trackdata{Nsessions} = trackdata_temp{n};
    end
end

% calculate activity and mean speed for 45 min
for s = 1:Nsessions
    % only use 45 min (relevant for concatenated files)
    if ~isempty(trackdata{s})
        ind = find(trackdata{s}.time>45*60);
        if ~isempty(ind)
            k = ind(1);
            speed = trackdata{s}.speed(1:k);
            activity(s) = numel(find(speed>Vthr))/k;
            meanspeed(s) = mean(trackdata{s}.speed(1:k));
        else
            speed = trackdata{s}.speed;
            activity(s) = numel(find(speed>Vthr))/length(speed);
            meanspeed(s) = mean(trackdata{s}.speed);
        end
    else
        % if file doesn't exist, assign a NaN value to activity and meanspeed
        activity(s) = NaN;
        meanspeed(s) = NaN;
    end
end
   
activity = activity';
meanspeed = meanspeed';

% save activity and meanspeed data
sdir = [data_locn 'Analysis/' mouseid '/training_data/'];
    if ~exist(sdir,'dir'), mkdir(sdir); end
sname_mat = [sdir mouseid '_trainingdata.mat'];
save(sname_mat, 'trackdata', 'activity', 'meanspeed');

fh1 = figure;
if any(isnan(activity))
    dotitle = 1;
else
    dotitle = 0;
end
activity(isnan(activity)) = -0.1;
plot(1:Nsessions,100*activity,'-bo'); 
axis([1 Nsessions -10 100]);
ylabel('Activity (%)'); xlabel('Session no.'); 
if dotitle, title('(Negative values correspond to missing files. Check logbook for notes.)'); end
saveas(fh1, [sdir mouseid '_activity'], 'fig');
saveas(fh1, [sdir mouseid '_activity'], 'png');

fh2 = figure;
if any(isnan(meanspeed))
    dotitle = 1;
else
    dotitle = 0;
end
meanspeed(isnan(meanspeed)) = -5;
plot(1:Nsessions,meanspeed,'-bo'); 
axis([1 Nsessions -5 100]);
ylabel('Activity (%)'); xlabel('Session no.'); 
if dotitle, title('(Negative values correspond to missing files.)'); end
saveas(fh2, [sdir mouseid '_meanspeed'], 'fig');
saveas(fh2, [sdir mouseid '_meanspeed'], 'png');

end