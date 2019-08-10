% Written by Ann Go

function frun_genTrainingPlots( mouseid, list, force )

if nargin<3, force = false; end

%% Load module folders and define data directory
test = false;                   % flag to use one of smaller files in test folder)
[data_locn,~,err] = load_neuroSEEmodules(test);
if ~isempty(err)
    beep
    cprintf('Errors',err);    
    return
end

files = extractFilenamesFromTxtfile( list );

for i = 1:size(files,1)
    file = files(i,:);
    fdate = [file(7:10) file(12:13) file(15:16)];
    filedir = [data_locn 'Data/' fdate '/Neurotar/' file '/'];
    
    % find tracking file
    matfile = subdir(fullfile(filedir,['*.','mat']));
            
    if force || isempty(matfile)
        csvfile = subdir(fullfile(filedir,['*.','csv'])); % Files from 2018
        if numel(csvfile) == 1
            fname_track = csvfile.name;
        else 
            tdmsfile = subdir(fullfile(filedir,['*.','tdms'])); % Files from 2019
            if numel(tdmsfile) == 1
                fname_track = tdmsfile.name;
            else
                fprintf('%s: Tracking file does not exist in folder\n', file);
            end
        end
    else
        fname_track = matfile.name;
    end
    
    % load tracking data
    [~,~,ext] = fileparts(fname_track);
    switch ext
        case('.mat')
            data = load(fname_track);
            trackdata.time = data.time;
            trackdata.r = data.r;
            trackdata.phi = data.phi;
            trackdata.alpha = data.alpha;
            trackdata.x = data.x;
            trackdata.y = data.y;
            trackdata.w = data.w;
            trackdata.TTLout = data.TTLout;
            trackdata.speed = data.speed;
        case('.csv') % 2018 files
            save_fname = [filedir fname_track(end-33:end-4) '.mat'];
            trackdata = csv_import(fname_track,save_fname);
        case('.tdms') % 2019 files            
            [channelData,~] = TDMS_readChannelOrGroup(fname_track,'Data',[]);
            
            reltime = datetime(channelData{2},'Format','HH:mm:ss.SSS');
                [H,M,S] = hms(reltime(2:end));
                time = S + M*60 + H*60*60;     trackdata.time   = time;

            r      = channelData{3}(2:end);    trackdata.r      = r;
            phi    = channelData{4}(2:end);    trackdata.phi    = phi;
            alpha  = channelData{5}(2:end);    trackdata.alpha  = alpha;
            x      = channelData{6}(2:end);    trackdata.x      = x;
            y      = channelData{7}(2:end);    trackdata.y      = y;
            w      = channelData{8}(2:end);    trackdata.w      = w;
            speed  = channelData{9}(2:end);    trackdata.speed  = speed;
            TTLout = channelData{12}(2:end);   trackdata.TTLout = TTLout;
            try
                save_fname = [filedir fname_track(end-34:end-5) '.mat'];
                save(save_fname, 'time','r','phi','alpha','x','y','w','speed','TTLout');
            catch
                save_fname = [filedir fname_track(end-29:end-5) '.mat'];
                save(save_fname, 'time','r','phi','alpha','x','y','w','speed','TTLout');
            end
    end

    M(i).x =  

end