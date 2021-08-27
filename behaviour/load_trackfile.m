% Written by Ann Go
%
% Loads the tracking file depending on its type
% INPUTS
%   fname_track : tracking file name
% OUPUT
% trackdata is a structure with fields
%   x : mouse's x position                                (mm)
%   y : mouse's y position                                (mm)
%   r : mouse's radial position                           (mm)
%   phi: mouse's angular position                         (degree)
%   time : time from start of recording                   (s)
%   speed : mouse speed                                   (mm/s)
%   TTLout : TTL output signal
%   alpha 
%   w 

function trackdata = load_trackfile( data_locn, file, fname_track, force, BT )
if nargin<5, BT = false; end

    str = sprintf('%s: Loading tracking data\n', file);
    cprintf(str)
    
    if ~BT
        dir_processed = [data_locn 'Data/' file(1:8) '/Processed/' file '/behaviour/'];
        if ~exist(dir_processed,'dir'), mkdir(dir_processed); fileattrib(dir_processed,'+w','g','s'); end
    else
        timestamp  = extractTimeFromFilename( file, 'Track_yyyy-mm-dd-HH-MM-SS' );
        dir_processed = [data_locn 'Data/' datestr(timestamp,'yyyymmdd') '/Neurotar/' file '/'];
    end
    
    
    [~,~,ext] = fileparts(fname_track);
    fname_fig = [dir_processed file '_mtrajectory.fig'];
    fname_png = [dir_processed file '_mtrajectory.png'];
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
            if isfield(trackdata,'zone')
                trackdata.zone = data.zone;
            end
        case('.csv') % 2018 files
            save_fname = [dir_processed fname_track(end-33:end-4) '.mat'];
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
            zone   = channelData{10}(2:end);   trackdata.zone   = zone;
            TTLout = channelData{12}(2:end);   trackdata.TTLout = TTLout;
            try
                save_fname = [dir_processed fname_track(end-34:end-5) '.mat'];
                save(save_fname, 'time','r','phi','a lpha','x','y','w','speed','TTLout');
            catch
                save_fname = [dir_processed fname_track(end-29:end-5) '.mat'];
                save(save_fname, 'time','r','phi','alpha','x','y','w','speed','TTLout');
            end
    end
    
    yn_fname_fig = exist(fname_fig,'file');
    yn_fname_png = exist(fname_png,'file');
    if any([ force, ~yn_fname_fig, ~yn_fname_png ])
        fig = figure; plot(trackdata.x,trackdata.y); axis square; axis off; 
        tmax = (trackdata.time(end))/60; % min
            if ~BT
                titletext = [file(1:8) '-' file(10:11) '.' file(13:14) '.' file(16:17) ' (' num2str(round(tmax)) ' min)'];
            else
                titletext = [datestr(timestamp,'yyyymmdd HH:MM:SS') ' (' num2str(round(tmax)) ' min)'];
            end
            title(titletext);
        savefig(fig,fname_fig);
        saveas(fig,fname_png);
        close(fig);
    end
    
    newstr = sprintf( '%s: Tracking data loaded\n', file );
    refreshdisp( newstr, str );
end
