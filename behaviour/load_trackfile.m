% Written by Ann Go
%
% Loads the tracking file depending on its type
% INPUTS
%   file        : image file name (optional, specify as [] if NA)
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

function trackdata = load_trackfile(data_locn, file, fname_track, force)
    
    if ~isempty(file)
        dir_processed = [data_locn 'Data/' file(1:8) '/Processed/' file '/behaviour/'];
        if ~exist(dir_processed,'dir'), mkdir(dir_processed); end
        fname_mat = [dir_processed file '.mat'];
        fname_fig = [dir_processed file '_mtrajectory.fig'];
        fname_png = [dir_processed file '_mtrajectory.png'];
        txt = file;
    else
        format = 'Track_yyyy-mm-dd-HH-MM-SS';
        [filepath,name,ext] = fileparts(fname_track);
        timestamp  = extractTimeFromFilename( name, format );
        fname_mat = [filepath '/' name '.mat'];
        fname_fig = [filepath '/' name '_mtrajectory.fig'];
        fname_png = [filepath '/' name '_mtrajectory.png'];
        txt = name;
    end
    
    str = sprintf('%s: Loading tracking data\n', txt);
    cprintf(str)
        
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
            trackdata = csv_import(fname_track,fname_mat);
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
            save(fname_mat, 'time','r','phi','alpha','x','y','w','speed','TTLout');
    end
    
    yn_fname_fig = exist(fname_fig,'file');
    yn_fname_png = exist(fname_png,'file');
    if any([ force, ~yn_fname_fig, ~yn_fname_png ])
        fig = figure; plot(trackdata.x,trackdata.y); axis square; axis off; 
        tmax = (trackdata.time(end))/60; % min
            if ~isempty(file)
                titletext = [file(1:8) '-' file(10:11) '.' file(13:14) '.' file(16:17) ' (' num2str(round(tmax)) ' min)'];
                title(titletext);
            else
                titletext = datestr(timestamp, 'yyyy/mm/dd HH:MM:SS');
                title(titletext);
            end
        savefig(fig,fname_fig);
        saveas(fig,fname_png);
        close(fig);
    end
    
    newstr = sprintf( '%s: Tracking data loaded\n', txt );
    refreshdisp( newstr, str );
end
