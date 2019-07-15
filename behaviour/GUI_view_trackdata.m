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

% [fname_track, pathname, proceed] = uigetfile('*.csv;*.tdms', 'Choose tracking data to view');
% 
% if proceed
%     try
%         trackdata = load_trackfile(pathname,fname_track);
%     catch
%         beep;
%         errordlg('Not a valid data file!');
%         return
%     end
% end
pathname = '/Volumes/thefarm2/live/CrazyEights/AD_2PCa/Data/20190401/SavedTrack_2019-04-01-11-13-43/';
filename = 'SavedTrack-2019-04-01-11-13-43.tdms';
fname_track = fullfile(pathname,filename);
trackdata = load_trackfile(fname_track); 
    
    fname_fig = [fname_track(1:end-5) '.fig'];
    fig = figure; plot(trackdata.x,trackdata.y); axis square; axis off; 
        titletext = [filename(1:end-5) ' : ' num2str(round(trackdata.time(end)/60,2)) ' min'];
        title(titletext);
    savefig(fig,fname_fig);
    saveas(fig,fname_fig(1:end-5),'pdf');
    


function trackdata = load_trackfile(fname_track)
    [~,~,ext] = fileparts(fname_track);
    switch ext
%         case('.mat')
%             data = load(fname_track);
%             trackdata.time = data.time;
%             trackdata.r = data.r;
%             trackdata.phi = data.phi;
%             trackdata.alpha = data.alpha;
%             trackdata.x = data.x;
%             trackdata.y = data.y;
%             trackdata.w = data.w;
%             trackdata.TTLout = data.TTLout;
%             trackdata.speed = data.speed;
        case('.csv') % 2018 files
            trackdata = csv_import(fname_track);
        case('.tdms') % 2019 files
            [channelData,~] = TDMS_readChannelOrGroup(fname_track,'Data',[]);
            
            reltime = datetime(channelData{2},'Format','HH:mm:ss.SSS');
                [H,M,S] = hms(reltime(2:end));
                trackdata.time = S + M*60 + H*60*60;

            r      = channelData{3};    trackdata.r      = r(2:end);
            phi    = channelData{4};    trackdata.phi    = phi(2:end);
            alpha  = channelData{5};    trackdata.alpha  = alpha(2:end);
            x      = channelData{6};    trackdata.x      = x(2:end);
            y      = channelData{7};    trackdata.y      = y(2:end);
            w      = channelData{8};    trackdata.w      = w(2:end);
            speed  = channelData{9};    trackdata.speed  = speed(2:end);
            TTLout = channelData{12};   trackdata.TTLout = TTLout(2:end);
    end
    

end
