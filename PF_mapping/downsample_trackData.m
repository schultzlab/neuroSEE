function varargout = downsample_trackData( trackData, spikes, fr )

if iscell(trackData)
    r_all = [];
    Nfiles = numel(spikes);
    for jj = 1:Nfiles
        x = trackData{jj}.x;
        y = trackData{jj}.y;
        r = trackData{jj}.r;
        phi = trackData{jj}.phi;
        speed = trackData{jj}.speed;
        tracktime = trackData{jj}.time;

        % Pre-process tracking data
        t0 = tracktime(1);                  % initial time in tracking data
        Nt(jj) = size(spikes{jj},2);
        
        % Convert -180:180 to 0:360
        if min(phi)<0
           phi(phi<0) = phi(phi<0)+360;
        end

        % generate imaging timestamps using known image frame rate
        dt = 1/fr;
        t = (t0:dt:Nt(jj)*dt)';
        if length(t) ~= Nt(jj)
            t = (t0:dt:(Nt(jj)+1)*dt)';
        end

        % Downsample tracking to Ca trace
        downData{jj}.phi = interp1(tracktime,phi,t,'linear');
        downData{jj}.x = interp1(tracktime,x,t,'linear');
        downData{jj}.y = interp1(tracktime,y,t,'linear');
        downData{jj}.speed = interp1(tracktime,speed,t,'linear'); % mm/s
        downData{jj}.r = interp1(tracktime,r,t,'linear'); % mm/s
        downData{jj}.time = t;
        r_all = [r_all; downData{jj}.r];
    end
    varargout{1} = downData;
    varargout{2} = r_all;
    
else
    x = trackData.x;
    y = trackData.y;
    r = trackData.r;
    phi = trackData.phi;
    speed = trackData.speed;
    tracktime = trackData.time;

    % Pre-process tracking data
    t0 = tracktime(1);                  % initial time in tracking data
    Nt = size(spikes,2);                % number of timestamps for spikes

    % Convert -180:180 to 0:360
%     if min(phi)<0
%        phi(phi<0) = phi(phi<0)+360;
%     end
    phi = phi + 180;

    % generate imaging timestamps using known image frame rate
    dt = 1/fr;
    t = (t0:dt:Nt*dt)';
    if length(t) ~= Nt
        t = (t0:dt:(Nt+1)*dt)';
    end

    % Downsample tracking to Ca trace
    downData.phi = interp1(tracktime,phi,t,'nearest');
    downData.x = interp1(tracktime,x,t,'nearest');
    downData.y = interp1(tracktime,y,t,'nearest');
    downData.speed = interp1(tracktime,speed,t,'nearest'); % mm/s
    downData.r = interp1(tracktime,r,t,'nearest'); % mm/s
    downData.time = t;
    
    varargout{1} = downData;
    varargout{2} = [];
end
