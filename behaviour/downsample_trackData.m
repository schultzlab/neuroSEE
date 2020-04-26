function varargout = downsample_trackData( trackData, Nt, fr, imtime )

if nargin<4, imtime = []; end
Nfiles = numel(trackData); 

if iscell(trackData)
    r_all = [];
    for jj = 1:Nfiles
        x = trackData{jj}.x;
        y = trackData{jj}.y;
        r = trackData{jj}.r;
        phi = trackData{jj}.phi;
        speed = trackData{jj}.speed;
        tracktime = trackData{jj}.time;

        % Pre-process tracking data
        t0 = tracktime(1);                  % initial time in tracking data

        % Convert -180:180 to 0:360
        if min(phi)<0
            phi(phi<0) = phi(phi<0)+360;
        end
%         phi = phi + 180;

        % generate imaging timestamps using known image frame rate
        if isempty(imtime)
            dt = 1/fr;
            imtime = t0+(0:dt:(Nt-1)*dt)';
        end

        % Downsample tracking to Ca trace
        [tracktime, ind] = unique(tracktime); 
        downData{jj}.phi = interp1( tracktime, phi(ind), imtime, 'nearest' );
        downData{jj}.x = interp1( tracktime, x(ind), imtime, 'nearest' );
        downData{jj}.y = interp1( tracktime, y(ind), imtime, 'nearest' );
        downData{jj}.speed = interp1( tracktime, speed(ind), imtime, 'nearest' ); % mm/s
        downData{jj}.r = interp1( tracktime, r(ind), imtime, 'nearest' ); % mm/s
        downData{jj}.time = imtime;
        r_all = [r_all; downData{jj}.r];
        clear ind
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

    % Convert -180:180 to 0:360
    if min(phi)<0
       phi(phi<0) = phi(phi<0)+360;
    end
%     phi = phi + 180;

    % generate imaging timestamps using known image frame rate
    if isempty(imtime)
        dt = 1/fr;
        imtime = t0+(0:dt:(Nt-1)*dt)';
    end

    % Downsample tracking to Ca trace
    [tracktime, ind] = unique(tracktime); 
    downData.phi = interp1( tracktime, phi(ind), imtime, 'nearest');
    downData.x = interp1( tracktime, x(ind), imtime, 'nearest');
    downData.y = interp1( tracktime, y(ind), imtime, 'nearest');
    downData.speed = interp1( tracktime, speed(ind), imtime, 'nearest'); % mm/s
    downData.r = interp1( tracktime, r(ind), imtime, 'nearest'); % mm/s
    downData.time = imtime;
    % downData.time = interp1( tracktime, tracktime, imtime, 'nearest' );
    
    varargout{1} = downData;
    varargout{2} = [];
end
