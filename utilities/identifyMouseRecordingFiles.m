% [mouseexp, recordings] = identifyMouseRecordingDates( [serverdir] )
% Extract valid mouse experiments that have a Calcium folder and a tracking
% folder within a few minutes of the Calcium start time, and a valid entry
% in the logbook. 
% Inputs:
%   serverdir  - access to the farm, defaults to '/Volumes'
% Outputs:
%   mouseexp   - cell array with entry for each mouse that has valid data 
%                each cell has a cell array of valid experiments, where each is 
%                a struct containing valid Calcium directory names & tracking
%                directory names 
%   recordings - valid recordings from the log book, with matching calcium
%                and tracking recordings

function mouseDataList = identifyMouseRecordingFiles( data_locn, mousename )

   logbookdir  = fullfile( data_locn, 'Digital_Logbook', 'byAnimal' );
   gcampdir    = fullfile( data_locn, 'Data' );
   dateformat  = 'yyyymmdd'; % using matlab's formatting for datestr
   matchtype   = '[0-9a-zA-Z ]+'; % match dir & filenames based on alphanumeric chars (but not '_')
      
   % months must be capital to separate from minutes
   imgformat   = 'yyyymmdd_HH_MM_SS_2P';
   trackformat = 'Track_yyyy-mm-dd-HH-MM-SS';
   photondir   = '2P';
   trackdirname= 'Neurotar';
   
   mouseDataList.name     = mousename;
   mouseDataList.num      = mousename(2:end);

  % Now go through the log book files for file matching mousename
   lognames  = dir( logbookdir );
   for li=1:length( lognames )
      % we're looking at filenames here
      filename = lognames(li).name;
      % if filename is name of a mouse, proceed
      if strcmp(filename(1:length(mousename)),mousename)
         [~,sheets] = xlsfinfo( fullfile( logbookdir, filename) );
         % cycle through each worksheet name & extract ones with fam1, fam1fam2,
         % fam1nov, fam1fam1rev, fam2, open
         for si=1:length(sheets)
            if strcmp(sheets{si}(1:3),'fam') || strcmp(sheets{si}(1:3),'ope')
               exp = sheets{si};
               num = xlsread(fullfile( logbookdir, filename), sheets{si}, 'B1');
               expdate = datestr( datetime( num, 'ConvertFrom', 'excel' ), 'yyyymmdd' ); 
               clear num

               % go through recording times & extract the recording start times
               % - data starts on row 16 of column A, but ranges have to
               % have an end number too, so just pick something large
               % because the function cuts it off if empty anyways
               % - the date data will be into num rather than text or raw
               num = xlsread( fullfile( logbookdir, filename), sheets{si}, 'A16:A40' );
               [~,~,raw] = xlsread(fullfile( logbookdir, filename), sheets{si}, 'B16:B40');
               % for some reason num can have extra columns - perhaps if
               % columns are merged for some rows? Also it starts with
               % column A even for raw, which is supposed to be column B,
               % because
               if size(num, 2)>1
                  % assumes all columns are merged in at least 1 row so
                  % that columns are read starting from A
                  mergedcells = true; 
                  num = num(1:end,1);
                  raw = raw(1:end,2);
               else
                  mergedcells = false; 
               end
               % get hours, minutes, seconds from num
               [h,m,s] = hms( datetime( num, 'ConvertFrom', 'excel' ) ); 
               invalid = cellfun( @(c) any(isnan(c)), raw );
               valid   = cellfun( @(c) length(c)<4 && c(1)=='X', raw );
               valid   = find( ~invalid & valid );
               h = h(valid); m = m(valid); s = s(valid);

               % now get environment type for each valid row
               [~,~,raw] = xlsread(fullfile( logbookdir, filename), sheets{si}, 'J16:J40');
               if mergedcells
                  raw = raw(1:end,10);
               end
               env = raw( valid );

               % add experiment for current date, including recording times 
               mouseDataList.exp{si}.exp      = exp;
               mouseDataList.exp{si}.date     = expdate;
               mouseDataList.exp{si}.logtimes = [h m s];
               mouseDataList.exp{si}.env      = env;
            end
         end
     end
   end % end for logbook names

   % We need to loop through each experiment, and then look at each time
   % and find the corresponding Ca imaging and Neurotar tracking files
   for ei = 1:numel( mouseDataList.exp )
        logtimes = mouseDataList.exp{ei}.logtimes;
        recorddates = dir( gcampdir );
        for ri = 1:length( recorddates )
            if strcmp(mouseDataList.exp{ei}.date, recorddates{ri}.name)
                % 2 photon directory - will have a list of directory names from times
                % of image recording, e.g. 20181016_10_34_37_2P
                imgdir      = fullfile( gcampdir, recorddates{ri}.name, photondir );
                rawimgtimes    = getTimesFromDirectories( imgdir, imgformat, matchtype );
                % tracking directory - will have a list of directory names from times
                % of track recording, e.g. Track_2018-10-16-10-16-01
                trackdir    = fullfile( gcampdir, dirname, trackdirname );
                rawtracktimes  = getTimesFromDirectories( trackdir, trackformat, matchtype );
                if ~isempty( rawimgtimes ) && ~isempty( rawtracktimes )
                    % get times for this date where we have valid calcium imaging & track data
                    [imgtimes, tracktimes, ok] = matchImgAndTrackTimes(rawimgtimes, rawtracktimes);
                    if ~ok
                        return;
                    end
                    if ~isempty( imgtimes )
                        % some of the recordings will have been aborted or invalid, so
                        % check times with the experiment times for the mouse
                        %  - hour & min of logtime matches img times
                        if length( imgtimes ) >= length( logtimes )
                           [keepCa, keeptime] = ismember( imgtimes(:,1:2), logtimes(:,1:2), 'rows' );
                           Catimes    = Catimes(keepCa,:);
                           tracktimes = tracktimes(keepCa,:);
                           exptimes   = exptimes( keeptime(keeptime>0), : );
                        elseif length( Catimes ) <= length( exptimes )
                           [keeptrack, keeptime] = ismember( exptimes(:,1:2), tracktimes(:,1:2), 'rows' );
                           exptimes   = exptimes(keeptrack,:);
                           exptimes   = exptimes( keeptime(keeptime>0), : );
                        end
                    end
                end
            end
        end % end for each recording date folder
   end % end for each experiment
   

   
   % now you need to compare logbook recordings with mouse data - where
   % they're a match, we're good to go!
   
   % for each mouse, & each experiment for the mouse, find the imaging &
   % tracking files
   recorddates = cell2mat( getStructFieldFromCell(recordings, 'date','cell')' );
   mi = 1; % use while loop cuz removing mice with no valid exp so num mice changes
   while mi<=length(mouseexp)
      mousename = mouseexp{mi};
      mname = mousename.name;
      mnum  = mousename.num;
      numexp= length( mousename.exp );
      ii = 1; % use while loop coz we delete invalid experiments
      while ii<=length(mousename.exp)
         exp     = mousename.exp{ii};
         expdate = exp.date;
         exptimes= exp.times;
         exptype = exp.type;
         
         recind  = find( ismember( recorddates, expdate, 'rows' ) );

         % if there are no valid img/track datasets for this experiment, ditch it
         if isempty( recind ) || isempty( exptimes )
            mousename.exp(ii) = [];
            
         % if there are valid data files for this logged recording, add it
         % to the mouse's video names & track names
         else
            Catimes = recordings{recind}.Catimes;
            tracktimes = recordings{recind}.tracktimes;
            
            % some of the recordings will have been aborted or invalid, so
            % check times with the experiment times for the mouse
            %  - hour & min of exp matches track times
            if length( Catimes ) >= length( exptimes )
               [keepCa, keeptime] = ismember( tracktimes(:,1:2), exptimes(:,1:2), 'rows' );
               Catimes    = Catimes(keepCa,:);
               tracktimes = tracktimes(keepCa,:);
               exptimes   = exptimes( keeptime(keeptime>0), : );
            elseif length( Catimes ) <= length( exptimes )
               [keeptrack, keeptime] = ismember( exptimes(:,1:2), tracktimes(:,1:2), 'rows' );
               exptimes   = exptimes(keeptrack,:);
               exptimes   = exptimes( keeptime(keeptime>0), : );
            end
            
            % at this point we have all calcium and track recording times
            % for all valid experiments for this mouse - get filenames
            if ~isempty( exptimes )
               imgsubdirs    = date2str( expdate, imgformat, Catimes );
               imgdir        = fullfile( gcampdir, recordings{recind}.dirname, photondir );
               imgsubdirs    = [repmat([imgdir filesep], [size(imgsubdirs,1) 1])    imgsubdirs];
               tracksubdirs  = date2str( expdate, trackformat, tracktimes );
               imgdir        = fullfile( gcampdir, recordings{recind}.dirname, trackdirname );
               tracksubdirs  = [repmat([imgdir filesep], [size(tracksubdirs,1) 1])  tracksubdirs];

               exp.times     = exptimes;
               exp.ddir      = fullfile( gcampdir, recordings{recind}.dirname );
               exp.cadirs    = imgsubdirs;
               exp.trackdirs = tracksubdirs;
               mousename.exp{ii} = exp;
               ii = ii+1; % increment to get next experiment
               
            % no valid track & imaging times so ditch exp
            else
               mousename.exp(ii) = []; % shift remaining exp down so don't increment i
            end
            
         end       
      end % end for each experiment for this mouse
      
      % before we save this mouse make sure that there are valid experiments
      if isfield( mousename, 'exp' ) && length( mousename.exp ) > 0
         mouseexp{mi} = mousename;
         mi = mi + 1; % increment to get next mouse
      else
         mouseexp(mi) = []; % don't increment cuz we've shifted mice down
      end
   end % end for each mouse
end

function [finalimgtimes, finaltracktimes, ok] = matchImgAndTrackTimes(imgtimes, tracktimes)
   maxtrackdelay = 11; % track recording must begin within this many mins of Ca
   % now go through image & tracking files, & find the files that
   % have both within a few minutes of each other, else ditch
   % - there could be more of either image or tracking files
   % - as the files have been processed, remove from the lists
   finalimgtimes = []; finaltracktimes = []; % gotta be the same length!!
   ok = true; 
   if any( isnan( imgtimes(:) ) ) || any( isnan( tracktimes(:) ) )
      fprintf( 'img and track times are NaN, exiting...\n' );
      ok = false;
      return;
   end
   
   % remove any track recordings that come before the next img
   % (in case there are additional track recordings we don't want)
   [ih, im] = matsplit( imgtimes(1,:) );   % get img hour mins secs
   [th, tm] = matsplit( tracktimes(1,:) ); % get track hour mins secs
   % if track hour is less, or same hour but mins are less, ditch it coz
   % there's no corresponding calcium recording
   while th<ih || (th==ih && tm<im) || (th==ih && tm==im && tm<im)
      tracktimes(1,:) = [];
      [th, tm] = matsplit( tracktimes(1,:) ); % update track hour mins secs
   end
   maxiter = 1e2;
   numiter = 0;
   while ~isempty( imgtimes )
      numiter = numiter + 1;
      if numiter>=maxiter
         fprintf('stuck in loop trying to match calium & tracking times\n');
         return;
      end
%       [ih, im] = matsplit( imgtimes(1,:) );   % get img hour mins secs
%       [th, tm] = matsplit( tracktimes(1,:) ); % get track hour mins secs
      % check that track recording started within 3 mins of img recording
      if ih==th && im<=tm && tm<=(im+maxtrackdelay) ...
      || ( (ih+1)==th && im>=(60-maxtrackdelay) && tm<maxtrackdelay)
         finalimgtimes(end+1,:)   = imgtimes(1,:);
         finaltracktimes(end+1,:) = tracktimes(1,:);
         imgtimes(1,:) = [];
         tracktimes(1,:) = [];
      end
      if isempty(   imgtimes ), break; end
      if isempty( tracktimes ), break; end
      % now remove any track recordings that come before the next img
      % (in case there are additional track recordings we don't want)
      [ih, im] = matsplit( imgtimes(1,:) );   % get img hour mins secs
      [th, tm] = matsplit( tracktimes(1,:) ); % get track hour mins secs
      % if track hour is less, or same hour but mins are less, ditch it coz
      % there's no corresponding calcium recording
      while th<ih || (th==ih && tm<im) || (th==ih && tm==im && tm<im)
         if isempty( tracktimes), break; end
         tracktimes(1,:) = [];
         if isempty( tracktimes ), break; end
         [th, tm] = matsplit( tracktimes(1,:) ); % update track hour mins secs
      end
      % if there's no track time within maxtrackdelay of this Ca time, then
      % we need to remove it from the list, otherwise we get stuck on it
      % - gotta account for circularity of time!
      if ih==th && tm>(im + maxtrackdelay) ...
      || ( th > ih+1) ...  % track hour is more than an hour later
      || ( th == (ih+1) && ( im < (60-maxtrackdelay) || tm>=maxtrackdelay ) )
      % last condition - track time is less than an hour later, but more than 
      % maxtrackdelay mins
         imgtimes(1,:) = [];
         if isempty( imgtimes ), break; end
         [ih, im] = matsplit( imgtimes(1,:) );   % get img hour mins secs
      end
   end % end for each calcium img time
end

% Convert a date into the provided format. Can't use matlab's datestr
% because that's expecting one of their date vectors or a serial number.
% Date has format: [yyyy mm dd] (can be a matrix with rows of dates)
% Time (optional): [  HH MM SS] (can be a matrix with rows of dates)
% The format string can have bits other than date & time, but those bits
% can't use any of the date & time indiciators, such as yyyy, MM etc
function str = date2str( datenum, formatstr, timenum )
   % sometimes a single date may be provided for a list of times - repmat date
   if size(datenum,1) < size(timenum,1)
      if size(datenum,1) > 1
         str = sprintf( 'Number of dates (%d) and times (%d) must match, else provide a single date\n', ...
                         size(datenum,1), size(timenum,1) );
         fprintf( str );
         return;
      end
      datenum = repmat( datenum, [size(timenum,1) 1] );
   end
      
   stind = @( sym ) strfind( formatstr(1,:), sym ); % start index of date/time symbol
   ind   = @( sym ) ternaryOp( isempty(stind(sym)), [], stind(sym):(stind(sym)+length(sym)-1) );
   resz  = @(  st ) reshape( st, [length(st)/size(datenum,1) size(datenum,1)] )';
   % initialise string to format, & then replace each of the date/time
   % symbols with the actual date & time - make enough strings for number
   % of dates provided
   str = repmat( formatstr, [size(datenum, 1) 1] ); 
   str( :, ind('yyyy') ) = resz( sprintf( '%4d', datenum(:,1) ) );
   str( :, ind(  'mm') ) = resz( sprintf( '%2d', datenum(:,2) ) );
   str( :, ind(  'dd') ) = resz( sprintf( '%2d', datenum(:,3) ) );

   % if format string has times in it also, populate these too
   if nargin>=3 && ~isempty( timenum )
      str( :, ind('HH') ) = resz( sprintf( '%2d', timenum(:,1) ) );
      str( :, ind('MM') ) = resz( sprintf( '%2d', timenum(:,2) ) );
      str( :, ind('SS') ) = resz( sprintf( '%2d', timenum(:,3) ) );
   end
   
   % where a number was, for example, 1 digit, but string format was for 2, an
   % empty space would have been returned from sprintf --> change to 0s
   str( str==' ') = '0';
end

% get times of recordings from directory names using provided format
function times = getTimesFromDirectories( ddir, format, matchtype )
   times = []; % 
   if ~exist( ddir, 'dir' ) 
      return;
   end
   imgregexp = genRegExpFromFormat( format, matchtype );
   timesdir  = dir( ddir );
   for ii=1:length(timesdir)
      timedir  = timesdir(ii).name;
      if ~isempty( regexp( timedir, imgregexp, 'once' ) )
         time  = extractImgTimeFromFilename( timedir, format );
         times = [times; time]; 
      end
   end
end

function ddate = extractDateFromFilename( dirname, dateformat )
   % we don't know how many Y chars there'll be so look for them
   % specifically, which assumes that there's not other Y characters in the
   % directory name format
   year    = dirname( dateformat == 'y' );
   year    = str2double( year );
   % for day & month we know it'll have 2 chars so look for that
   % specifically, to avoid any problems from other text in the directory
   % name, such as 'Experiment Data Sheet - yyyy.mm.dd' creating a problem
   % with an extra d (using matlab's formatting for datestr)
   monthind= strfind( dateformat, 'mm' );
   month   = dirname( monthind:monthind+1 );
   month   = str2double( month );
   dayind  = strfind( dateformat, 'dd' );
   day     = dirname( dayind:dayind+1 );
   day     = str2double( day );
   
   ddate   = [year month day];
end

function imgtime = extractImgTimeFromFilename( imgtimedir, imgformat )
   if length( imgtimedir ) ~= length( imgformat )
      imgtime = [];
      return;
   end
   stind = @( sym ) strfind( imgformat(1,:), sym ); % start index of date/time symbol
   ind   = @( sym ) ternaryOp( isempty(stind(sym)), [], stind(sym):(stind(sym)+length(sym)-1) );

   hour = str2double( imgtimedir( ind( 'HH' ) ) );
   min  = str2double( imgtimedir( ind( 'MM' ) ) );
   sec  = str2double( imgtimedir( ind( 'SS' ) ) );
   
   imgtime = [hour min sec];
end

% generate a regular expression from date formats
function dateregexp = genRegExpFromFormat( format, matchtype )
   if nargin<2, matchtype = '[0-9a-zA-Z\ ]+'; end
   dateparts   = regexp( format, matchtype );
   datesep     = unique( format( dateparts(2:end)-1 ) );
   % if we seem to have more than 1 separator, we'll have to construct the
   % format manually
   if length( datesep ) > 1
      % we start by expecting chars since most OSs require alpha char 1st 
      dateregexp = matchtype;
      for ii=1:length( format )
         if any( format(ii)==datesep )
            sep = ['[' datesep( format(ii)==datesep ) ']'];
            dateregexp = [dateregexp  sep  matchtype]; 
         end
      end
      
   % only 1 separator found so create a reg exp from it 
   else
      datesep     = ['[' datesep ']']; % require 1 separator btwn each component
      % gen reg exp for date as characters separated by separator
      dateregexp  = repmat( [matchtype datesep], [1 length(dateparts)] );
      dateregexp  = dateregexp(1:end-length(datesep)); % remove final separator
   end
end



