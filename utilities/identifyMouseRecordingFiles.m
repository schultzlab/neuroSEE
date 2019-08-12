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
function [mouseexp, recordings] = identifyMouseRecordingFiles( data_locn )

   logbookdir  = fullfile( data_locn, 'Digital Logbook', 'by animal' );
   gcampdir    = fullfile( data_locn, 'Data' );
   dateformat  = 'yyyymmdd'; % using matlab's formatting for datestr
   matchtype   = '[0-9a-zA-Z ]+'; % match dir & filenames based on alphanumeric chars (but not '_')
      
   % months must be capital to separate from minutes
   imgformat   = 'yyyymmdd_HH_MM_SS_2P';
   trackformat = 'Track_yyyy-mm-dd-HH-MM-SS';
   photondir   = '2P';
   trackdirname= 'Neurotar';

   % We need to loop through each logbook date, and then look at each time to see what
   % mice are in the file - identify them by the name of the worksheets in the files
   % 2018.10.16/Neurotar/Track_2018-10-16-08-59-05/SavedTrack-2018-10-16-08-59-05.csv
   % 2018.10.16/2P/20181016_09_14_03_2P/20181016_09_14_03_2P_XYT.raw

   recorddates = dir( gcampdir );
   numdates    = length( recorddates );

   % generate regexp for date without knowing in advance what the format is -
   % - identify date separator, e.g. 10.09 or 10/09 etc
   % in farm2 there's no separator, so create manually
   dateregexp = '^\d{4}\d{2}\d{2}$'; % ^ is start of line, $ is end
   logregexp  = 'ExptNotes_\d{4}\d{2}\d{2}.xlsm'; % starts with ExptNotes, ends with .xlsm
   
   % After cycling through each date, we'll have a cell array containing
   % the dates that valid data was generated, and the image & tracking
   % times for each
   numrec = 0;
   for ri=1:numdates
      % make sure directory is a valid date directory  
      dirname = recorddates(ri).name;
      % if we find the date regular expression in the dir name, then proceed
      if ~isempty( regexp( dirname, dateregexp, 'once' ) ) 
         % extract date that data was recorded
         dirdate     = extractDateFromFilename( dirname, dateformat );
         % 2 photon directory - will have a list of directory names from times
         % of image recording, e.g. 20181016_10_34_37_2P
         imgdir      = fullfile( gcampdir, dirname, photondir );
         imgtimes    = getTimesFromDirectories( imgdir, imgformat, matchtype );
         % tracking directory - will have a list of directory names from times
         % of track recording, e.g. Track_2018-10-16-10-16-01
         trackdir    = fullfile( gcampdir, dirname, trackdirname );
         tracktimes  = getTimesFromDirectories( trackdir, trackformat, matchtype );
         if ~isempty( imgtimes ) && ~isempty( tracktimes )
            % get times for this date where we have valid calcium imaging & track data
            [finalimgtimes, finaltracktimes, ok] = matchImgAndTrackTimes(imgtimes, tracktimes);
            if ~ok
               return;
            end
            if ~isempty( finalimgtimes )
               numrec = numrec + 1;
               recordings{numrec}.date      = dirdate;
               recordings{numrec}.Catimes   = finalimgtimes;
               recordings{numrec}.tracktimes= finaltracktimes;
            end
         end         
      end % end if valid gcamp directory
   end % end for each gcamp directory
   
   % Now go through the log book files & see what mice were active for each
   % of the valid data recordings
   logn  = dir( logbookdir );
   currnames = []; nummice = 0;
   for li=1:length( logdates )
      % we're looking at filenames here
      filename = logdates(li).name;
     [status,sheets] = xlsfinfo( fullfile( logbookdir, filename) );
     % cycle through each worksheet name & extract ones with fam1fam1,
     % fam1nov, fam1fam1rev, open
     for si=1:length(sheets)
        if sheets{si}(1)=='m' && length(sheets{si})==3
           mousenum = str2double( sheets{si}(2:end) );
           % go through recording times & extract the recording start times
           % - data starts on row 15 of column A, but ranges have to
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
           valid   = cellfun( @(c) length(c)==2 && c(1)=='X', raw );
           valid   = find( ~invalid & valid );
           h = h(valid); m = m(valid); s = s(valid);

           % now get experiment type for each valid row
           if farm==1
              [~,~,raw] = xlsread(fullfile( logbookdir, filename), sheets{si}, 'J15:J40');
           elseif farm==2
              [~,~,raw] = xlsread(fullfile( logbookdir, filename), sheets{si}, 'J16:J40');
           end
           if mergedcells
              raw = raw(1:end,10);
           end
           exp = raw( valid );

           % add mouse data to struct - add new mouse if not seen before
           name = sprintf('m%2d', mousenum); name(name==' ') = 0;

           % if mouse not already known, add to struct
           if ~any( strcmpi( name, currnames ) )
              nummice   = nummice + 1; 
              currmouse = nummice; 
              mouseexp{nummice}.name     = name;
              mouseexp{nummice}.num      = mousenum;
              mouseexp{nummice}.exp      = cell(0);
              mouseexp{nummice}.numexp   = 0;
              currnames{nummice}         = name;
           else
              currmouse = find( strcmpi( name, currnames ) );
           end
           % add experiment for current date, including recording times 
           numexp = mouseexp{currmouse}.numexp + 1;
           mouseexp{currmouse}.exp{numexp}.exp   = lower( title2Str( cell2mat( unique( exp )' ) ) );
           mouseexp{currmouse}.exp{numexp}.date  = extractDateFromFilename( filename, logformat );
           mouseexp{currmouse}.exp{numexp}.times = [h m s];
           mouseexp{currmouse}.exp{numexp}.type  = exp;
           mouseexp{currmouse}.numexp = numexp;
        end
     end
   end % end for logbook dates
   
   % now you need to compare logbook recordings with mouse data - where
   % they're a match, we're good to go!
   
   % for each mouse, & each experiment for the mouse, find the imaging &
   % tracking files
   recorddates = cell2mat( getStructFieldFromCell(recordings, 'date','cell')' );
   mi = 1; % use while loop cuz removing mice with no valid exp so num mice changes
   while mi<=length(mouseexp)
      mouse = mouseexp{mi};
      mname = mouse.name;
      mnum  = mouse.num;
      numexp= length( mouse.exp );
      ii = 1; % use while loop coz we delete invalid experiments
      while ii<=length(mouse.exp)
         exp     = mouse.exp{ii};
         expdate = exp.date;
         exptimes= exp.times;
         exptype = exp.type;
         
         recind  = find( ismember( recorddates, expdate, 'rows' ) );

         % if there are no valid img/track datasets for this experiment, ditch it
         if isempty( recind ) || isempty( exptimes )
            mouse.exp(ii) = [];
            
         % if there are valid data files for this logged recording, add it
         % to the mouse's video names & track names
         else
            catimes = recordings{recind}.Catimes;
            tracktimes = recordings{recind}.tracktimes;
            
            % some of the recordings will have been aborted or invalid, so
            % check times with the experiment times for the mouse
            %  - hour & min of exp matches track times
            if length( catimes ) >= length( exptimes )
               [keepca, keeptime] = ismember( tracktimes(:,1:2), exptimes(:,1:2), 'rows' );
               catimes    = catimes(keepca,:);
               tracktimes = tracktimes(keepca,:);
               exptimes   = exptimes( keeptime(keeptime>0), : );
            elseif length( catimes ) <= length( exptimes )
               [keeptrack, keeptime] = ismember( exptimes(:,1:2), tracktimes(:,1:2), 'rows' );
               exptimes   = exptimes(keeptrack,:);
               exptimes   = exptimes( keeptime(keeptime>0), : );
            end
            
            % at this point we have all calcium and track recording times
            % for all valid experiments for this mouse - get filenames
            if ~isempty( exptimes )
               imgsubdirs    = date2str( expdate, imgformat, catimes );
               imgdir        = fullfile( gcampdir, recordings{recind}.dirname, photondir );
               imgsubdirs    = [repmat([imgdir filesep], [size(imgsubdirs,1) 1])    imgsubdirs];
               tracksubdirs  = date2str( expdate, trackformat, tracktimes );
               imgdir        = fullfile( gcampdir, recordings{recind}.dirname, trackdirname );
               tracksubdirs  = [repmat([imgdir filesep], [size(tracksubdirs,1) 1])  tracksubdirs];

               exp.times     = exptimes;
               exp.ddir      = fullfile( gcampdir, recordings{recind}.dirname );
               exp.cadirs    = imgsubdirs;
               exp.trackdirs = tracksubdirs;
               mouse.exp{ii} = exp;
               ii = ii+1; % increment to get next experiment
               
            % no valid track & imaging times so ditch exp
            else
               mouse.exp(ii) = []; % shift remaining exp down so don't increment i
            end
            
         end       
      end % end for each experiment for this mouse
      
      % before we save this mouse make sure that there are valid experiments
      if isfield( mouse, 'exp' ) && length( mouse.exp ) > 0
         mouseexp{mi} = mouse;
         mi = mi + 1; % increment to get next mouse
      else
         mouseexp(mi) = []; % don't increment cuz we've shifted mice down
      end
   end % end for each mouse
end

function [finalimgtimes, finaltracktimes, ok] = matchImgAndTrackTimes(imgtimes, tracktimes)
   maxtrackdelay = 3; % track recording must begin within this many mins of Ca
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



