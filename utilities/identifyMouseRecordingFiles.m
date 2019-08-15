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
      
   % months must not be capital to separate from minutes
   imgformat   = 'yyyymmdd_HH_MM_SS_2P';
   trackformat = 'Track_yyyy-mm-dd-HH-MM-SS';
   photondir   = '2P';
   trackdirname= 'Neurotar';
   matchtype   = '[0-9a-zA-Z ]+'; % match dir & filenames based on alphanumeric chars (but not '_')

   mouseDataList.name     = mousename;
   mouseDataList.num      = mousename(2:end);

  % Now go through the log book files for file matching mousename
   lognames  = dir( logbookdir );
   for li=1:length( lognames )
      % we're looking at filenames here
      filename = lognames(li).name;
      % if filename is name of a mouse, proceed
      if length(filename) > 2
          if strcmp(filename(1:length(mousename)),mousename)
             [~,sheets] = xlsfinfo( fullfile( logbookdir, filename) );
             % cycle through each worksheet name & extract ones with fam1, fam1fam2,
             % fam1nov, fam1fam1rev, fam2, open
             numexp = 0;
             for si=1:length(sheets)
                if strcmp(sheets{si}(1:3),'fam') || strcmp(sheets{si}(1:3),'ope')
                   numexp = numexp + 1; 
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
                   [y,m,d,H,M,S] =  datevec(datetime( num, 'ConvertFrom', 'excel' )); 
                   invalid = cellfun( @(c) any(isnan(c)), raw );
                   valid   = cellfun( @(c) length(c)>1 && c(1)=='X', raw );
                   valid   = find( ~invalid & valid );
                   y = y(valid); m = m(valid); d = d(valid);
                   H = H(valid); M = M(valid); S = S(valid);

                   % now get environment type for each valid row
                   [~,~,raw] = xlsread(fullfile( logbookdir, filename), sheets{si}, 'J16:J40');
                   if mergedcells
                      raw = raw(1:end,10);
                   end
                   env = raw( valid );

                   % add experiment for current date, including recording times 
                   mouseDataList.exp{numexp}.exp      = exp;
                   mouseDataList.exp{numexp}.date     = expdate;
                   mouseDataList.exp{numexp}.logtimes = [y m d H M S];
                   mouseDataList.exp{numexp}.env      = env;
                end
             end
          end % end for matching logbook name
      end % end for if length(filename)>2
   end % end for logbook names

   % We need to loop through each experiment, and then look at each time
   % and find the corresponding Ca imaging and Neurotar tracking files
   recorddates = dir( gcampdir );
   for ei = 1:numel( mouseDataList.exp )
        logtimes = mouseDataList.exp{ei}.logtimes;
        for ri = 1:length( recorddates )
            if strcmp(mouseDataList.exp{ei}.date, recorddates(ri).name)
                % 2 photon directory - will have a list of directory names from times
                % of image recording, e.g. 20181016_10_34_37_2P
                imgdir      = fullfile( gcampdir, recorddates(ri).name, photondir );
                imgtimes    = getTimesFromDirectories( imgdir, imgformat, matchtype );
                % tracking directory - will have a list of directory names from times
                % of track recording, e.g. Track_2018-10-16-10-16-01
                trackdir    = fullfile( gcampdir, recorddates(ri).name, trackdirname );
                tracktimes  = getTimesFromDirectories( trackdir, trackformat, matchtype );
                if ~isempty( imgtimes ) && ~isempty( tracktimes )
                    % get times for this date where we have valid calcium imaging & track data
                    [finalimgtimes, finaltracktimes] = matchImgAndTrackTimes(logtimes, imgtimes, tracktimes);
                    if ~isempty( imgtimes )
                       mouseDataList.exp{ei}.imgtimes   = finalimgtimes(:,4:6);
                       mouseDataList.exp{ei}.tracktimes = finaltracktimes(:,4:6);
                    end
                end
            end
        end % end for each recording date folder
   end % end for each experiment
end

function [finalimgtimes, finaltracktimes] = matchImgAndTrackTimes(logtimes, imgtimes, tracktimes)
    maxtrackdelay = 12; % track recording must begin within this many mins of Ca image
    fi = 0;
    for i = 1:length(logtimes)
        logtime = logtimes(i,:);
        logtime(1:3) = imgtimes(1,1:3);
        dt1 = etime(imgtimes,logtime);
        imgtime = imgtimes( abs(dt1)==min(abs(dt1)), : );
        if ~isempty(imgtime) && size(imgtime,1)==1
            timediff = etime(tracktimes,imgtime);
            valid_ind =  timediff>0 & timediff<=maxtrackdelay*60;
            if ~isempty(valid_ind)
                fi = fi + 1;
                finaltracktimes(fi,:) = tracktimes(timediff==min(timediff(valid_ind)),:);
                finalimgtimes(fi,:) = imgtime;
            end
        end
    end
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

function imgtime = extractImgTimeFromFilename( imgtimedir, imgformat )
   if length( imgtimedir ) ~= length( imgformat )
      imgtime = [];
      return;
   end
   stind = @( sym ) strfind( imgformat(1,:), sym ); % start index of date/time symbol
   ind   = @( sym ) ternaryOp( isempty(stind(sym)), [], stind(sym):(stind(sym)+length(sym)-1) );

   year = str2double( imgtimedir( ind( 'yyyy' ) ) );
   month= str2double( imgtimedir( ind( 'mm' ) ) );
   day  = str2double( imgtimedir( ind( 'dd' ) ) );
   hour = str2double( imgtimedir( ind( 'HH' ) ) );
   min  = str2double( imgtimedir( ind( 'MM' ) ) );
   sec  = str2double( imgtimedir( ind( 'SS' ) ) );
   
   imgtime = [year month day hour min sec];
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

