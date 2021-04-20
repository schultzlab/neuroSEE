% Written by Ann Go (adapted from Katie's identifyMouseExperimentDates.m)
% 
% This function finds the tracking file in dir_track that matches the 2P image 
% INPUTS
%   data_locn   : directory where GCaMP6 data is
%   file        : part of file name of 2P image in the format
%                   yyyymmdd_HH_MM_SS
%   force       : ignores existing mat file for matching tracking file
% OUTPUT
%   fname_track : filename of tracking file matching 2P image

function fname_track = findMatchingTrackingFile(data_locn, file, force)
    if nargin<3, force = 0; end
        
    dir_track = [data_locn 'Data/' file(1:8) '/Neurotar/'];
    dir_processed = [data_locn 'Data/' file(1:8) '/Processed/' file '/behaviour/'];
        if ~exist(dir_processed,'dir'), mkdir(dir_processed); end
    
    imgformat   = 'yyyymmdd_HH_MM_SS';
    imgtime     = datevec(file,imgformat);
    matchtype   = '[0-9a-zA-Z ]+'; % match dir & filenames based on alphanumeric chars (but not '_')
    [tracktimes, trackformat, Sind]  = getTrackTimesFromDirectory( dir_track, matchtype );
    finaltracktime = matchImgAndTrackTimes( imgtime, tracktimes );
    trackfolder = [dir_track trackformat(1:Sind-1) datestr( finaltracktime, trackformat(Sind:end) )];
    fname_track = findTrackFileName( trackfolder, dir_processed, force );

    function [times,format,Sind] = getTrackTimesFromDirectory( ddir, matchtype )
       times = [];  
       if ~exist( ddir, 'dir' ) 
          return;
       end
       timesdir  = dir( ddir );
       for i=1:length(timesdir)
          timedir  = timesdir(i).name;
          if contains(timedir,'Saved')
            format = 'SavedTrack_yyyy-mm-dd-HH-MM-SS'; Sind = 12;
          else
              format = 'Track_yyyy-mm-dd-HH-MM-SS'; Sind = 7;
          end
          trackregexp = genRegExpFromFormat( format, matchtype );
          if ~isempty( regexp( timedir, trackregexp, 'once' ) )
             time  = extractTimeFromFilename( timedir, format );
             times = [times; time]; 
          end
       end
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

    function filetime = extractTimeFromFilename( timedir, format )
       if length( timedir ) ~= length( format )
          filetime = [];
          return;
       end
       stind = @( sym ) strfind( format(1,:), sym ); % start index of date/time symbol
       ind   = @( sym ) ternaryOp( isempty(stind(sym)), [], stind(sym):(stind(sym)+length(sym)-1) );

       year = str2double( timedir( ind( 'yyyy' ) ) );
       month= str2double( timedir( ind( 'mm' ) ) );
       day  = str2double( timedir( ind( 'dd' ) ) );
       hour = str2double( timedir( ind( 'HH' ) ) );
       min  = str2double( timedir( ind( 'MM' ) ) );
       sec  = str2double( timedir( ind( 'SS' ) ) );

       filetime = [year month day hour min sec];
    end

    function finaltracktime = matchImgAndTrackTimes(imgtime, tracktimes)
       maxtrackdelay = 12; % track recording must begin within this many mins of Ca image

       % now go through tracking files, & find the files that
       % are within a few minutes of imaging time
       timediff = etime(tracktimes,imgtime);
       valid_ind =  timediff>0 & timediff<=maxtrackdelay*60;
       finaltracktime = tracktimes(timediff==min(timediff(valid_ind)),:);
    end

    function fname_track = findTrackFileName(trackfolder,dir_processed,force)
        if nargin<2, force = 0; end
        if force
            csvfile = subdir(fullfile(trackfolder,['*.','csv'])); % Files from 2018
            if numel(csvfile) == 1
                fname_track = csvfile.name;
            else 
                tdmsfile = subdir(fullfile(trackfolder,['*.','tdms'])); % Files from 2019
                if numel(tdmsfile) == 1
                    fname_track = tdmsfile.name;
                else
                    fprintf('%s: Tracking file does not exist in folder\n', file);
                end
            end
        else
            matfile = subdir(fullfile(dir_processed,['*.','mat']));
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
                    end
                end
            else
                fname_track = matfile.name;
            end            
        end
    end

end % function