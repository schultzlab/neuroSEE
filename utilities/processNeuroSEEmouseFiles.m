%% process neuroSEE mouse files

%% Note
% - allow saving a tif file after motion correction so can write frame by frame
% - assumes that each tif file starts with an odd numbered frame, so that
%   it's a green frame (the red frames are even numbered)
% - assume that experiments are in date order, so can grab new experiments
%   from the end of the list 
% - changed Seig's spike extraction because it hard coded AR(1) coeff to 0.9
% - re-wrote tif reading function to be much faster
% - changed mode shift to include a histogram coz mode of a real number is nothing
%

%% TO DO:
% - check that number of motion corrected files = number of images file &,
%   if not, try motion correcting the missing files again
% - remove data from running stages when saving mouse output file, else the
%   data's saved in individual stage files, as well as pipeline file

% Unfinished directories
% - list directories that don't have all necessary tif files

server_dir = 

init     = false;
% mice     = 'all';  % cell array or string of which mice to gen tif files for (or all)
mice     = {'m62'}; % {'m66', 'm62', 'm69', 'm70'}; % 'all' % 
maxerr   = 30;
if contains( computer, 'win', 'ignorecase', true )
   out_dir = fullfile( 'W:','live','CrazyEights','AD_2PCa','katie');
else
   out_dir = fullfile( '/Volumes','Schultz_group_data','Crazy Eights','Ann','katie');
end
out_file = fullfile( out_dir, 'mouse_pipelines.mat' );

% if we already have mouse pipelines, load them & merge them with mice
% recordings (i.e. experimental data) to make sure we use previous processing 
% info + new experiments
if exist( out_file, 'file' )
   load( out_file, 'mouseexp' );
   mousenames = getStructFieldFromCell( mouseexp, 'name' );
else
   mousenames = [];
end

% extract recorded mouse experiment list from log book & Gcamp directory
if init || ~exist('micerec','var')   
   % Get mouse info
   try
      micerec = identifyMouseRecordingDates;
      micerec_names = getStructFieldFromCell( micerec, 'name' );
   catch ME
      str = sprintf('Error getting mouse experiment details (%s, line %d, function %s)\n',...
                     ME.message, ME.stack(1).line, ME.stack(1).name);
      return;
   end
end

% now we gotta see which mice we already know from previous processing, &
% copy over any new mouse info
if isempty( mousenames )
   % if we have no mouse names yet, all mice are knew, so just rename mice
   % recording struct to mouse experiment struct
   mouseexp = micerec;
else
   % we already know some mice, so check for new ones
   for mi=1:length( micerec_names )
      know = strcmpi( micerec_names{mi}, mousenames );
      if ~any( know )
         % if don't know this mouse yet, add to cell array of mice
         mouseexp(end+1) = micerec( mi );
      else
         % if we already know this mouse, check for new experiments
         if length( mouseexp{mi}.exp ) ~= length( micerec{know}.exp )
            % new experiments, so copy them over - easiest to copy old ones
            % across to mouse recordings, & then copy all of them back to mouseexp
            nold = length( mouseexp{mi}.exp );
            micerec{know}.exp(1:nold) = mouseexp{mi}.exp;
            mouseexp{mi}.exp = micerec{know}.exp; 
         end
      end
   end
end


% if user's requested to process all mice, get a list of their names
if strcmpi( mice, 'all' )
   mice = getStructFieldFromCell( mouseexp( randperm( length( mouseexp ) ) ), 'name' );
   % just so we process different mice each time, randomise the order
end

% For each mouse that user's interested in, process their valid experiments
totalerrors = 0;
for mi=1:length( mouseexp )
   mouse = mouseexp{ mi };
   
   % if user's interested in this mouse, process it
   if any( strcmpi( mouse.name, mice ) )
      % For each experiment extract the valid experiment times
      for exp_ind=1:length( mouse.exp )
         if isfield(mouse.exp{exp_ind}, 'ddir')
            % Process each directory with valid calcium & tracking files
            cadir     = mouse.exp{exp_ind}.cadirs;
            numtrials = size( cadir, 1 );

            for trial_ind=1:numtrials 
               % get mouse pipeline if it exists (can't create empty array
               % cuz then it's assumed to be type double, instead of pipeline)
               if ( isfield( mouse.exp{exp_ind}, 'pipelines' ) ...
                 && length( mouse.exp{exp_ind}.pipelines ) >= trial_ind )
              
                  mousepl = mouse.exp{exp_ind}.pipelines( trial_ind ); 
                  % initialise pipeline if not already done
                  if ~mousepl.IsPipelineInitialised
                     mousepl = mousepl.InitialisePipeline( mouse, exp_ind, trial_ind, out_dir );
                  end
               else
                  % no mouse pipeline for this trial - create
                  mousepl = MousePipeline( mouse, exp_ind, trial_ind );
               end
               
               err = processTrial( mousepl ); 
               mouse.exp{exp_ind}.pipelines( trial_ind ) = mousepl;
               mouseexp{ mi } = mouse; 
               save( out_file, 'mouseexp' );

               % if we need to do any of the steps in the pipeline, process
               % the trial~
               totalerrors = totalerrors + double( err );
               if totalerrors > maxerr
                  disp( 'Too many errors, exiting...' );
                  return;
               end

            end % end for each Calcium img directory
         end
      end % end if user requested to process this mouse
   end % end for each mouse experiment
end % end for each mouse

% % Write the directory of where you want to save the processed files
% path_name = [pwd, '/m62/fam1_nov/fam1/'];
% 
% % Write the details of the session and the mouse
% sessiondate = '151018'; %DDMMYY format
% 
