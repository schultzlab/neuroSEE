% UNFINISHED
% Written by Ann Go (adapted from Katie's processNeuroSEEmouseFiles.m)
%
% This generates .mat file (mouseDatasets.mat) containing structure of mouse datasets -
% mDatasets.
% Each cell of mDatasets is a structure with the following fields
%   name : 'm62'
%   exp  : {[1x1 struct] [1x1 struct]}

mice     = {'m62'}; % {'m66', 'm62', 'm69', 'm70', 'all'} 
                    % or {'20190310', '20190311'}
out_dir  = fullfile( '/Volumes','RDS','project','thefarm2','live','CrazyEights','AD_2PCa','Processed');
out_file = fullfile( out_dir, 'm62_pipeline.mat' );

% if we already have mouse dataset list, load it & merge it with new list
if exist( out_file, 'file' )
   load( out_file, 'mouseDset' );
   mousenames = getStructFieldFromCell( mouseDset, 'name' );
else
   mousenames = [];
end

% extract recorded mouse experiment list from log book & Gcamp directory
if ~exist('micerec','var')   
   % Get mouse info
   try
      miceexp = identifyMouseExperimentDates;
      miceexp_names = getStructFieldFromCell( miceexp, 'name' );
   catch ME
      str = sprintf('Error getting mouse experiment details (%s, line %d, function %s)\n',...
                     ME.message, ME.stack(1).line, ME.stack(1).name);
      return;
   end
end

