% UNFINISHED
% Written by Ann Go (adapted from Katie's processNeuroSEEmouseFiles.m)
%
% This generates .mat file (mDataList.mat) containing structure of mouse datasets -
% mDataList.
% Each cell of mDataList is a structure with the following fields
%   name : 'm62'
%   exp  : with fields exp, data, Cafiles, trackfiles, env

% mice     = {'m66', 'm62', 'm69', 'm70', 'all'}; 

function miceDataList = generate_DataList( mice )

[data_locn,~,err] = load_neuroSEEmodules(false);
if ~isempty(err)
    beep
    cprintf('Errors',err);    
    return
end

out_dir  = [data_locn 'Digital_Logbook/lists/'];
out_file = [out_dir 'mDataList.mat'];

% if we already have master mouse dataset list, load it & merge it with new list
if exist( out_file, 'file' )
   load( out_file, 'mDataList' );
   mDataList_mousenames = getStructFieldFromCell( mDataList, 'name' );
   mDataList_micenum = numel( mDataList_mousenames );
else
   mDataList_mousenames = [];
   mDataList_micenum = 0;
end

% if mouse is not on the list, extract recorded mouse experiment list from
% log book, Gcamp directory and Neurotar directories
for i = 1:numel(mice)
    mouse = mice{i};
    if ~any( strcmp( mDataList_mousenames, mouse ) )
       % Get mouse info
       try
          miceDataList{i} = identifyMouseRecordingFiles( data_locn, mouse );
          mDataList{ mDataList_micenum + 1 } = miceDataList{i};
          mDataList_micenum = mDataList_micenum + 1;
       catch ME
          str = sprintf('Error getting mouse experiment details (%s, line %d, function %s)\n',...
                         ME.message, ME.stack(1).line, ME.stack(1).name);
          return;
       end  
    else
        ind = find( strcmp(mDataList_mousenames, mice) );
        miceDataList{i} = mDataList{ ind };
    end
end

