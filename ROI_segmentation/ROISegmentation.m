% [err, roi_fname] = ROISegmentation( session, path_name, user_path, video_names )
% 
% Segment two-photon images into individual neurons
% Inputs:
%  session     - date/time session info of trial
%  path_name   - directory to save files to
%  user_path   - user path (backup save directory if saving to farm fails)
%  video_names - cell array of full path to all calcium image files
% Outputs:
%  err         - true if managed to save ROI masks without error, else false
%  roi_fname   - output filename for masks
function [err, roi_fname, masks] = ROISegmentation( session, path_name, user_path, img_names, logfh )
   err = true; roi_fname = []; masks = [];
   maxframes   = 500; % max number of frames to use from each tif video in the segmentation
   allstack_g  = [];
   allmean_r   = [];
   mc_found    = false;  
   
   for i = 1:length( img_names )
      mc_fname = img_names{i};
      
      if ~exist( mc_fname, 'file' )
         str = sprintf( '\tMotion corrected file not found, ignoring (%s)\n', mc_fname );
         cprintf('Errors', str);
         if logfh>0, fprintf( logfh, str ); end
         
      else
         str = sprintf( '\tMotion corrected file found (%s)\n', mc_fname );
         cprintf('Keywords', str);
         if logfh>0, fprintf( logfh, str ); end
         
         % try loading & concatenating green stack, & mean red images from
         % each video (i.e. from each motion corrected tif file)
         try
            load( mc_fname, 'stack_g', 'mean_r' );
            mc_found = true;
            nstack   = min( maxframes, size(stack_g, 3) );
            allstack_g  = cat( 3, allstack_g, stack_g(:,:,1:nstack) );
            allmean_r   = cat( 3, allmean_r,  mean_r );
            
         catch ME
            str = sprintf( '\tError opening motion correction file (%s)\n', mc_fname );
            str = getCatchMEstring( ME, str, true );
            cprintf('Errors', str);
            if logfh>0, fprintf( logfh, str ); end
         end
      end
   end
   if ~mc_found
      str = sprintf( '\tNo motion corrected files found in %s\n', path_name );
      if logfh>0, fprintf( logfh, str ); end
      cprintf('Errors*', str);
      return;
   end      

   clearvars stack_g mean_r;

   % Perform segmentation on these stacks to get the masks
   cellRadius  = 10;
   maxCells    = 200;
   try
      allmean_r = mean( allmean_r, 3 ); % expecting a single red frame so get mean across videos
      [ ~, ~, masks, fig] = neuroSEE_segment(allstack_g, allmean_r, cellRadius, maxCells );
      
      if isempty( masks )
         error( 'NeuroSEE:ROISegmentation:MaskError: empty mask returned');
      end
      if any( isnan( masks(:) ) ) || any( isinf( masks(:) ) )
         error( 'NeuroSEE:ROISegmentation:MaskError: invalid mask with NaN or Inf returned');
      end
      
   catch ME
      str = getCatchMEstring( ME, '\tError running ROISegmentation', true );
      cprintf( 'Errors', str );
      if logfh>0, fprintf( logfh, str ); end
      return;
   end

   roi_fname = fullfile(path_name, [session, '_ROImask.mat']);
   roi_fig   = fullfile(path_name, [session, '_ROImask']);
   
   % try saving to server, if that fails save to user's local directory
   try
      save( roi_fname, 'masks', '-v7.3' );
      str = sprintf( '\tROI masks saved in %s\n', roi_fname );
      cprintf( 'Keywords', str );
      if logfh>0, fprintf( logfh, str ); end
      
   catch 
      % ran out of space so it failed - try saving to my local dir
      if ~exist( user_path, 'dir' )
         mkdir( user_path );
      end
      roi_fname = fullfile(user_path, [session, '_ROImask.mat']);
      roi_fig   = fullfile(user_path, [session, '_ROImask']);
      
      try
         save( roi_fname, 'masks', '-v7.3' );
         str = sprintf( 'Failed saving roi on farm, saving locally (%s)\n', roi_fname );
         cprintf( 'Errors', str );
         if logfh>0, fprintf( logfh, str ); end
         
      catch ME
         str = getCatchMEstring( ME, 'Error saving ROI segmentation masks', true );
         cprintf('Errors*', str);
         if logfh>0, fprintf( logfh, str ); end
      end
   end
   
   % if we've saved the masks consider segmentation successful
   err = false;
   
   % try saving figure if 'fig' exists & is not empty
   if ~exist('fig','var') || isempty(fig)
      str = 'Invalid fig - cannot save ROI segmentation figure';
      cprintf( 'Errors', str );      
      return;
   end
   
   try
      saveFigure( fig, roi_fig, true, {'fig','pdf'} );
      
   catch ME
      str = getCatchMEstring( ME, 'Error saving ROI segmentation figure', true );
      if logfh>0, fprintf( logfh, str ); end
      cprintf('Errors*', str);
      
      return;
   end
   
end