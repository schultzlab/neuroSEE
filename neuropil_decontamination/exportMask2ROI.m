% err = exportImageJABC( fname, roi )
% OR
% [err, fname, coords] = exportImageJABC( fname, roi )
%
% Write a binary file with the roi polygon info for ImageJ

% Note: if on a mac, DON'T use the compress function in finder to create
% the ROI zip file, as the read imagej function in python won't process it.
% Heaven only knows why! But if you use the zip function in the terminal
% it's fine. Bah!
function varargout = exportMask2ROI( ddir, fname, mask )
   err     = true; 
   
   if ~any( strcmpi( fname(end-2:end), {'bin','roi'} ) )
      fname = [fname '.bin'];
   end
   file = fullfile( ddir, fname );
   
   % if file exists make sure user wants to write over it
%    if exist( file, 'file' )
%       str = sprintf( 'File %s exists, write over it [y/n]? ', fname ); 
%       in = input( str, 's' );
%       if ~strcmpi( in, {'y', 'yes', 'ys'} )
%          return;
%       end
%    end
   fh = fopen( file, 'w' );
   
   if fh <= 0
      cprintf( 'Keywords*', 'Invalid file handle, exiting...');
      outputs = setOutputs( nargout, err, file );
      return;
   end
   
   [x, y] = getPolygon( mask );
   sz     = size( mask );
   rect   = [0 0 sz];
   rect   = [min(y) min(x) max(y) max(x)];
   err    = writeROI( fh, rect, x, y, sz, fname );
   
   outputs = setOutputs( nargout, err, file, x, y );
   varargout = outputs;
end

function outputs = setOutputs( nargout, err, fname, x, y )
   outputs = cell( nargout, 1 );
   switch nargout
      case 0
         outputs    = cell(0);
      case 1
         outputs{1} = err;
      case 2
         outputs{1} = err;
         outputs{2} = fname;
      case 3
         outputs{1} = err;
         outputs{2} = fname;
         outputs{3} = [x(:) y(:)];
      otherwise
         outputs{1} = err;
         outputs{2} = fname;
   end   
end

function [x, y] = getPolygon( mask )
   poly = mask2poly( mask );
   if length( poly ) > 1
      l = arrayfun( @(i) poly(i).Length, 1:length( poly ) );
      [~, i] = max( l );
      poly = poly(i);
   end
   x = poly.X;
   y = poly.Y;
   
%    boundaries = bwboundaries( mask );
%    % if more than one polygon found, assume the largest is the cell
%    if length( boundaries ) > 1
%       l = cellfun( @length, boundaries );
%       [~,i] = max( l );
%       poly = boundaries{ i };
%    else
%       % only 1 polygon so no brainer
%       poly = boundaries{1};
%    end
%    x    = poly(:,2);
%    y    = poly(:,1);
end

