% writeROI( fh, rect, x, y )
%
% Writes a binary ImageJ ROI file of type polygon from a list of coordinates.
%
% Inputs:
%  fh    - valid file handle ready for writing, or fuyll filename with path
%  rect  - rectangular bounds of the ROI
%  x     - x coordinates of the polygon
%  y     - y coordinates of the polygon
%  sz    - can save size in roi header
%  fname - can save filename string in roi header
function err = writeROI( fh, rect, x, y, sz, fname )
   if nargin<5 || isempty(sz), sz = 0; end
   err = true;
   defext    = '.roi'; % default file extension
   
   machine   = 'ieee-be';
   prec      = 'int16'; % precision is signed shorts
   vers      = 227; % I've no idea so making it the oldest 
   roitype   = 0;
   top       = rect(1); 
   left      = rect(2); 
   bottom    = rect(3);
   right     = rect(4);
   ncoords   = length( x );
   skip      = 0;
   h2offset  = 64 + length( x ) * 2 * 2; % header + 2 bytes for each x & y coord
   offset    = 3; % to make my file match fissa roi exactly, need an offset of this many 32 long ints
   nameoffset= h2offset + 52 + offset*4; % fissa roi files had space for 3 extra long ints
   
   % you're supposed to pass in a handle to the file, but can optionally
   % input the filename string in fname so it can be written to the header
   % - but if user's input an empty file handle or a string instead of a
   % file handle, then open a file for writing & changet fname from full
   % path to just the name string for saving in the header
   if nargin<6 || isempty( fname )
      % if no name provided & fh is a string rather than a file handle, get
      % name string for the header from fh
      if ischar( fh )
         [~,fname,~] = fileparts( fh );
      % if no name string, & also fh is a file handle, use a default name
      else
         fname = 'unnamed';
      end
   end
   
   % if file handle is empty get filename from fname & open a file handle
   % - if fname has full path to it, extract just name for header string
   if isempty( fh )
      [path, fname, ext] = fileparts( fname ); 
      if isempty( ext )
         ext = defext;
      end
      fh = fopen( fullfile( path, [fname ext]), 'w' );
      
   % if file handle is a string array, use it to open a file handle
   elseif ischar( fh )
      % check there's a file extension provided
      [path, fh, ext] = fileparts( fh );
      if isempty( ext )
         ext = defext;
      end
      if ~exist( path, 'file' )
         mkdir( path )
      end
      fh = fopen( fullfile( path, [fh ext] ), 'w' );
      % make sure name string doesn't have path & extension included
      [path, fname, ext] = fileparts( fname ); 
      
   % file handle provided but it's not valid 
   elseif ishandle( fh ) && fh<0 
      % fh is a file handle so if it's valid we can write to it directly
      str = 'Invalid file handle, exiting...';
      cprintf( 'Keywords*', str );
      return
      
   % valid file handle provided
   else      
      % make sure name string doesn't have path & extension included
      [path, fname, ext] = fileparts( fname ); 
      
   end
   namelen   = ceil( length( fname ) / 2 ) * 2;
   
   try
      % fwrite( fh, A, prec, skip, machine );
      fwrite( fh,  "Iout",   'char',   skip, machine ); % 0-3
      fwrite( fh,    vers,     prec,   skip, machine ); % 4-5
      fwrite( fh, roitype,     prec,   skip, machine ); % 6-7
      fwrite( fh,     top,     prec,   skip, machine ); % 8-9
      fwrite( fh,    left,     prec,   skip, machine ); % 10-11
      fwrite( fh,  bottom,     prec,   skip, machine ); % 12-13
      fwrite( fh,   right,     prec,   skip, machine ); % 14-15
      fwrite( fh, ncoords,     prec,   skip, machine ); % 16-17
      fwrite( fh,       0,'float32',   skip, machine ); % 18-21: x1
      fwrite( fh,       0,'float32',   skip, machine ); % 22-25: y1
      fwrite( fh,       0,'float32',   skip, machine ); % 26-29: x2
      fwrite( fh,       0,'float32',   skip, machine ); % 30-33: y2
      fwrite( fh,       0,     prec,   skip, machine ); % 34-35: stroke width
      fwrite( fh,       0, 'uint32',   skip, machine ); % 36-39: ShapeRoi size
      fwrite( fh,       0, 'uint32',   skip, machine ); % 40-43: stroke colour
      fwrite( fh,       0, 'uint32',   skip, machine ); % 44-47: fill colour
      fwrite( fh,       0,     prec,   skip, machine ); % 48-49: subtype (not supported yet, so 0)
      fwrite( fh,       0,     prec,   skip, machine ); % 50-51: options
      fwrite( fh,       0,  'uint8',   skip, machine ); % 52-52: arrow style
      fwrite( fh,       0,  'uint8',   skip, machine ); % 53-53: arrow head style
      fwrite( fh,       0,     prec,   skip, machine ); % 54-55: rounded rect arc size
      fwrite( fh,      41, 'uint32',   skip, machine ); % 56-59: position
      fwrite( fh,h2offset, 'uint32',   skip, machine ); % 60-63: header 2 offset

      % 64 onwards
      for i=1:length(x)
         fwrite( fh, x(i)-rect(2), prec,   skip, machine ); 
      end
      for i=1:length(y) 
         fwrite( fh, y(i)-rect(1), prec,   skip, machine );       
      end
      
      fwrite( fh,      0, 'uint32',   skip, machine ); %  0- 3: blank
      fwrite( fh,      0, 'uint32',   skip, machine ); %  4- 7: C position
      fwrite( fh,      0, 'uint32',   skip, machine ); %  8-11: Z position
      fwrite( fh,      0, 'uint32',   skip, machine ); % 12-15: T position
      fwrite( fh,nameoffset, 'uint32',skip, machine ); % 16-19: name offset
      fwrite( fh,namelen, 'uint32',   skip, machine ); % 20-23: name length
      fwrite( fh,      0, 'uint32',   skip, machine ); % 24-27: overlay label colour
      fwrite( fh,      0,     prec,   skip, machine ); % 28-29: overlay font size
      fwrite( fh,      0,  'uint8',   skip, machine ); % 30-30: empty byte
      fwrite( fh,      0,  'uint8',   skip, machine ); % 31-31: image opacity
      fwrite( fh,      0, 'uint32',   skip, machine ); % 32-35: image size
      fwrite( fh,      0,'float32',   skip, machine ); % 36-39: stroke width
      fwrite( fh,      0, 'uint32',   skip, machine ); % 40-43: roi properties offset
      fwrite( fh,      0, 'uint32',   skip, machine ); % 44-47: roi properties params
      fwrite( fh,      0, 'uint32',   skip, machine ); % 48-51: counters offset

      % fissa roi files have 3 extra long ints but I've no idea why ?!?!
      while offset>0
         fwrite( fh,      0, 'uint32',   skip, machine ); % 52-55: counters offset
         offset = offset - 1;
      end
      while ~isempty( fname )
         l = length( fname );
         fwrite( fh, fname(1:min(l,2)),  'int16',   skip, machine ); % 0-3
         if length( fname ) <=2
            fname = [];
         else
            fname = fname(3:end);
         end
      end

      fclose( fh );
      
   catch ME
      fclose( fh );
      
      str = getCatchMEstring( ME, 'Error writing ROI: ', true );
      cprintf( 'Keywords*', str );
      return
   end
   
   err = false;

end
 
   % ## ImageJ/NIH Image 64 byte ROI outline header
   % ## - 2 byte numbers are big-endian signed shorts
   % ##     0-3     "Iout"
   % ##     4-5     version (>=217)
   % ##     6-7     roi type
   % ##     8-9     top
   % ##     10-11   left
   % ##     12-13   bottom
   % ##     14-15   right
   % ##     16-17   NCoordinates
   % ##     18-33   x1,y1,x2,y2 (straight line) | x,y,width,height (double rect) | size (npoints)
   % ##     34-35   stroke width (v1.43i or later)
   % ##     36-39   ShapeRoi size (type must be 1 if this value>0)
   % ##     40-43   stroke color (v1.43i or later)
   % ##     44-47   fill color (v1.43i or later)
   % ##     48-49   subtype (v1.43k or later)
   % ##     50-51   options (v1.43k or later)
   % ##     52-52   arrow style or aspect ratio (v1.43p or later)
   % ##     53-53   arrow head size (v1.43p or later)
   % ##     54-55   rounded rect arc size (v1.43p or later)
   % ##     56-59   position
   % ##     60-63   reserved (zeros)
   % ##     64-     x-coordinates (short), followed by y-coordinates (shortz
   % ##
   
   % fwrite options: https://au.mathworks.com/help/matlab/ref/fwrite.html
   %  precision: 
   %    'float', 'float32', 'float64'
   %    'integer', 'integer*4' for 4-byte integers (*i for i=1,2,3,4)
   %    'real', 'real*4', 'real*8'
   %    'single', 'double'
   %    'uint', 'uint8', 'uint16', 'int8', 'int16', 'int32', 'int64'
   %    'uchar', 'schar', 'char', 'charl'
   %    'ushort', 'ulong', 'ubit<n>', 'short', 'long'
   %  skip: if precision is in bits, specify skip in bits
   %  machinefmt: 
   %    'ieee-be' - big endian
   %    'ieee-be.l64' - big endian, 64-bit, long data-type
   %    'ieee-le.l64' - little endian, 64-bit, long data-type
   %    'n'    = native
   %    'b'    = big endian
   %    'l'    = little endian

