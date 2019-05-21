tic
tiffile = fname_tif_gr;
images = imG;
datatype = 'uint16';
InfoImage = imfinfo( tiffile );
   wdth   = InfoImage(1).Width;
   hgt    = InfoImage(1).Height;
   numimg = length(InfoImage);
   fileh  = tifflib( 'open', tiffile, 'rw' );
   rps    = tifflib( 'getField', fileh, Tiff.TagID.RowsPerStrip );

   for i=1:numimg
      tifflib( 'setDirectory', fileh, i-1 ); % 0 indexing
      % Go through each strip of data.
      rps = min(rps, hgt );
      for r = 1:rps:hgt
         row_inds = r:min( hgt, r+rps-1 );
         stripnum = tifflib( 'computeStrip', fileh, r ) - eps; % failing when exactly 1
         tifflib( 'writeEncodedStrip', fileh, stripnum, images(row_inds,:,i) );
      end
   end
   tifflib( 'close', fileh );
   toc
   
   