% Written by Ann Go

function globalregmasks = registerMasks( masks, dx, dy)
   Nrows = size(masks,1);
   Ncols = size(masks,2);
   % register rows
   if dy > 0
      shift = dy;
   	  masks( :, (1+shift):Ncols, : ) = masks( :, 1:(Ncols-shift), : );
      masks( :, 1:shift,         : ) = 0;
   else   % if dx < 0
      shift = abs( dy );
      masks( :, 1:(Ncols-shift)    , : ) = masks( :, (1+shift):Ncols, : );
      masks( :, Ncols-shift+1:Ncols, : ) = 0;
   end
   % register columns
   if dx > 0
      shift = dx;
   	  masks( (1+shift):Nrows, :, : ) = masks( 1:(Nrows-shift), :, : );
      masks( 1:shift,         :, : ) = 0;
   else   % if dy < 0
      shift = abs( dx );
      masks( 1:(Ncols-shift)    , :, : ) = masks( (1+shift):Ncols, :, : );
      masks( Ncols-shift+1:Ncols, :, : ) = 0;
   end
   globalregmasks = masks;
end