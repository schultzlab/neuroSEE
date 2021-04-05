function filetime = extractTimeFromFilename( timedir, format )
   if length( timedir ) ~= length( format )
%       filetime = [];
%       return;
        format = ['Saved' format];
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