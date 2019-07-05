function [stack_g, stack_r] = doDezippering( stack_g, stack_r, Navg )
% Function for generating image stacks corrected for Start Delay (pixel
% shift) in every second line
% INPUTS:
% stack_g : mat file containing image stack from green channel
% stack_r : mat file containing image stack from red channel
% Nave : number of images to be averaged for calculating pixel shift in
%        every second line of image stacks. I suggest using 10 or 14. 
% Note: Number of images in stack divided by Nave must be a whole number

% Written by Ann Go
% Edited by Katie Davey

   % determine if file is memory mapped or not
   usememmap = ternaryOp( strcmpi( 'memmapfile', class( stack_g ) ), true, false );

   % Load the dimensions of the green image stack 
   if ~usememmap
      [szX, szY, Nimg] = size( stack_g );
   else
      info = imfinfo( stack_g.Filename );
      szX  = info.Width;
      szY  = info.Height;
      Nimg = numel( info );    % 'info' is an array with 1 entry per tiff frame
      stack_g.Writable = true; % make sure we can write to the mem mapped file
      stack_r.Writable = true; 
   end

   % make sure that Navg is a factor of num images so that there's no
   % frames left out at the end
   factors = divisors( Nimg );
   if ~any( factors == Navg )
      % get the factor that's closest in value to Navg
      Ndiff = factors - Navg; 
      [~,i] = min( abs( Ndiff ) ); % get index of factor that's closest
      Navg  = factors( i );
      % if Navg is now 1, get the next factor up as that's stupid
      if Navg==1
         Navg = factors( 2 ); % i.e. find( closest) + 1 = 2 since 1 must be closest(1)
      end
      str   = sprintf( '\nNumber of images to be averaged is not a multiple of number of images, changing to %d\n', ...
                        Navg );
      cprintf( 'Comments', str );
   end

   % Take the average of Navg frames of stack_g and use them to calculate the
   % pixel shift in every second line
   pixPerFrame = szX * szY; pixPerStack = pixPerFrame * Navg; 
   
   % for each of the Navg mini-stacks, calc avg & then calc the pixel shift
   for i=1:Nimg/Navg
      % get the mini stack of Navg frames
      if usememmap % calc indices in num bytes
         % we want Navg frames from the g stack each time through the loop
         sind = pixPerStack*(i-1) + 1;
         eind = pixPerStack*i; 
         g = double( reshape( stack_g.Data(sind:eind), [szX szY Navg] ) );
         r = double( reshape( stack_r.Data(sind:eind), [szX szY Navg] ) );

      else % calc indices in num frames
         sind = (i-1)*Navg+1;
         eind = i*Navg;
         g = stack_g(:,:,sind:eind);
         r = stack_r(:,:,sind:eind);
      end
      
       % calc shift & correct the frames in the mini stack of Navg frames
      [g, r] = correctFrames( g, r, szX, szY );
      
      % write corrections to the big stacks
      if usememmap
         stack_g.Data(sind:eind) = g;
         stack_r.Data(sind:eind) = r;
      else
         stack_g(:,:,sind:eind)  = g; 
         stack_r(:,:,sind:eind)  = r;
      end
   end
end

function [ministack_g, ministack_r] = correctFrames( ministack_g, ministack_r, Nrows, Ncols )
   mean_g = mean( ministack_g, 3 ); % calculate the stack average
   % calculate shift and apply it to correct the mini stack
   shift = findBestShift( mean_g );
  if shift > 0
   	  ministack_g( 1:2:Nrows, (1+shift):Ncols, : ) = ministack_g( 1:2:Nrows, 1:(Ncols-shift), : );
      ministack_g( 1:2:Nrows, 1:shift,         : ) = 0;
      ministack_r( 1:2:Nrows, (1+shift):Ncols, : ) = ministack_r( 1:2:Nrows, 1:(Ncols-shift), : );
      ministack_r( 1:2:Nrows, 1:shift,         : ) = 0;
   end
   if shift < 0
      shift = abs( shift );
      ministack_g( 2:2:Nrows, (1+shift):Ncols, : ) = ministack_g( 2:2:Nrows, 1:(Ncols-shift), : );
      ministack_g( 2:2:Nrows, 1:shift,         : ) = 0;
      ministack_r( 2:2:Nrows, (1+shift):Ncols, : ) = ministack_r( 2:2:Nrows, 1:(Ncols-shift), : );
      ministack_r( 2:2:Nrows, 1:shift,         : ) = 0;
   end  
end

function shift = findBestShift( mean_g )
    maxshift  = 50; % max no. of pixel shift we will consider
    img1 = mean_g(1:2:end,:); 
    img2  = mean_g(2:2:end,:); 
    err_i = immse(img1,img2); 
    i = 1;
    
    % Find the pixel shift for which mean-squared error of the 2 sub-images 
    % is minimum
    % First, anchor img1 & shift img2 to the right
    while i < maxshift
        err_f = immse(img1(:,1:end-i),img2(:,i+1:end));
        if err_f > err_i
            i = i - 1;
            break
        else
            err_i = err_f;
            i = i + 1;
        end
    end
    
    if i == 0 % that is, if shifting img2 to the right didn't work
        i = -1;
        % try, anchoring img2 & shifting img1 to the right
        while abs(i) < maxshift
        err_f = immse(img2(:,1:end-abs(i)),img1(:,abs(i)+1:end));
            if err_f > err_i
                i = i + 1;
                break
            else
                err_i = err_f;
                i = i - 1;
            end
        end
    end
    shift = i;
end