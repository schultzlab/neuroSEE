% [stack_g, stack_r, output, fh] = motionCorrectToNearestPixel(stack_g, stack_r, filedir, file, scale, Nimg_ave)
% Register both red & green channels of 2-photon calcium imaging data to a
% mean image to reduce motion artifact. Uses the first 200 frames to make
% an initial mean image, which is then updated as frames are motion
% corrected. The first 400 images are then re-corrected after the first pass 
% through, once the template image has improved. 
%
% INPUTS
%   stack_g     : matrix of green image stack
%   stack_r     : matrix of red image stack
%   file        : part of file name of image stacks in the format
%                   yyyymmdd_HH_MM_SS
%   Nimg_ave    : number of images to be averaged when correcting for pixel
%                           shift (zippering)
%   scale       : downsampling factor for motion correction
%   force       : if =1, motion correction will be done even though motion
%                   corrected images already exist
% OUTPUTS
%   stack_g     : matrix of motion corrected green image stack
%   stack_r     : matrix of motion corrected red image stack
%   output      : cell array containing
%                   output.settings.imscale & output.settings.Nimg_ave
%                   output.green & output.red
%                   each contains
%                       output.<colour>.template : template used for registration
%                       output.<colour>.meanframe : mean frame for original stack 
%                       output.<colour>.meanregframe : mean frame for registered stack
%                       output.<colour>.shift : matrix of x&y shift for each frame
%   fh          : figure summarising comparison between mean frames for original and
%                   registered stacks for green and red channels

% Written by Katie Davey
% Edit by Ann: 

function [stack_g, stack_r, out_g, out_r, col_shift, shifts, template, fh] = motionCorrectToNearestPixel(stack_g, stack_r, file, scale, Nimg_ave, refChannel, redoT, doplot)
   if nargin<8, doplot = 0; end
   if nargin<7, redoT = 300; end
   if nargin<6, refChannel = 'green'; end
   if nargin<5, Nimg_ave = 14; end
   if nargin<4, scale  = 1;    end

   tempfn = @(x) mean(x,3); % fn to convert buffer to template with - can try using median or mode (mode is bad!)
   initT  = 200; % use this many frames to make initial buffer with

   % if memory mapping, the memory mapped variable is in red, & green's empty
   % - make a copy of the memory mapped file in green to call for each channel
   usememmap = ternaryOp( strcmpi( 'memmapfile', class( stack_g ) ), true, false );
   
   % first dezipper images
   [ stack_g, stack_r, col_shift ] = doDezippering( stack_g, stack_r, Nimg_ave );
   str = sprintf( '%s: Dezippering done\n', file );
   cprintf( 'Text', str );

   % gotta register for green & for red
   %tic; 
   
   % use green shifts to register red frames since movement must be the
   % same for both channels, but green snr is much better, even though the
   % green channel is non-stationary while the red channel is supposed to
   % be stationary 
   if strcmpi(refChannel,'green')
       str = sprintf( '%s: Registering green channel\n', file );
       cprintf( 'Text', str );
       [stack_g, out_g, shifts, template] = registerChannel( initT, redoT, stack_g, tempfn, scale, usememmap  );
       str = sprintf( '%s: Registering red channel\n', file );
       cprintf( 'Text', str );
       [stack_r, out_r, ~, ~] = registerChannel( initT, redoT, stack_r, tempfn, scale, usememmap, shifts );
   else % refChannel = 'red'
       str = sprintf( '%s: Registering red channel\n', file );
       cprintf( 'Text', str );
       [stack_r, out_r, shifts, template] = registerChannel( initT, redoT, stack_r, tempfn, scale, usememmap  );
       str = sprintf( '%s: Registering green channel\n', file );
       cprintf( 'Text', str );
       [stack_g, out_g, ~, ~] = registerChannel( initT, redoT, stack_g, tempfn, scale, usememmap, shifts );
   end
      
   %printTime( toc, 'Running registration took' );

   if doplot
       fh = figure; 
       subplot(221), 
         imagesc( out_g.meanframe ); 
         axis image; colorbar; axis off;
         title( 'Mean frame for raw green' );
       subplot(222), 
         imagesc( out_g.meanregframe ); 
         axis image; colorbar; axis off; 
         title( 'Mean frame for corrected green' );
       subplot(223), 
         imagesc( out_r.meanframe ); 
         axis image; colorbar; axis off; 
         title( 'Mean frame for raw red' );
       subplot(224), 
         imagesc( out_r.meanregframe ); 
         axis image; colorbar; axis off;
         title( 'Mean frame for corrected red' );
       axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0  1],'Box','off',...
         'Visible','off','Units','normalized', 'clipping' , 'off');
   else
       fh = [];
   end
end

% register channel by either calculating highest correlation between each
% frame & the current mean frame, or else by applying given shifts
function [regch, output, shifts, template] = registerChannel( initT, redoT, channel, tempfn, scale, usememmap, shifts )
   % to apply green shift to red channel you can input your own shift
   % matrix, which is [dx dy] vectors merged 
   if nargin<7 || isempty(shifts)
      calcshift = true; 
   else
      calcshift = false;  
   end
   zeromean  = true; 
   
   % reading whole calcium image into memory at once
   if ~usememmap
      [szX, szY, maxT] = size( channel );
      
   % memory mapped file --> only read one frame into memory at one time
   else
      info   = imfinfo( channel.Filename );
      szX    = info.Width;
      szY    = info.Height;
      maxT   = numel( info ) / 2; % info is an array with 1 entry per tiff stack (interleaved green/red) 
   end
   pixPerFrame = szX * szY; 

   % if we're calculating shift amounts in x & y direction get initial
   % buffer full of frames to calculate initial template from
   if calcshift
      shifts = zeros( maxT, 2 );   % record change in x & y for each frame
      
      % get initial frames to make initial template from
      if usememmap
      % only need a buffer to calc template from if calculating shifts
         % for mem mapped file, get only enough for initial template image
         buff = zeros( szX, szY, initT );
         for i=1:initT
            buff(:,:,i) = double( reshape( channel.Data( (pixPerFrame*(i-1))+1 : (pixPerFrame*i) ), [szX szY] )' );
         end
         
      else
         buff = channel(:,:,1:initT); % frame buffer to calculate template from
      end
      
      template  = tempfn( buff );     % get init template to register frames to
      tempmean  = ternaryOp( zeromean, mean( template(:) ), 0 ); % zero mean or not
      R         = fft2( template - tempmean ); % fourier of registration template
      
   end

   % get indices with image centre at (0,0) so shift is calculated relative to centre
   indx      = ifftshift( -fix(szX/2) : ceil(szX/2)-1 ); % 2D freq indices at original resolution
   indy      = ifftshift( -fix(szY/2) : ceil(szY/2)-1 );
   [freqy, freqx] = meshgrid(indy,indx); % mesh of freq indices centred at 0
   freqx     = freqx / szX; % proportion of image traversed in x direction
   freqy     = freqy / szY;
   indx_us   = ifftshift( -fix(szX/2*scale) : ceil(szX/2*scale)-1 ); % upsampled indices
   indy_us   = ifftshift( -fix(szY/2*scale) : ceil(szY/2*scale)-1 );
   tvec      = [1:maxT 1:redoT];         % redo initial frames once out registration img is better

   regch     = zeros( szX, szY, maxT ); % registered channel   
   meanframe = zeros( szX, szY ); 
   meanregframe = zeros( szX, szY ); 
   ii=0; tt = length(tvec); prevstr = [];
   for ti=tvec
      ii = ii+1;
      if ( mod( round((ii-1)/tt*10000), 1000 ) == 0 ) || ( ii==length(tvec) )
          % total times through loop
         str = sprintf( '\t%g percent complete...\n', round( ii/tt*100 ) );
         refreshdisp( str, prevstr );
         prevstr = str; 
      end

      % register green frame, update green registration buffer & calculate template
      if usememmap
      % only need a buffer to calc template from if calculating shifts
         frame     = double( reshape( channel.Data( (pixPerFrame*(ti-1))+1 : (pixPerFrame*ti) ), [szX szY] )' );
      else
         frame     = channel(:,:,ti);
      end

      meanF        = ternaryOp( zeromean, mean( frame(:) ), 0 );
      F            = fft2( frame - meanF );    % max cross-corr calculated in freq domain
      if calcshift
         % calculate shift in x & y directions that maximises cross-correlation
         [dx, dy]  = registerFrame( R, F, scale, indx_us, indy_us ); 
         shifts(ti,:) = [dx dy];
      else
         % use shifts provided 
         dx = shifts(ti,1); dy = shifts(ti,2);
      end
      % apply shifts to frame to get the registered frame
      regF         = calcRegisteredFrame( F, freqx, freqy, dx, dy );

      regframe     = real( ifft2( regF ) + meanF );
      regch(:,:,ti)= regframe;          % populate corrected green frames
      
      % if calculating shift amounts update the buffer with the new frame
      if calcshift
         buff      = updateBuffer( initT, ti, buff, regframe );
         template  = tempfn( buff );         % update green template
         tempmean  = ternaryOp( zeromean, mean( template(:) ), 0 );
         R         = fft2( template - tempmean );% get fourier of green registration template
      end
      
      % calculate mean frames & corrected frames as we go along to display
      % to user at the end for integrity testing
      meanframe    =    meanframe +    frame / maxT; 
      meanregframe = meanregframe + regframe / maxT; 

   end % end for each image frame
   
   output.meanframe    = meanframe; 
   output.meanregframe = meanregframe;
   if ~calcshift
       template = [];
   end
end

function buffer = updateBuffer( initT, ti, buffer, regframe )
   % update mean image to register to
   mind = mod( ti, initT ); % index into mean img to update - only keep last initT images
   mind = ternaryOp( mind, mind, initT ); % if mean index is 0, set to last img
   buffer(:,:,mind) = regframe;  % update buffer to calc template from
end

% input fourier domain version of registration template (R), and frame to
% be registered (F), and return fourier domain version of registered frame
% (regF)
function [dx, dy] = registerFrame( R, F, scale, indx_us, indy_us )
   C = R .* conj(F); % cross-corr in fourier domain
   if scale>1
      C  = fftInterpolate( C, scale*[szX szY], true );
   end
   c = ifft2( C ); % cross-corr in spatial domain

   % find position of highest cross-corr btwn reg img & frame
   [dx, dy] = find( abs(c) == max(abs(c(:))), 1, 'first' );

   % get position of highest cross-corr with centre of img at (0,0)
   dx   = indx_us(dx) / scale;
   dy   = indy_us(dy) / scale;
end

function regF = calcRegisteredFrame( F, freqx, freqy, dx, dy )
   % calculate registered frame
   regF = F .* exp( 1i*2*pi * (-dx*freqx - dy*freqy) );
end

