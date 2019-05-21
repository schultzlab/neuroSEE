% Written by Ann Go

function [ shift, regImage ] = globalregisterImage(template, imagetoreg, scale )
   if nargin<3, scale  = 1;    end
    
    % template = double( reshape( channel.Data( (pixPerFrame*(ti-1))+1 : (pixPerFrame*ti) ), [szX szY] )' );
    zeromean  = true; 
    tempmean  = ternaryOp( zeromean, mean( template(:) ), 0 ); % zero mean or not
    R         = fft2( template - tempmean ); % fourier of registration template

    szX    = size(imagetoreg,1);
    szY    = size(imagetoreg,2);

    indx      = ifftshift( -fix(szX/2) : ceil(szX/2)-1 ); % 2D freq indices at original resolution
    indy      = ifftshift( -fix(szY/2) : ceil(szY/2)-1 );
    [freqy, freqx] = meshgrid(indy,indx); % mesh of freq indices centred at 0
    freqx     = freqx / szX; % proportion of image traversed in x direction
    freqy     = freqy / szY;
    indx_us   = ifftshift( -fix(szX/2*scale) : ceil(szX/2*scale)-1 ); % upsampled indices
    indy_us   = ifftshift( -fix(szY/2*scale) : ceil(szY/2*scale)-1 );

    meanF        = ternaryOp( zeromean, mean( imagetoreg(:) ), 0 );
    F            = fft2( imagetoreg - meanF );    % max cross-corr calculated in freq domain
    % calculate shift in x & y directions that maximises cross-correlation
    [dx, dy]  = registerFrame( R, F, scale, indx_us, indy_us ); 
    shift = [dx dy];
    regF         = calcRegisteredFrame( F, freqx, freqy, dx, dy );
    regImage     = real( ifft2( regF ) + meanF );      
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

