function [ regstack_g, regstack_r, regmasks] = globalregisterStack( stack_g, stack_r, dx, dy )
    szX     = size(stack_g,1);
    szY     = size(stack_g,2);
    nT      = size(stack_g,3);
    zeromean  = true; 


    indx      = ifftshift( -fix(szX/2) : ceil(szX/2)-1 ); % 2D freq indices at original resolution
    indy      = ifftshift( -fix(szY/2) : ceil(szY/2)-1 );
    [freqy, freqx] = meshgrid(indy,indx); % mesh of freq indices centred at 0
    freqx     = freqx / szX; % proportion of image traversed in x direction
    freqy     = freqy / szY;
    
    regstack_g  = zeros( szX, szY, nT );
    regstack_r  = zeros( szX, szY, nT );
    for i=1:nT
        frame_g             = stack_g(:,:,i);
        meanF_g             = ternaryOp( zeromean, mean( frame_g(:) ), 0 );
        F_g                 = fft2( frame_g - meanF_g );    % max cross-corr calculated in freq domain
        regF_g              = calcRegisteredFrame( F_g, freqx, freqy, dx, dy );
        regstack_g(:,:,i)   = real( ifft2( regF_g ) + meanF_g );  
        
        frame_r             = stack_r(:,:,i);
        meanF_r             = ternaryOp( zeromean, mean( frame_r(:) ), 0 );
        F_r                 = fft2( frame_r - meanF_r );    % max cross-corr calculated in freq domain
        regF_r              = calcRegisteredFrame( F_r, freqx, freqy, dx, dy );
        regstack_r(:,:,i)   = real( ifft2( regF_r ) + meanF_r );          
    end
end
    
function regF = calcRegisteredFrame( F, freqx, freqy, dx, dy )
   % calculate registered frame
   regF = F .* exp( 1i*2*pi * (-dx*freqx - dy*freqy) );
end
