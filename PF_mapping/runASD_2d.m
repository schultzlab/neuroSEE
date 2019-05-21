function [ASD,ASDstats] = runASD_2d(x,z,dims,mask,ml)
    % -- inputs --
    % x:        the flattened and binned 2d directtory
    % z:        the reponse array, sampled like x
    % dims:     a 2x1 array containing the dimensions of the 2D env
    % mask:     spatial mask of the environment
    % ml:       the minimum lengthscale for ASD, default=100
    % -- outputs --
    % kasd:     smooth ASD estimaion
    % ASDstats: extra prameters and stats of ASD estimation
    
    if nargin<4; MASK = 0; ml = 100; end
    if nargin<5; MASK = 1; ml = 100;
    elseif sum(mask(:)); MASK = 1;
    else; MASK = 0; end
    
    xx = zeros(length(x),prod(dims));
    for it=1:length(x)
        xx(it,x(it)) = 1;
    end
    % run ASD
    minlens = [dims(1)/ml,dims(2)/ml];  % minimum length scale along each dimension
    [kasd,ASDstats] = fastASD(xx,z,dims,minlens);
    %kasd(kasd<0) = 0; % get rid of possible negative estimates  
    if MASK; kasd(~mask) = 0; end
    ASD = reshape(kasd,dims);
end