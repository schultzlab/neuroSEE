function [kasd,ASDstats] = runASD_1d(x,z,dims)
    % z = z - mean(z);
    xx = zeros(length(x),dims);
    for it=1:length(x)
        xx(it,x(it)) = 1;
        % xx(it,:) = xx(it,:) - mean(xx(it,:));
    end
    % run ASD
    % minlens = round(dims/100);  % minimum length scale along each dimension
    minlens = ceil(dims/100);
    [kasd,ASDstats] = fastASD(xx,z,[dims,1],minlens);
    %kasd(kasd<0) = 0; % get rid of possible negative estimates    
end