function SpatialEdges = getEnvEdgePrior(occMapR,mode)
    SpatialPriorEnv = occMapR>0; % initialise prior to be the shape of occupied pixels
    % fill the holes in the explored regions
    if mode==1 % environment is doughnut-shaped
        SpatialPriorEnv = bwmorph(SpatialPriorEnv,'close',2); % do it twice to
        SpatialPriorEnv = bwmorph(SpatialPriorEnv,'dilate',2); % fill everything
        % Potentially manually select the holes to be filled?
        % SpatialPriorEnv = imfill(SpatialPriorEnv);
    elseif mode==2 % environment is a closed-fill space, highly sampled (GP)
        SpatialPriorEnv = bwmorph(SpatialPriorEnv,'dilate',1);
        SpatialPriorEnv = imfill(SpatialPriorEnv,'Holes');
        SpatialPriorEnv = bwmorph(SpatialPriorEnv,'close',1);
     else % environment is a closed-fill space, downsampled (hist)
        SpatialPriorEnv = imfill(SpatialPriorEnv,'Holes');
    end
    SpatialEdges = SpatialPriorEnv;
end