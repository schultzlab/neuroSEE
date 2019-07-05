    globalmasks = zeros(512,512,162);
    [cvsROIs] = ReadImageJROI('RoiSet.zip');
    [x] = ROIs2Regions(cvsROIs, [512 512]);
    
    for i = 1:162
        globalmasksInd{i} = masks.PixelIdxList{i};
    end
    
    for i = 1:180
        masks()
    end
    
    [sRegions] = ROIs2Regions(Coor, [512 512]);