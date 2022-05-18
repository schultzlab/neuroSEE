numcells =  numel(m62a.pfData.spkRaster) + ...
            numel(m62b.pfData.spkRaster) + ...
            numel(m66.pfData.spkRaster) + ...
            numel(m70.PFdata.spkRaster) + ...
            numel(m111.PFdata.spkRaster) + ...
            numel(m116.PFdata.spkRaster) + ...
            numel(m118a.PFdata.spkRaster) + ...
            numel(m118b.PFdata.spkRaster) + ...
            numel(m118c.PFdata.spkRaster) + ...
            numel(m125.PFdata.spkRaster);

%m62a
for i = 1:numel(m62a.pfData.spkRaster) 
    cell.SIsec{ i } = m62a.hist.infoMap(i,1);
    cell.SIspk{ i } = m62a.hist.infoMap(i,2);
    cell.mouseID{ i } = 'm62';
    cell.genotype{ i } = 'WT';
    cell.ageMOS{ i } = 5.6;
    cell.sex{ i } = 'm';
    if any(ismember(i,m62a.hist.SIsec.pcIdx))
        cell.pc_SIsec{ i } = 'pc';
    else
        cell.pc_SIsec{ i } = 'nonpc';
    end
    if any(ismember(i,m62a.hist.SIspk.pcIdx))
        cell.pc_SIspk{ i } = 'pc';
    else
        cell.pc_SIspk{ i } = 'nonpc';
    end
end
%m62b
k = numel(cell.sex);
for i = 1:numel(m62b.pfData.spkRaster) 
    cell.SIsec{ k + i } = m62b.hist.infoMap(i,1);
    cell.SIspk{ k + i } = m62b.hist.infoMap(i,2);
    cell.mouseID{ k + i } = 'm62';
    cell.genotype{ k + i } = 'WT';
    cell.ageMOS{ k + i } = 5.6;
    cell.sex{ k + i } = 'm';
    if any(ismember(i,m62b.hist.SIsec.pcIdx))
        cell.pc_SIsec{ k + i } = 'pc';
    else
        cell.pc_SIsec{ k + i } = 'nonpc';
    end
    if any(ismember(i,m62b.hist.SIspk.pcIdx))
        cell.pc_SIspk{ k + i } = 'pc';
    else
        cell.pc_SIspk{ k + i } = 'nonpc';
    end
end
%m66
k = numel(cell.sex);
for i = 1:numel(m66.pfData.spkRaster) 
    cell.SIsec{ k + i } = m66.hist.infoMap(i,1);
    cell.SIspk{ k + i } = m66.hist.infoMap(i,2);
    cell.mouseID{ k + i } = 'm66';
    cell.genotype{ k + i } = 'WT';
    cell.ageMOS{ k + i } = 9.3;
    cell.sex{ k + i } = 'f';
    if any(ismember(i,m66.hist.SIsec.pcIdx))
        cell.pc_SIsec{ k + i } = 'pc';
    else
        cell.pc_SIsec{ k + i } = 'nonpc';
    end
    if any(ismember(i,m66.hist.SIspk.pcIdx))
        cell.pc_SIspk{ k + i } = 'pc';
    else
        cell.pc_SIspk{ k + i } = 'nonpc';
    end
end
%m70
k = numel(cell.sex);
for i = 1:numel(m70.PFdata.spkRaster) 
    cell.SIsec{ k + i } = m70.hist.infoMap(i,1);
    cell.SIspk{ k + i } = m70.hist.infoMap(i,2);
    cell.mouseID{ k + i } = 'm70';
    cell.genotype{ k + i } = 'WT';
    cell.ageMOS{ k + i } = 8.6;
    cell.sex{ k + i } = 'm';
    if any(ismember(i,m70.hist.SIsec.pcIdx))
        cell.pc_SIsec{ k + i } = 'pc';
    else
        cell.pc_SIsec{ k + i } = 'nonpc';
    end
    if any(ismember(i,m70.hist.SIspk.pcIdx))
        cell.pc_SIspk{ k + i } = 'pc';
    else
        cell.pc_SIspk{ k + i } = 'nonpc';
    end
end
%m111
k = numel(cell.sex);
for i = 1:numel(m111.PFdata.spkRaster) 
    cell.SIsec{ k + i } = m111.hist.infoMap(i,1);
    cell.SIspk{ k + i } = m111.hist.infoMap(i,2);
    cell.mouseID{ k + i } = 'm111';
    cell.genotype{ k + i } = '5xFAD';
    cell.ageMOS{ k + i } = 7.5;
    cell.sex{ k + i } = 'f';
    if any(ismember(i,m111.hist.SIsec.pcIdx))
        cell.pc_SIsec{ k + i } = 'pc';
    else
        cell.pc_SIsec{ k + i } = 'nonpc';
    end
    if any(ismember(i,m111.hist.SIspk.pcIdx))
        cell.pc_SIspk{ k + i } = 'pc';
    else
        cell.pc_SIspk{ k + i } = 'nonpc';
    end
end
%m116
k = numel(cell.sex);
for i = 1:numel(m116.PFdata.spkRaster) 
    cell.SIsec{ k + i } = m116.hist.infoMap(i,1);
    cell.SIspk{ k + i } = m116.hist.infoMap(i,2);
    cell.mouseID{ k + i } = 'm116';
    cell.genotype{ k + i } = 'WT';
    cell.ageMOS{ k + i } = 7.6;
    cell.sex{ k + i } = 'f';
    if any(ismember(i,m116.hist.SIsec.pcIdx))
        cell.pc_SIsec{ k + i } = 'pc';
    else
        cell.pc_SIsec{ k + i } = 'nonpc';
    end
    if any(ismember(i,m116.hist.SIspk.pcIdx))
        cell.pc_SIspk{ k + i } = 'pc';
    else
        cell.pc_SIspk{ k + i } = 'nonpc';
    end
end
%m118a
k = numel(cell.sex);
for i = 1:numel(m118a.PFdata.spkRaster) 
    cell.SIsec{ k + i } = m118a.hist.infoMap(i,1);
    cell.SIspk{ k + i } = m118a.hist.infoMap(i,2);
    cell.mouseID{ k + i } = 'm118';
    cell.genotype{ k + i } = '5xFAD';
    cell.ageMOS{ k + i } = 9.1;
    cell.sex{ k + i } = 'm';
    if any(ismember(i,m118a.hist.SIsec.pcIdx))
        cell.pc_SIsec{ k + i } = 'pc';
    else
        cell.pc_SIsec{ k + i } = 'nonpc';
    end
    if any(ismember(i,m118a.hist.SIspk.pcIdx))
        cell.pc_SIspk{ k + i } = 'pc';
    else
        cell.pc_SIspk{ k + i } = 'nonpc';
    end
end
%m118b
k = numel(cell.sex);
for i = 1:numel(m118b.PFdata.spkRaster) 
    cell.SIsec{ k + i } = m118b.hist.infoMap(i,1);
    cell.SIspk{ k + i } = m118b.hist.infoMap(i,2);
    cell.mouseID{ k + i } = 'm118';
    cell.genotype{ k + i } = '5xFAD';
    cell.ageMOS{ k + i } = 9.1;
    cell.sex{ k + i } = 'm';
    if any(ismember(i,m118b.hist.SIsec.pcIdx))
        cell.pc_SIsec{ k + i } = 'pc';
    else
        cell.pc_SIsec{ k + i } = 'nonpc';
    end
    if any(ismember(i,m118b.hist.SIspk.pcIdx))
        cell.pc_SIspk{ k + i } = 'pc';
    else
        cell.pc_SIspk{ k + i } = 'nonpc';
    end
end
%m118c
k = numel(cell.sex);
for i = 1:numel(m118c.PFdata.spkRaster) 
    cell.SIsec{ k + i } = m118c.hist.infoMap(i,1);
    cell.SIspk{ k + i } = m118c.hist.infoMap(i,2);
    cell.mouseID{ k + i } = 'm118';
    cell.genotype{ k + i } = '5xFAD';
    cell.ageMOS{ k + i } = 9.1;
    cell.sex{ k + i } = 'm';
    if any(ismember(i,m118c.hist.SIsec.pcIdx))
        cell.pc_SIsec{ k + i } = 'pc';
    else
        cell.pc_SIsec{ k + i } = 'nonpc';
    end
    if any(ismember(i,m118c.hist.SIspk.pcIdx))
        cell.pc_SIspk{ k + i } = 'pc';
    else
        cell.pc_SIspk{ k + i } = 'nonpc';
    end
end
%m125
k = numel(cell.sex);
for i = 1:numel(m125.PFdata.spkRaster) 
    cell.SIsec{ k + i } = m125.hist.infoMap(i,1);
    cell.SIspk{ k + i } = m125.hist.infoMap(i,2);
    cell.mouseID{ k + i } = 'm125';
    cell.genotype{ k + i } = '5xFAD';
    cell.ageMOS{ k + i } = 6.5;
    cell.sex{ k + i } = 'm';
    if any(ismember(i,m125.hist.SIsec.pcIdx))
        cell.pc_SIsec{ k + i } = 'pc';
    else
        cell.pc_SIsec{ k + i } = 'nonpc';
    end
    if any(ismember(i,m125.hist.SIspk.pcIdx))
        cell.pc_SIspk{ k + i } = 'pc';
    else
        cell.pc_SIspk{ k + i } = 'nonpc';
    end
end


cell.SIsec = cell.SIsec';
cell.SIspk = cell.SIspk';
cell.mouseID = cell.mouseID';
cell.genotype = cell.genotype';
cell.ageMOS = cell.ageMOS';
cell.sex = cell.sex';
cell.pc_SIsec = cell.pc_SIsec';
cell.pc_SIspk = cell.pc_SIspk';

% 165
% 165+256 = 421
% 165+256+88 = 509
% 165+256+88+54 = 563
% 165+256+88+54+348 = 911
% 165+256+88+54+348+225 = 1136
% 165+256+88+54+348+225+190 = 1326
% 165+256+88+54+348+225+190+188 = 1514
% 165+256+88+54+348+225+190+188+61 = 1575
% 165+256+88+54+348+225+190+188+61+160 = 1735

