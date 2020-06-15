function Nrois = frun_getFileNrois( list, mcorr_method )

if nargin<2, mcorr_method = 'normcorre-nr'; end     % values: [normcorre, normcorre-r, normcorre-nr, fftRigid] 
segment_method = 'CaImAn';      % [ABLE,CaImAn]    

% Load module folders and define data directory
[data_locn, ~, err] = load_neuroSEEmodules;
if ~isempty(err)
    beep
    cprintf('Errors',err);    
    return
end

% list
listfile = [data_locn 'Digital_Logbook/lists_imaging/' list];
files = extractFilenamesFromTxtfile( listfile );

Nrois = zeros(size(files,1),1);
for n = 1:size(files,1)
    % image file
    file = files(n,:);

    dir_segment = [data_locn 'Data/' file(1:8) '/Processed/' file '/' 'mcorr_' mcorr_method '/' segment_method '/'];
    if exist([dir_segment file '_segment_output.mat'],'file')
        M = load([dir_segment file '_segment_output.mat']);
        masks_all = M.masks;
        % Eliminate very small rois and rois touching image border
        area = zeros(size(masks_all,3),1);
        borderpix = 3;
        for j = 1:size(masks_all,3)
            mask = masks_all(borderpix:size(masks_all,1)-borderpix,borderpix:size(masks_all,2)-borderpix,j);
            im = imclearborder(mask);
            c = regionprops(im,'area');
            if ~isempty(c)
                area(j) = c.Area;                    % area of each ROI
            end
        end
        Nrois(n) = length(area>roiarea_thr);
    end

end