% Written by Ann Go
% Script which takes as input a list of image files and outputs the
% number of ROIs detected in each file
% mcorr_method - the motion correction method used

function Nrois = frun_getFileNrois( list, mcorr_method )

if nargin<2, mcorr_method = 'normcorre'; end     % values: [normcorre, normcorre-r, normcorre-nr, fftRigid] 
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
        Nrois(n) = size(masks_all,3);
    else
        cprintf('Errors','No ROI segmentation output for %s.\n',file);
    end
end