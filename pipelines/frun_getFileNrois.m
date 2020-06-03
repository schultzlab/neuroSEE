function Nrois = frun_getFileNrois( list )

mcorr_method = 'normcorre-nr';     % values: [normcorre, normcorre-r, normcorre-nr, fftRigid] 
segment_method = 'CaImAn';      % [ABLE,CaImAn]    

% Load module folders and define data directory
[data_locn, ~, err] = load_neuroSEEmodules(false);
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
        Nrois(n) = size(M.masks,3);
    end

end