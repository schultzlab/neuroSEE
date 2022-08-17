%% USER INPUTS
list = 'list_m123_fam1fam2fam1.txt'; 
n = 5; % SPECIFY FILE # to register
refind = 11; % SPECIFY FILE # for reference 
max_dev = 7; % I normally tweak this. I found that 7-13 is usually good enough.
% ----------------------

[data_locn,comp,err] = load_neuroSEEmodules;
if ~isempty(err)
    beep
    cprintf('Errors',err);    
    return
end

listfile = [data_locn 'Digital_Logbook/lists_imaging/' list];
files = extractFilenamesFromTxtfile( listfile );

for i = 1:size(files,1)
    file = files(i,:);
    load([data_locn 'Data/' file(1:8) '/Processed/' file '/mcorr_normcorre/' file '_mcorr_output.mat'])
    imG(:,:,i) = green.meanregframe;
    imR(:,:,i) = red.meanregframe;
    clear red green shifts col_shift template params
end

g = 64; % this is half the default value for normcorre
params = neuroSEE_setparams(...
            'max_shift_r', 30,...       
            'max_shift_nr', 30,...
            'grid_size_nr', [g,g],...
            'iter',1,...
            'max_dev',max_dev,... 
            'overlap_pre', [g/4,g/4],... 
            'min_patch_size', [g/4,g/4],...      
            'min_diff', [g/8,g/8]);         
        
regchannel = 1; % 1:green, 2:red

if regchannel == 1
    Y = imG(:,:,n); X = imR(:,:,n);
    ref_ch1 = imG(:,:,refind); ref_ch2 = imR(:,:,refind);
    titlestr1 = ['Image ' num2str(n) ' reg to ' num2str(refind) ' (green ch)'];
    titlestr2 = ['Image ' num2str(n) ' reg to ' num2str(refind) ' (red ch)'];
else
    Y = imR(:,:,n); X = imG(:,:,n);
    ref_ch1 = imR(:,:,refind); ref_ch2 = imG(:,:,refind);
    titlestr1 = ['Image ' num2str(n) ' reg to ' num2str(refind) ' (red ch)'];
    titlestr2 = ['Image ' num2str(n) ' reg to ' num2str(refind) ' (green ch)'];
end

[YY, shifts_r, ~, options_r, col_shift_r] = normcorre(Y,params.mcorr.normcorre_r,ref_ch1);
[YY, shifts_nr, ~,options_nr,col_shift_nr] = normcorre(YY,params.mcorr.normcorre_nr,ref_ch1);
figure;
fy = imfuse( YY, ref_ch1, 'falsecolor', 'Scaling', 'joint', 'ColorChannels', [1 2 0]);
imshow(fy); title(titlestr1);

