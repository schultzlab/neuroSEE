clear 
[data_locn,comp,err] = load_neuroSEEmodules;
if ~isempty(err)
    beep
    cprintf('Errors',err);    
    return
end

list = 'list_m70_fam1fam2.txt'; % SPECIFY LIST
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
            'max_dev',9,... % I normally tweak this. I found that 7-9 is usually good enough.   
            'overlap_pre', [g/4,g/4],... 
            'min_patch_size', [g/4,g/4],...      
            'min_diff', [g/8,g/8]);         
        
%%           
n = 13; % SPECIFY FILE # to register
regchannel = 1; % 1:green, 2:red
refind = 1; % SPECIFY FILE # for reference 
%%

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
% [ YY, XX, ~, ~, col_shift_r, shifts_r, ~, ~ ] = normcorre_2ch( Y, X, params.mcorr.normcorre_r, ref_ch1 );
% [ YY, XX, ~, ~, col_shift_nr, shifts_nr, ~, ~ ] = normcorre_2ch( YY, XX, params.mcorr.normcorre_nr, ref_ch1 );
figure;
fy = imfuse( YY, ref_ch1, 'falsecolor', 'Scaling', 'joint', 'ColorChannels', [1 2 0]);
imshow(fy); title(titlestr1);

% clear params
% params.shifts = shifts_r; params.options = options_r; params.col_shift = col_shift_r;
% save('m133_r.mat','params')
% params.shifts = shifts_nr; params.options = options_nr; params.col_shift = col_shift_nr;
% save('m133_nr.mat','params')

