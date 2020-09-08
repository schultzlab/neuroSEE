g = 64;
params = neuroSEE_setparams(...
            'mcorr_method', mcorr_method, ...
            'segment_method', segment_method,...
            'runpatches', runpatches,...
            'dofissa', dofissa, ...
            'doasd', doasd,...
            'max_shift_r', 80,...       
            'max_shift_nr', 30,...
            'grid_size_nr', [g,g],...
            'iter',1,...
            'max_dev',8,...     
            'overlap_pre', [g/4,g/4],... 
            'min_patch_size', [g/4,g/4],...      
            'min_diff', [g/8,g/8]);             
           

Y = A1; ref = A;     
Ystr = 'A1'; refstr = 'B';
titlestr = 'A1 registered to B';

[Yar,A1toBshifts_r,~,A1toBoptions_r,A1toBcol_shift_r] = normcorre(Y,params.mcorr.normcorre_r,ref);
[Ya,A1toBshifts_nr,~,A1toBoptions_nr,A1toBcol_shift_nr] = normcorre(Yar,params.mcorr.normcorre_nr,ref);
% [Ca2,~,~,~,~] = normcorre(C,params.mcorr.normcorre_nr,A);
%     figure;
%     CA1 = imfuse( Car, A, 'falsecolor', 'Scaling', 'joint', 'ColorChannels', [1 2 0]);
%     imshow(CA1); title('C registered to A (rigid)');
    figure;
    CY2 = imfuse( Ya, ref, 'falsecolor', 'Scaling', 'joint', 'ColorChannels', [1 2 0]);
    imshow(CY2); title(titlestr);
    figure;
    CY2 = imfuse( Ya, ref, 'falsecolor', 'Scaling', 'joint', 'ColorChannels', [2 1 0]);
    imshow(CY2); title(titlestr);

%     figure;
%     CA3 = imfuse( Ca2, A, 'falsecolor', 'Scaling', 'joint', 'ColorChannels', [1 2 0]);
%     imshow(CA3);  title('C registered to A (nonrigid)');

figure;
subplot(221); imagesc(Y); axis square; axis off; colormap(gray); title(Ystr);
subplot(222); imagesc(ref); axis square; axis off; colormap(gray); title(refstr);
subplot(223); CY2 = imfuse( Ya, ref, 'falsecolor', 'Scaling', 'joint', 'ColorChannels', [1 2 0]);
    imshow(CY2);  title(titlestr);
subplot(224); CY2 = imfuse( Ya, ref, 'falsecolor', 'Scaling', 'joint', 'ColorChannels', [2 1 0]);
    imshow(CY2);  title(titlestr);
