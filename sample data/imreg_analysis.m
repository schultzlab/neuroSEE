mcorr_method = 'normcorre';
segment_method = 'CaImAn';      % [ABLE,CaImAn]    
runpatches = false;            % for CaImAn processing, flag to run patches (default: false)
dofissa = true;                 % flag to implement FISSA (when false, overrides force(3) setting)
doasd = false;                  % flag to do asd pf calculation

% Processing parameters (any parameter that is not set gets a default value)
% Add any parameters you want to set after FOV. See neuroSEE_setparams for
% full set of parameters
params = neuroSEE_setparams(...
            'mcorr_method', mcorr_method, ...
            'segment_method', segment_method,...
            'runpatches', runpatches,...
            'dofissa', dofissa, ...
            'doasd', doasd,...
            'grid_size_r', [64 64],...
            'grid_size_nr', [32 32]);

% A
[Abr,shifts_r,~,~,~] = normcorre(A,params.mcorr.normcorre_r,B);
[Ab,shifts,~,~,~] = normcorre(Abr,params.mcorr.normcorre_nr,B);
    figure;
    AB1 = imfuse( Abr, B, 'falsecolor', 'Scaling', 'joint', 'ColorChannels', [1 2 0]);
    imshow(AB1); title('A registered to B (rigid)');
    figure;
    AB2 = imfuse( Ab, B, 'falsecolor', 'Scaling', 'joint', 'ColorChannels', [1 2 0]);
    imshow(AB2); title('A registered to B');

[Acr,shifts_r,~,~,~] = normcorre(A,params.mcorr.normcorre_r,C);
[Ac,shifts,~,~,~] = normcorre(Acr,params.mcorr.normcorre_nr,C);
[Ac2,shifts,~,~,~] = normcorre(A,params.mcorr.normcorre_nr,C);
    figure;
    AC1 = imfuse( Acr, C, 'falsecolor', 'Scaling', 'joint', 'ColorChannels', [1 2 0]);
    imshow(AC1); title('A registered to C (rigid)');
    figure;
    AC2 = imfuse( Ac, C, 'falsecolor', 'Scaling', 'joint', 'ColorChannels', [1 2 0]);
    imshow(AC2);  title('A registered to C');
    figure;
    AC3 = imfuse( Ac2, C, 'falsecolor', 'Scaling', 'joint', 'ColorChannels', [1 2 0]);
    imshow(AC3);  title('A registered to C (nonrigid)');

[Adr,shifts_r,~,~,~] = normcorre(A,params.mcorr.normcorre_r,D);
[Ad,shifts,~,~,~] = normcorre(Adr,params.mcorr.normcorre_nr,D);
    figure;
    AD1 = imfuse( Adr, D, 'falsecolor', 'Scaling', 'joint', 'ColorChannels', [1 2 0]);
    imshow(AD1); title('A registered to D (rigid)');
    figure;
    AD2 = imfuse( Ad, D, 'falsecolor', 'Scaling', 'joint', 'ColorChannels', [1 2 0]);
    imshow(AD2);  title('A registered to D');
    
% B
[Bar,shifts_r,~,~,~] = normcorre(B,params.mcorr.normcorre_r,A);
[Ba,shifts,~,~,~] = normcorre(Bar,params.mcorr.normcorre_nr,A);
    figure;
    BA1 = imfuse( Bar, A, 'falsecolor', 'Scaling', 'joint', 'ColorChannels', [1 2 0]);
    imshow(BA1); title('B registered to A (rigid)');
    figure;
    BA2 = imfuse( Ba, A, 'falsecolor', 'Scaling', 'joint', 'ColorChannels', [1 2 0]);
    imshow(BA2); title('B registered to A');

[Bcr,shifts_r,~,~,~] = normcorre(B,params.mcorr.normcorre_r,C);
[Bc,shifts,~,~,~] = normcorre(Bcr,params.mcorr.normcorre_nr,C);
    figure;
    BC1 = imfuse( Bcr, C, 'falsecolor', 'Scaling', 'joint', 'ColorChannels', [1 2 0]);
    imshow(BC1); title('B registered to C (rigid)');
    figure;
    BC2 = imfuse( Bc, C, 'falsecolor', 'Scaling', 'joint', 'ColorChannels', [1 2 0]);
    imshow(BC2); title('B registered to C');

[Bdr,shifts_r,~,~,~] = normcorre(B,params.mcorr.normcorre_r,D);
[Bd,shifts,~,~,~] = normcorre(Bdr,params.mcorr.normcorre_nr,D);
    figure;
    BD1 = imfuse( Bdr, D, 'falsecolor', 'Scaling', 'joint', 'ColorChannels', [1 2 0]);
    imshow(BD1); title('B registered to D (rigid)');
    figure;
    BD2 = imfuse( Bd, D, 'falsecolor', 'Scaling', 'joint', 'ColorChannels', [1 2 0]);
    imshow(BD2); title('B registered to D');
    
% C
[Car,shifts_r,~,~,~] = normcorre(C,params.mcorr.normcorre_r,A);
[Ca,shifts,~,~,~] = normcorre(Car,params.mcorr.normcorre_nr,A);
    figure;
    CA1 = imfuse( Car, A, 'falsecolor', 'Scaling', 'joint', 'ColorChannels', [1 2 0]);
    imshow(CA1); title('C registered to A (rigid)');
    figure;
    CA2 = imfuse( Ca, A, 'falsecolor', 'Scaling', 'joint', 'ColorChannels', [1 2 0]);
    imshow(CA2); title('C registered to A');

[Cbr,shifts_r,~,~,~] = normcorre(C,params.mcorr.normcorre_r,B);
[Cb,shifts,~,~,~] = normcorre(Cbr,params.mcorr.normcorre_nr,B);
    figure;
    CB1 = imfuse( Cbr, B, 'falsecolor', 'Scaling', 'joint', 'ColorChannels', [1 2 0]);
    imshow(CB1); title('C registered to B (rigid)');
    figure;
    CB2 = imfuse( Cb, B, 'falsecolor', 'Scaling', 'joint', 'ColorChannels', [1 2 0]);
    imshow(CB2); title('C registered to B');

[Cdr,shifts_r,~,~,~] = normcorre(C,params.mcorr.normcorre_r,D);
[Cd,shifts,~,~,~] = normcorre(Cdr,params.mcorr.normcorre_nr,D);
    figure;
    CD1 = imfuse( Cdr, D, 'falsecolor', 'Scaling', 'joint', 'ColorChannels', [1 2 0]);
    imshow(CD1); title('C registered to D (rigid)');
    figure;
    CD2 = imfuse( Cd, D, 'falsecolor', 'Scaling', 'joint', 'ColorChannels', [1 2 0]);
    imshow(CD2); title('C registered to D');
    
% D
[Dar,shifts_r,~,~,~] = normcorre(D,params.mcorr.normcorre_r,A);
[Da,shifts,~,~,~] = normcorre(Dar,params.mcorr.normcorre_nr,A);
    figure;
    DA1 = imfuse( Dar, A, 'falsecolor', 'Scaling', 'joint', 'ColorChannels', [1 2 0]);
    imshow(DA1); title('D registered to A (rigid)');
    figure;
    DA2 = imfuse( Da, A, 'falsecolor', 'Scaling', 'joint', 'ColorChannels', [1 2 0]);
    imshow(DA2); title('D registered to A');

[Dbr,shifts_r,~,~,~] = normcorre(D,params.mcorr.normcorre_r,B);
[Db,shifts,~,~,~] = normcorre(Dbr,params.mcorr.normcorre_nr,B);
    figure;
    DB1 = imfuse( Dbr, B, 'falsecolor', 'Scaling', 'joint', 'ColorChannels', [1 2 0]);
    imshow(DB1); title('D registered to B (rigid)');
    figure;
    DB2 = imfuse( Db, B, 'falsecolor', 'Scaling', 'joint', 'ColorChannels', [1 2 0]);
    imshow(DB2); title('D registered to B');

[Dcr,shifts_r,~,~,~] = normcorre(D,params.mcorr.normcorre_r,C);
[Dc,shifts,~,~,~] = normcorre(Dcr,params.mcorr.normcorre_nr,C);
    figure;
    DC1 = imfuse( Dcr, C, 'falsecolor', 'Scaling', 'joint', 'ColorChannels', [1 2 0]);
    imshow(DC1); title('D registered to C (rigid)');
    figure;
    DC2 = imfuse( Dc, C, 'falsecolor', 'Scaling', 'joint', 'ColorChannels', [1 2 0]);
    imshow(DC2); title('D registered to C');