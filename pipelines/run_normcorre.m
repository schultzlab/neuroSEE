clear; close all;
tic

%% Basic setup
test = 0;       % set to 1 if testing, this will use one of smaller files in ../test
force = 0;      % set to 1 to overwrite saved processed files. This will 
                %    force pipeline to redo all steps incl. raw to tif
                %    conversion. If force = 0, processing step is skipped
                %    if output of said step already exists in individual
                %    file directory.

                
%% Load module folders and define data directory
data_locn = load_neuroSEEmodules(test);


%% Read files

file = '20181015_10_04_26';

%% Start processing

    
% Check if file has been processed. If not, continue processing unless forced to overwrite 
% existing processed data
fname_allData = fullfile(data_locn,'Data/',file(1:8),'/Processed/',file,'/',file,'_allData.mat');
if exist(fname_allData,'file')
    if ~force
        str = sprintf( '%s: File has been processed. Skipping processing\n', file );
        cprintf(str)
    end
end

if ~exist(fname_allData,'file') || force

    % Load raw and save tif files or load tif files if they exist
    [imG,imR] = load_imagefile( data_locn, file, force );

    % Do motion correction
    params.nonrigid = NoRMCorreSetParms(...
                        'd1',size(imG,1),...
                        'd2',size(imG,2),...
                        'grid_size',[32,32],...
                        'overlap_pre',[32,32],...
                        'overlap_post',[32,32],...
                        'iter',1,...
                        'use_parallel',false,...
                        'max_shift',30,...
                        'mot_uf',4,...
                        'bin_width',200,...
                        'max_dev',3,...
                        'us_fac',50,...
                        'init_batch',200);
    [imG, imR, shifts, template, options, col_shift, green, red] = normcorre_2ch( imG, imR, params.nonrigid);

    filedir = fullfile(data_locn,'Data/',file(1:8),'/Processed/',file,'/');
        if ~exist(filedir,'dir'), mkdir(filedir); end
        fname_mat_mcorr = [filedir file '_2P_mcorr_output.mat'];
        save(fname_mat_mcorr,'shifts','template','options','col_shift','green','red');

    fh = figure; 
    subplot(221), 
      imagesc( green.meanframe ); 
      axis image; colorbar; axis off;
      title( 'Mean frame for raw green' );
    subplot(222), 
      imagesc( green.meanregframe ); 
      axis image; colorbar; axis off; 
      title( 'Mean frame for corrected green' );
    subplot(223), 
      imagesc( red.meanframe ); 
      axis image; colorbar; axis off; 
      title( 'Mean frame for raw red' );
    subplot(224), 
      imagesc( red.meanregframe ); 
      axis image; colorbar; axis off;
      title( 'Mean frame for corrected red' );
    axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0  1],'Box','off',...
      'Visible','off','Units','normalized', 'clipping' , 'off');
%           titletext = [file(1:8) '-' file(10:11) '.' file(13:14) '.' file(16:17)];
%           text(0.5, 0.98,titletext);

        fname_fig = [filedir file '_2P_mcorr_summary.fig'];
        savefig( fh, fname_fig );
        saveas( fh, fname_fig(1:end-4), 'pdf' );
        close( fh );


    fname_tif_gr_mcorr = [filedir file '_green_mcorr.tif'];
    fname_tif_red_mcorr = [filedir file '_red_mcorr.tif'];

    prevstr = sprintf( '%s: Saving motion corrected tif images...\n', file );
            cprintf('Text',prevstr);
                writeTifStack( imG,fname_tif_gr_mcorr );
                writeTifStack( imR,fname_tif_red_mcorr );
            str = sprintf( '%s: Motion corrected tif images saved\n', file );
            refreshdisp( str, prevstr );
end

t = toc;
str = sprintf('Processing done in %g hrs',round(t/3600,2));
cprintf(str)

