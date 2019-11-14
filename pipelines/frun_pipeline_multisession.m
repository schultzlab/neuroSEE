

function frun_pipeline_multisession( list, force )

    if nargin<2, force = false; end

    %% Load module folders and define data directory
    tic
    [data_locn,~,err] = load_neuroSEEmodules(false);
    if ~isempty(err)
        beep 
        cprintf('Errors',err);    
        return
    end

    %% Find mouseid and experiment name
    tf = zeros(size(list));
    for ii = 1:numel(list)  
        tf(ii) = strcmpi(list(ii),'_');
    end
    ind = find(tf);
    mouseid = list(ind(1)+1:ind(2)-1);

    switch numel(ind)
        case 2
            expname = list(ind(2)+1:end-4);
        case 3
            expname = list(ind(2)+1:end-4);
        case 4
            expname = list(ind(2)+1:ind(3)-1);
        otherwise
            expname = list(ind(2)+1:ind(3)-1);
    end

    %% Load ROIs and templates from each session 
    listfile = [data_locn 'Digital_Logbook/lists/' list];
    files = extractFilenamesFromTxtfile(listfile);
    Nsessions = size(files,1);

    fdir = [ data_locn 'Analysis/' mouseid '/summaries based on registered ROIs/' mouseid '_' expname '_ref_' files(1,:) '/'];
        if ~exist(fdir,'dir'), mkdir(fdir); end
    fname = [fdir mouseid '_' expname '_ref_' files(1,:) '_registered_rois.mat' ];

    if ~exist(fname,'file') || force
        fprintf('%s: loading ROIs and templates\n',[mouseid '_' expname]);
        for jj = 1:Nsessions
            file = files(jj,:); 
            M = load([data_locn 'Data/' file(1:8) '/Processed/' file '/mcorr_normcorre/CaIman/' file '_segment_output.mat']);
            masks = M.masks;
            A_temp = zeros(size(masks,1)*size(masks,2),size(masks,3));
            for ii = 1:size(masks,3)
                masks_temp = masks(:,:,ii);
                A_temp(:,ii) = masks_temp(:);
                outlines{ii,jj} = bwboundaries(masks(:,:,ii));    % boundary of each ROI
            end
            A{jj} = sparse(A_temp); % masks

            t = load([data_locn 'Data/' file(1:8) '/Processed/' file '/mcorr_normcorre/' file '_mcorr_output.mat']);
            templates{jj} = t.template;
        end

        %% Register ROIs across sessions
        params.d1 = size(masks,1);
        params.d2 = size(masks,2);
        params.maxthr = [];
%         params.dist_maxthr = 0.1;
%         params.dist_exp = 1;
%         params.dist_thr = 0.5;
%         params.dist_overlap_thr = 0.8;
        params.dist_maxthr = 0.1;
        params.dist_exp = 0.5;
        params.dist_thr = 0.7;
        params.dist_overlap_thr = 0.5;


        params.plot_reg = true;
        params.print_msg = true;

        params_mc = NoRMCorreSetParms(...
                    'd1',params.d1,...        % width of image [default: 512]  *Regardless of user-inputted value, neuroSEE_motioncorrect reads this 
                    'd2',params.d2,...        % length of image [default: 512] *value from actual image    
                    'grid_size',[32,32],...     % default: [32,32]
                    'overlap_pre',[32,32],...   % default: [32,32]
                    'overlap_post',[32,32],...  % default: [32,32]
                    'iter',1,...                % default: 1
                    'use_parallel',false,...    % default: false
                    'max_shift',50,...          % default: 50
                    'mot_uf',4,...              % default: 4
                    'bin_width',200,...         % default: 200
                    'max_dev',3,...             % default: 3
                    'us_fac',50,...             % default: 50
                    'init_batch',200);          % default: 200
        params_mc.print_msg = false;

        fprintf('%s: registering ROIs\n',[mouseid '_' expname]);
        [A_union, assignments, matchings] = register_multisession(A, params, templates, params_mc);
        masks = reshape(full(A_union), params.d1, params.d2, size(A_union,2));
        Nrois = size(masks,3);

        % save masks
        registered_rois.masks = masks;
        registered_rois.outlines = outlines;
        registered_rois.assignments = assignments;
        registered_rois.matchings = matchings;
        registered_rois.params = params;
        registered_rois.params_mc = params_mc;

        % fprintf('%s: saving registered ROIs\n',[mouseid '_' expname]);
        % save(fname, '-struct', 'registered_rois')
            
        % superposition plots of matched ROIs
        [nRow, nCol] = getnRownCol(Nrois);
        nPlot = nRow*nCol;

        if ~exist([fdir 'registered_rois/'],'dir')
            mkdir([fdir 'registered_rois/']);
        end
        for n = 0:Nrois/nPlot
            fh = figure;
            ha = tight_subplot(nRow,nCol,[.01 .01],[.01 .05],[.01 .01]);
            for ii = 0:nPlot-1
                if (n*nPlot+ii+1) <= Nrois
                    axes(ha(ii+1));
                    imshow(zeros(512,512)); hold on
                    for jj = 1:Nsessions
                        k = assignments(n*nPlot+ii+1,jj);
                        if ~isnan(k)
                            plot(outlines{k,jj}{1}(:,2),outlines{k,jj}{1}(:,1),'Linewidth',1.5); hold on
                        end
                    end
                    axis off; title(['Cell ' num2str(n*nPlot+jj+1)],'fontsize',15);
                end
            end
%             if Nrois/nPlot <= 1
%                 fname_fig = [fdir 'registered_rois/' mouseid '_' expname '_ref_' files(1,:) '_registered_rois'];
%             else
%                 fname_fig = [fdir 'registered_rois/' mouseid '_' expname '_ref_' files(1,:) '_registered_rois_' num2str(n+1)];
%             end
%             savefig( fh, fname_fig );
%             saveas( fh, fname_fig, 'png' );
%             close( fh );
        end 

    else
        fprintf('%s: loading registered ROIs\n',[mouseid '_' expname]);
        load(fname)
        Nrois = size(masks,3);
    end
    
    %% Collate the data for tsG, dtsG, df_f, ddf_f, spikes for each roi in A_union
    fname = [fdir mouseid '_' expname '_ref_' files(1,:) '_spikes_tracking_data.mat' ];

    if ~exist(fname,'file') || force
        fprintf('%s: concatinating timeseries and tracking data\n',[mouseid '_' expname]);
        downr_all = [];
        for jj = 1:Nsessions
            file = files(jj,:);
            ts{jj} = load([data_locn 'Data/' file(1:8) '/Processed/' file '/mcorr_normcorre/CaIman/FISSA/' file '_spikes_output.mat']);
            Nt(jj) = size(ts{jj}.spikes,2);

            trackfile = dir([data_locn 'Data/' file(1:8) '/Processed/' file '/behaviour/*.mat']);
            td = load([data_locn 'Data/' file(1:8) '/Processed/' file '/behaviour/' trackfile.name]);
            x = td.x;
            y = td.y;
            r = td.r;
            phi = td.phi;
            speed = td.speed;
            tracktime = td.time;

            % Pre-process tracking data
            t0 = tracktime(1);                  % initial time in tracking data

            % Convert -180:180 to 0:360
            if min(phi)<0
               phi(phi<0) = phi(phi<0)+360;
            end

            % generate imaging timestamps using known image frame rate
            fr = 30.9;
            dt = 1/fr;
            t = (t0:dt:Nt(jj)*dt)';
            if length(t) ~= Nt(jj)
                t = (t0:dt:(Nt(jj)+1)*dt)';
            end

            % Downsample tracking to Ca trace
            downData{jj}.phi = interp1(tracktime,phi,t,'linear');
            downData{jj}.x = interp1(tracktime,x,t,'linear');
            downData{jj}.y = interp1(tracktime,y,t,'linear');
            downData{jj}.speed = interp1(tracktime,speed,t,'linear'); % mm/s
            downData{jj}.r = interp1(tracktime,r,t,'linear'); % mm/s
            downData{jj}.time = t;
            r_all = [downr_all; downData{jj}.r];
        end
        clear phi x y speed r t

        % Initialise cell arrays
        tsG = cell(Nrois,1);        dtsG = cell(Nrois,1); 
        df_f = cell(Nrois,1);       ddf_f = cell(Nrois,1); 
        spikes = cell(Nrois,1);
        trackData.x = cell(Nrois,1);      trackData.y = cell(Nrois,1);
        trackData.phi = cell(Nrois,1);    trackData.r = cell(Nrois,1); 
        trackData.time = cell(Nrois,1);   trackData.speed = cell(Nrois,1);

        for ii = 1:Nrois
            for jj = 1:Nsessions
                k = assignments(ii,jj);
                if ~isnan(k)
                    tsG{ii} = [tsG{ii}; ts{jj}.tsG(k,:)'];
                    dtsG{ii} = [dtsG{ii}; ts{jj}.dtsG(k,:)'];
                    df_f{ii} = [df_f{ii}; ts{jj}.df_f(k,:)'];
                    ddf_f{ii} = [ddf_f{ii}; ts{jj}.ddf_f(k,:)'];
                    spikes{ii} = [spikes{ii}; ts{jj}.spikes(k,:)'];

                    trackData.phi{ii} = [trackData.phi{ii}; downData{jj}.phi];
                    trackData.x{ii} = [trackData.x{ii}; downData{jj}.x];
                    trackData.y{ii} = [trackData.y{ii}; downData{jj}.y];
                    trackData.speed{ii} = [trackData.speed{ii}; downData{jj}.speed];
                    trackData.r{ii} = [trackData.r{ii}; downData{jj}.r];
                    trackData.r_all = r_all;
                    if ~isempty(trackData.time{ii})
                        tt = trackData.time{ii};
                    else
                        tt = 0;
                    end
                    trackData.time{ii} = [trackData.time{ii}; tt(end) + downData{jj}.time];
                end
            end
        end
                
        % save timeseries and track data
        spike_track_data.tsG = tsG;
        spike_track_data.dtsG = dtsG;
        spike_track_data.df_f = df_f;
        spike_track_data.ddf_f = ddf_f;
        spike_track_data.spikes = spikes;

        spike_track_data.trackData = trackData;

        % save(fname, '-struct', 'spike_track_data')
    else
        fprintf('%s: loading timeseries and tracking data\n',[mouseid '_' expname]);
        load(fname)
    end

    %% Calculate place field maps

    % settings
    if any(trackData.r_all < 100)
        params.mode_dim = '2D';         % open field
        params.PFmap.Nbins = [16, 16];  % number of location bins in [x y]               
    else 
        params.mode_dim = '1D';         % circular linear track
        params.PFmap.Nbins = 30;        % number of location bins               
    end
    params.PFmap.Nepochs = 1;
    params.PFmap.Vthr = 20;
    params.PFmap.histsmoothFac = 7;
    params.fr = 30.9;

    Nepochs = params.PFmap.Nepochs;
    fname = [fdir mouseid '_' expname '_ref_' files(1,:) '_PFmap_output.mat'];
    if ~exist(fname,'file') || force
        fprintf('%s: generating PFmaps\n', [mouseid '_' expname]);
        if strcmpi(params.mode_dim,'1D')
            % Generate place field maps
            [ occMap, hist, asd, activeData ] = generatePFmap_1d_multisession( spikes, trackData, params );

            % If 1D, sort place field maps 
            [ hist.sort_pfMap, hist.sortIdx ] = sortPFmap_1d( hist.pfMap, hist.infoMap, Nepochs );
            [ asd.sort_pfMap, asd.sortIdx ] = sortPFmap_1d( asd.pfMap, asd.infoMap, Nepochs );
            for en = 1:Nepochs
                hist.sort_pfMap_sm(:,:,en) = hist.pfMap_sm(hist.sortIdx(:,en),:,en);
                hist.sort_normpfMap(:,:,en) = hist.normpfMap(hist.sortIdx(:,en),:,en);
                hist.sort_normpfMap_sm(:,:,en) = hist.normpfMap_sm(hist.sortIdx(:,en),:,en);
                if numel(asd.pcIdx) ~=0
                    asd.sort_pfMap(:,:,en) = asd.pfMap(asd.sortIdx(:,en),:,en);
                    asd.sort_normpfMap(:,:,en) = asd.normpfMap(asd.sortIdx(:,en),:,en);
                else
                    asd.sort_pfMap(:,:,en) = zeros(size(asd.pfMap(:,:,en)));
                    asd.sort_normpfMap(:,:,en) = zeros(size(asd.pfMap(:,:,en)));
                end
            end

            % Make plots
            plotPF_1d(occMap, hist, asd);

            % Save output
            output.occMap = occMap;
            output.hist = hist;
            output.asd = asd;
            output.activeData = activeData;
            output.params = params.PFmap;
            % save(fname,'-struct','output');
        else % '2D'
            [occMap, spkMap, spkIdx, hist, asd, ~, activeData] = generatePFmap_2d(spikes, [], trackData, params, false);

             % Make plots
            plotPF_2d(spkMap, activeData, hist, asd);

            % Save output
            output.occMap = occMap;
            output.spkMap = spkMap;
            output.spkIdx = spkIdx;
            output.hist = hist;
            output.asd = asd;
            output.activeData = activeData;
            output.params = params.PFmap;
            save(fname,'-struct','output');
        end
    end

    t = toc;
    str = sprintf('%s: Processing done in %g hrs\n', [mouseid '_' expname], round(t/3600,2));
    cprintf(str)
end

function plotPF_1d(occMap, hist, asd)
    Npcs = length(hist.pcIdx);
    Npcs_asd = length(asd.pcIdx);

    if ~exist([fdir 'PFmaps/'],'dir'), mkdir([fdir 'PFmaps/']); end
        
    % summary of occMap, spkMaps, pfMaps
    for e = 1:Nepochs
        fh = figure('Position',[1087 648 800 800]);
        subplot(10,8,2:4); imagesc(occMap(e,:));
            xticks([]); yticks([]); ylabel('Occ');
            title('Histogram estimation'); colorbar;
        
        if Npcs > 0     
            subplot(10,8,[10:12,18:20,26:28]);
                imagesc(hist.spkMap(hist.sortIdx(:,e),:,e));
                xticks([]);
                yticks([1 Npcs]); ylabel('Cell #'); 
                title('Spike map'); colorbar;
            subplot(10,8,[33,41,49]); 
                imagesc(hist.infoMap(hist.sortIdx,1,e));
                xticks([]);
                yticks([1 Npcs]); ylabel('Cell #'); 
                title('Max MI'); colorbar;
            subplot(10,8,[34:36,42:44,50:52]);    
                imagesc(hist.sort_pfMap(:,:,e)); 
                xticks([]); yticks([1 Npcs]);
                title('Place field map'); colorbar;
            subplot(10,8,[58:60,66:68,74:76]);    
                imagesc(hist.sort_pfMap_sm(:,:,e)); 
                yticks([1 Npcs]); ylabel('Cell #');
                xticks([1 15 30]); xticklabels([1 50 100]);
                xlabel('Position (cm)');
                title('Smoothened pf map'); colorbar;
        end
            
        if Npcs_asd > 0    
            subplot(10,8,6:8); imagesc(occMap(e,:));
                xticks([]); yticks([]); ylabel('Occ');
                title('ASD'); colorbar;
            subplot(10,8,[14:16,22:24,30:32]);
                imagesc(asd.spkMap(asd.sortIdx(:,e),:,e));
                xticks([]);
                yticks([1 Npcs_asd]);  ylabel('Cell #'); 
                title('Spike map'); colorbar;
            subplot(10,8,[37,45,53]); imagesc(asd.infoMap(asd.sortIdx,1,e));
                xticks([]);
                yticks([1 Npcs_asd]); 
                title('Max MI'); colorbar;
            subplot(10,8,[38:40,46:48,54:56]);    
                imagesc(asd.sort_pfMap(:,:,e)); 
                yticks([1 Npcs_asd]);
                xticks([1 15 30]); xticklabels([1 50 100]);
                xlabel('Position (cm)');
                title('Place field map'); colorbar;
        end

        if Nepochs == 1
            fname_fig = [fdir 'PFmaps/' mouseid '_' expname '_ref_' files(1,:)  '_PFmaps'];
        else
            fname_fig = [fdir 'PFmaps/' mouseid '_' expname '_ref_' files(1,:)  '_PFmaps_' num2str(e) 'of' num2str(Nepochs) 'ep'];
        end
        savefig( fh, fname_fig );
        saveas( fh, fname_fig, 'png' );
        close( fh );
    end

    % summary of occMap, spkMaps, normpfMaps
    for e = 1:Nepochs
        fh = figure('Position',[1087 648 800 800]);
        subplot(10,8,2:4); imagesc(occMap(e,:));
            xticks([]); yticks([]); ylabel('Occ');
            title('Histogram estimation'); colorbar;
            
        if Npcs > 0
            subplot(10,8,[10:12,18:20,26:28]);
                imagesc(hist.normspkMap(hist.sortIdx(:,e),:,e));
                xticks([]);
                yticks([1 Npcs]); ylabel('Cell #'); 
                title('Normalised spk map'); colorbar;
            subplot(10,8,[33,41,49]); 
                imagesc(hist.infoMap(hist.sortIdx,1,e));
                xticks([]);
                yticks([1 Npcs]); ylabel('Cell #'); 
                title('Max MI'); colorbar;
            subplot(10,8,[34:36,42:44,50:52]);    
                imagesc(hist.sort_normpfMap(:,:,e)); 
                xticks([]); yticks([1 Npcs]);
                title('Normalised pf map'); colorbar;
            subplot(10,8,[58:60,66:68,74:76]);    
                imagesc(hist.sort_normpfMap_sm(:,:,e)); 
                yticks([1 Npcs]); ylabel('Cell #');
                xticks([1 15 30]); xticklabels([1 50 100]);
                xlabel('Position (cm)');
                title('Norm smooth pf map'); colorbar;
        end
        
        if Npcs_asd > 0
            subplot(10,8,6:8); imagesc(occMap(e,:));
                xticks([]); yticks([]); ylabel('Occ');
                title('ASD'); colorbar;
            subplot(10,8,[14:16,22:24,30:32]);
                imagesc(asd.normspkMap(asd.sortIdx(:,e),:,e));
                xticks([]);
                yticks([1 Npcs_asd]);  ylabel('Cell #'); 
                title('Normalised spk map'); colorbar;
            subplot(10,8,[37,45,53]); 
                imagesc(asd.infoMap(asd.sortIdx,1,e));
                xticks([]);
                yticks([1 Npcs_asd]); 
                title('Max MI'); colorbar;
            subplot(10,8,[38:40,46:48,54:56]);    
                imagesc(asd.sort_normpfMap(:,:,e)); 
                yticks([1 Npcs_asd]);
                xticks([1 15 30]); xticklabels([1 50 100]);
                xlabel('Position (cm)');
                title('Normalised pf map'); colorbar;
        end

        if Nepochs == 1
            fname_fig = [fdir 'PFmaps/' mouseid '_' expname '_ref_' files(1,:) '_normPFmaps'];
        else
            fname_fig = [fdir 'PFmaps/' mouseid '_' expname '_ref_' files(1,:) '_normPFmaps_' num2str(e) 'of' num2str(Nepochs) 'ep'];
        end
        savefig( fh, fname_fig );
        saveas( fh, fname_fig, 'png' );
        close( fh );
    end

    % per trial spike maps for all cells
    [nRow, nCol] = getnRownCol(Nrois);
    nPlot = nRow*nCol;
    
    for ii=0:Nrois/nPlot
        fh = figure;
        ha = tight_subplot(nRow,nCol,[.01 .01],[.01 .05],[.01 .01]);
        for jj=0:nPlot-1
            if (ii*nPlot+jj+1) <= Nrois
                Ntrials = size(activeData.normspkMap_pertrial{ii*nPlot+jj+1},1);
                axes(ha(+jj+1));
                imagesc(activeData.normspkMap_pertrial{ii*nPlot+jj+1}); 
                yticks(1:Ntrials:Ntrials); yticklabels([Ntrials,1]); % ylabel('Trial #');
                xticks([]); % xticklabels([1 50 100]); xlabel('Position (cm)');
                title(['Cell ' num2str(ii*nPlot+jj+1)],'fontsize',15);
            end
        end
%         if Npcs/nPlot <= 1
%             fname_fig = [fdir 'PFmaps/' mouseid '_' expname '_ref_' files(1,:) '_spk_pertrial_hist'];
%         else
%             fname_fig = [fdir 'PFmaps/' mouseid '_' expname '_ref_' files(1,:) '_spk_pertrial_hist_' num2str(ii+1)];
%         end
%         savefig( fh, fname_fig );
%         saveas( fh, fname_fig, 'png' );
%         close( fh );
    end 


    % per trial spike maps for place cells
    [nRow, nCol] = getnRownCol(Npcs);
    nPlot = nRow*nCol;
    
    % histogram
    for ii=0:Npcs/nPlot
        fh = figure;
        ha = tight_subplot(nRow,nCol,[.01 .01],[.01 .05],[.01 .01]);
        for jj=0:nPlot-1
            if (ii*nPlot+jj+1) <= Npcs
                axes(ha(+jj+1));
                imagesc(hist.spkMap_pertrial{ii*nPlot+jj+1}); 
                yticks(1:Ntrials:Ntrials); yticklabels([Ntrials,1]); % ylabel('Trial #');
                xticks([]); % xticklabels([1 50 100]); xlabel('Position (cm)');
                title(['PC_h_i_s_t ' num2str(ii*nPlot+jj+1)],'fontsize',15);
            end
        end
%         if Npcs/nPlot <= 1
%             fname_fig = [fdir 'PFmaps/' mouseid '_' expname '_ref_' files(1,:) '_spk_pertrial_hist'];
%         else
%             fname_fig = [fdir 'PFmaps/' mouseid '_' expname '_ref_' files(1,:) '_spk_pertrial_hist_' num2str(ii+1)];
%         end
%         savefig( fh, fname_fig );
%         saveas( fh, fname_fig, 'png' );
%         close( fh );
    end 

    for ii=0:Npcs/nPlot
        fh = figure;
        ha = tight_subplot(nRow,nCol,[.01 .01],[.01 .05],[.01 .01]);
        for jj=0:nPlot-1
            if (ii*nPlot+jj+1) <= Npcs
                axes(ha(+jj+1));
                imagesc(hist.normspkMap_pertrial{ii*nPlot+jj+1}); 
                yticks(1:Ntrials:Ntrials); yticklabels([Ntrials,1]); %ylabel('Trial #');
                xticks([]); % xticklabels([1 50 100]); xlabel('Position (cm)');
                title(['PC_h_i_s_t ' num2str(ii*nPlot+jj+1)],'fontsize',15);
            end
        end
%         if Npcs/nPlot <= 1
%             fname_fig = [fdir 'PFmaps/' mouseid '_' expname '_ref_' files(1,:)e '_normspk_pertrial_hist'];
%         else
%             fname_fig = [fdir 'PFmaps/' mouseid '_' expname '_ref_' files(1,:) '_normspk_pertrial_hist_' num2str(ii+1)];
%         end
%         savefig( fh, fname_fig );
%         saveas( fh, fname_fig, 'png' );
%         close( fh );
    end 

    % asd
    [nRow, nCol] = getnRownCol(Npcs_asd);
    nPlot = nRow*nCol;

    Ntrials = size(asd.spkMap_pertrial,1);
    for ii=0:Npcs_asd/nPlot
        fh = figure;
        ha = tight_subplot(nRow,nCol,[.01 .01],[.01 .05],[.01 .01]);
        for jj=0:nPlot-1
            if (ii*nPlot+jj+1) <= Npcs_asd
                axes(ha(+jj+1));
                imagesc(asd.spkMap_pertrial{ii*nPlot+jj+1}); 
                yticks(1:Ntrials:Ntrials); yticklabels([Ntrials,1]); ylabel('Trial #');
                xticks([]); % xticklabels([1 50 100]); xlabel('Position (cm)');
                title(['PC_a_s_d ' num2str(ii*nPlot+jj+1)],'fontsize',15);
            end
        end
        if Npcs_asd/nPlot <= 1
            fname_fig = [fdir 'PFmaps/' mouseid '_' expname '_ref_' files(1,:) '_spk_pertrial_asd'];
        else
            fname_fig = [fdir 'PFmaps/' mouseid '_' expname '_ref_' files(1,:) '_spk_pertrial_asd_' num2str(ii+1)];
        end
        savefig( fh, fname_fig );
        saveas( fh, fname_fig, 'png' );
        close( fh );
    end 

    for ii=0:Npcs_asd/nPlot
        fh = figure;
        ha = tight_subplot(nRow,nCol,[.01 .01],[.01 .05],[.01 .01]);
        for jj=0:nPlot-1
            if (ii*nPlot+jj+1) <= Npcs_asd
                axes(ha(+jj+1));
                imagesc(asd.normspkMap_pertrial{ii*nPlot+jj+1}); 
                yticks(1:Ntrials:Ntrials); yticklabels([Ntrials,1]); ylabel('Trial #');
                xticks([]); % xticklabels([1 50 100]); xlabel('Position (cm)');
                title(['PC_a_s_d ' num2str(ii*nPlot+jj+1)],'fontsize',15);
            end
        end
        if Npcs_asd/nPlot <= 1
            fname_fig = [fdir 'PFmaps/' mouseid '_' expname '_ref_' files(1,:) '_normspk_pertrial_asd'];
        else
            fname_fig = [fdir 'PFmaps/' mouseid '_' expname '_ref_' files(1,:) '_normspk_pertrial_asd_' num2str(ii+1)];
        end
        savefig( fh, fname_fig );
        saveas( fh, fname_fig, 'png' );
        close( fh );
    end 

    % remapping within a session
    if Nepochs > 1
        fh = figure;
        for ei = 1:Nepochs % rows: sorting
            for ej = 1:Nepochs % cols: epochs 
                subplot(Nepochs, Nepochs, (ei-1)*Nepochs + ej); imagesc(hist.normpfMap(hist.sortIdx(:,ei),:,ej)); 
                title(['Epoch ' num2str(ej)]); ylabel(['Epoch' num2str(ei) ' sorting']);
            end
        end
        fname_fig = [fdir 'PFmaps/' mouseid '_' expname '_ref_' files(1,:) '_remapping_hist'];
        savefig( fh, fname_fig );
        saveas( fh, fname_fig, 'png' );
        close( fh );

        fh = figure;
        for ei = 1:Nepochs % rows: sorting
            for ej = 1:Nepochs % cols: epochs 
                subplot(Nepochs, Nepochs, (ei-1)*Nepochs + ej); imagesc(asd.normpfMap(asd.sortIdx(:,ei),:,ej)); 
                title(['Epoch ' num2str(ej)]); ylabel(['Epoch' num2str(ei) ' sorting']);
            end
        end
        fname_fig = [fdir 'PFmaps/' mouseid '_' expname '_ref_' files(1,:) '_remapping_asd'];
        savefig( fh, fname_fig );
        saveas( fh, fname_fig, 'png' );
        close( fh );
    end
end

function plotPF_2d(spkMap, activeData, hist, asd)
    Nspk = size(spkMap,3);
    nPlot = 4;
    for e = 1:Nepochs
        for ii=0:(Nspk/nPlot)-1 
            fh = figure; 
            ha = tight_subplot(nPlot,4,[.01 .005],[.01 .07],[.01 .01]);
            for jj=0:3
                if (ii*nPlot+jj) <= Nspk
                    axes(ha(jj*nPlot+1));
                    z = activeData.spikes(spkIdx(ii*nPlot+jj+1),:);
                    hold on; axis off;
                    plot(activeData.x,-activeData.y); plot(activeData.x(z>0),-activeData.y(z>0),'r.','markersize',10);
                    title(['Cell ',num2str(ii*nPlot+jj+1)],'fontsize',15);
                    axes(ha(jj*nPlot+2));
                    imagesc(squeeze(hist.pfMap(:,:,ii*nPlot+jj+1,e))');
                    axis off; colorbar; % caxis([0 0.06]);
                    if Nepochs >1 
                        title(['Epoch ',num2str(e)],'fontsize',15);
                    end
                    axes(ha(jj*nPlot+3)); 
                    imagesc(squeeze(hist.pfMap_sm(:,:,ii*nPlot+jj+1,e))');
                    axis off; colorbar; % caxis([0 0.005]);
                    axes(ha(jj*nPlot+4));
                    imagesc(squeeze(asd.pfMap(:,:,ii*nPlot+jj+1,e))');
                    axis off; colorbar; % caxis([0 0.003]);
                end
            end
            
            if ~exist([fdir 'PFmaps/'],'dir'), mkdir([fdir 'PFmaps/']); end
            if Nspk/nPlot <= 1
                if Nepochs == 1
                    fname_fig = [fdir 'PFmaps/' mouseid '_' expname '_ref_' files(1,:) '_PFmaps'];
                else
                    fname_fig = [fdir 'PFmaps/' mouseid '_' expname '_ref_' files(1,:) '_PFmaps_' num2str(e) 'of' num2str(Nepochs) 'ep' ];
                end
            else
                if Nepochs == 1
                    fname_fig = [fdir 'PFmaps/' mouseid '_' expname '_ref_' files(1,:) '_PFmaps_' num2str(ii+1)];
                else
                    fname_fig = [fdir 'PFmaps/' mouseid '_' expname '_ref_' files(1,:) '_PFmaps_' num2str(ii+1) '_' num2str(e) 'of' num2str(Nepochs) 'ep' ];
                end
            end
            savefig( fh, fname_fig );
            saveas( fh, fname_fig, 'png' );
            close( fh );
        end 
    end
end