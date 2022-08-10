tic
try
    [data_locn,comp,err] = load_neuroSEEmodules;
    list = 'list_m62_fov1_fam1fam2-fam1.txt';
    reffile = '20181011_15_10_39';
    [ mouseid, expname, fov ] = find_mouseIDexpname(list);

    grp_sdir = [data_locn 'Analysis/' mouseid '/' fov '/' mouseid '_' expname ...
                '/group_proc/imreg_normcorre_CaImAn/'...
                mouseid '_' expname '_imreg_ref' reffile '/'];
    segmentout = load([grp_sdir mouseid '_' expname '_ref' reffile '_segment_output.mat']);

    bl_prctile = 83;
    params = neuroSEE_setparams(...
                'FOV',490,...
                'virus','GCaMP6s',...
                'bl_prctile',bl_prctile);
    
    df_f = segmentout.df_f;
    spikes = extractSpikes( df_f, params.spkExtract );        
    output.spikes = spikes;
    output.params = params.spkExtract;
%     for i = 1:184
%         df_f_up(i,:) = upsample( df_f(i,:),5 );
%     end
%     spikes_caiman_up = extractSpikes( df_f_up, params.spkExtract );
    
    save([grp_sdir 'noFISSA/bl_prctile' num2str(bl_prctile) '/' mouseid '_' expname '_ref' reffile '_spikes_caiman.mat'],'-struct','output')
    disp(toc)
catch
    toc
end

% later check
% for i = 1:10
%     figure;
%     subplot(511); plot(df_f(i,:)); ylabel('df/f');
%     subplot(512); plot(spikes_caiman(i,:)); ylabel('spikes');
% %     subplot(513); plot(df_f_up(i,:)); ylabel('df/f (up)');
% %     subplot(514); plot(spikes_caiman_up(i,:)); ylabel('spikes (up)');
% %     subplot(515); plot(downsample(spikes_caiman_up(i,:),5)); ylabel('spikes (up-down)');
% end