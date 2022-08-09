tic
try
    [data_locn,comp,err] = load_neuroSEEmodules;
    list = 'list_m62_fov1_fam1fam2-fam1.txt';
    reffile = '20181011_15_10_39';
    [ mouseid, expname, fov ] = find_mouseIDexpname(list);

    grp_sdir = [data_locn 'Analysis/' mouseid '/' fov '/' mouseid '_' expname ...
                '/group_proc/imreg_normcorre_CaImAn/'...
                mouseid '_' expname '_imreg_ref' reffile '/'];
    fissa = load([grp_sdir 'FISSA/' mouseid '_' expname '_ref' reffile '_fissa_output.mat']);

    params = neuroSEE_setparams(...
                'FOV',490,...
                'virus','GCaMP6s',...
                'bl_prctile',87);
    
    ddf_f = fissa.ddf_f;
    spikes_fissa = extractSpikes( ddf_f, params.spkExtract );        

    for i = 1:184
        ddf_f_up(i,:) = upsample( ddf_f(i,:),5 );
    end
    spikes_fissa_up = extractSpikes( ddf_f_up, params.spkExtract );
    
    save([grp_sdir 'FISSA/' mouseid '_' expname '_ref' reffile '_spikes_fissa.mat'],'spikes_fissa')
    disp(toc)
catch
    toc
end

% later check
% for i = 1:10
%     figure;
%     subplot(511); plot(ddf_f(i,:)); ylabel('df/f');
%     subplot(512); plot(spikes(i,:)); ylabel('spikes');
%     subplot(513); plot(spikes_fissa(i,:)); ylabel('spikes');
% %     subplot(513); plot(ddf_f_up(i,:)); ylabel('df/f (up)');
%     subplot(514); plot(spikes_fissa_up(i,:)); ylabel('spikes (up)');
% %     subplot(515); plot(downsample(spikes_fissa_up(i,:),5)); ylabel('spikes (up-down)');
% end