%% USER INPUT
% 1) single file
file = '20181015_09_26_48';

% 2) single registered file
reffile = '20181015_09_37_54';                       % requires file above
mouseid = 'm62';                       % required if reffile is specified
expname = 'fam1nov-fam1';                       % required if reffile is specified

% 3) list
list = '';

% SETTINGS
groupreg_method = 'imreg';      % method for concatenating file data (either register images or rois)
imreg_method = 'normcorre';  % image registration method 
                                % values: [normcorre, normcorre-r, normcorre-nr, fftRigid] 
mcorr_method = 'normcorre';  % motion correction method used for reference file
                                % values: [normcorre, normcorre-r, normcorre-nr, fftRigid] 
                                    % CaImAn NoRMCorre method: 
                                    %   normcorre (rigid + nonrigid) 
                                    %   normcorre-r (rigid),
                                    %   normcorre-nr (nonrigid), 
                                    % fft-rigid method (Katie's)
segment_method = 'CaImAn';      % [ABLE,CaImAn]    
dofissa = true;                % flag if FISSA was implemented

%% file location
[data_locn,~,err] = load_neuroSEEmodules;
if ~isempty(err)
    beep
    cprintf('Errors',err);    
    return
end

% location of processed group data for list
if dofissa, str_fissa = 'FISSA'; else, str_fissa = 'noFISSA'; end
if ~isempty(list)
    [ mouseid, expname ] = find_mouseIDexpname(list);
    str_out = [ mouseid '_' expname ];

    if strcmpi(imreg_method, mcorr_method)
        savedir = [data_locn 'Analysis/' mouseid '/' mouseid '_' expname '/group_proc/'...
                    groupreg_method '_' imreg_method '_' segment_method '_' str_fissa '/'...
                    mouseid '_' expname '_imreg_ref' reffile '/'];
    else
        savedir = [data_locn 'Analysis/' mouseid '/' mouseid '_' expname '/group_proc/'...
                    groupreg_method '_' imreg_method '_' segment_method '_' str_fissa '/'...
                    mouseid '_' expname '_imreg_ref' reffile '_' mcorr_method '/'];
    end
    segment_outfile = [savedir mouseid '_' expname '_ref' reffile '_segment_output.mat'];
    if dofissa, fissa_outfile = [savedir mouseid '_' expname '_ref' reffile '_fissa_output.mat']; end
    spikes_outfile = [savedir mouseid '_' expname '_ref' reffile '_spikes.mat'];
    savename_pref = [savedir mouseid '_' expname '_ref' reffile];
else
    if ~isempty(reffile)
        str_out = [file '_ref' reffile];
        if strcmpi(imreg_method, mcorr_method)
            grpdir = [data_locn 'Analysis/' mouseid '/' mouseid '_' expname '/group_proc/'...
                        groupreg_method '_' imreg_method '_' segment_method '_' str_fissa '/'...
                        mouseid '_' expname '_imreg_ref' reffile '/'];
        else
            grpdir = [data_locn 'Analysis/' mouseid '/' mouseid '_' expname '/group_proc/'...
                        groupreg_method '_' imreg_method '_' segment_method '_' str_fissa '/'...
                        mouseid '_' expname '_imreg_ref' reffile '_' mcorr_method '/'];
        end
        segment_outfile = [grpdir mouseid '_' expname '_ref' reffile '_segment_output.mat'];
    
        if ~strcmpi(file,reffile)
            if strcmpi(imreg_method, mcorr_method)
                savedir = [data_locn 'Data/' file(1:8) '/Processed/' file '/imreg_' imreg_method '_ref' reffile '/'...
                    segment_method '_' mouseid '_' expname '/' str_fissa '/'];
            else
                savedir = [data_locn 'Data/' file(1:8) '/Processed/' file '/imreg_' imreg_method '_ref' reffile '_' mcorr_method '/'...
                    segment_method '_' mouseid '_' expname '/' str_fissa '/'];
            end
        else
            savedir = [data_locn 'Data/' file(1:8) '/Processed/' file '/mcorr_' mcorr_method '/'...
                    segment_method '_' mouseid '_' expname '/' str_fissa '/'];
        end
        fissa_outfile = [savedir file '_' mouseid '_' expname '_ref' reffile '_fissa_output.mat'];
        spikes_outfile = [savedir file '_' mouseid '_' expname '_ref' reffile '_spikes.mat'];
        savename_pref = [savedir file '_' mouseid '_' expname '_ref' reffile];
    else
        str_out = file;
        dir = [data_locn 'Data/' file(1:8) '/Processed/' file '/mcorr_' mcorr_method '/' segment_method '/'];
        segment_outfile = [ dir file '_segment_output.mat' ];
        savedir = [dir '/' str_fissa '/'];
        fissa_outfile = [ savedir  file '_spikes.mat' ];
        savename_pref = file;
    end
end

%% GUI
params = neuroSEE_setparams(...
            'groupreg_method', groupreg_method,...
            'imreg_method', imreg_method,...
            'mcorr_method', mcorr_method,...
            'segment_method', segment_method,...
            'dofissa', dofissa);
M = load(segment_outfile);
corr_image = M.corr_image;
masks = M.masks;
if isempty(list) && ~isempty(reffile)
    tsG = [];
    df_f = [];
else
    tsG = M.tsG;
    df_f = M.df_f;
end
if dofissa
    if exist(fissa_outfile,'file')
        M = load(fissa_outfile);
        dtsG = M.dtsG;
        ddf_f = M.ddf_f;
    else
        fprintf('%s: FISSA step was ticked but FISSA output does not exist. Cannot proceed.\n', str_out);
        return
    end
end
if exist(spikes_outfile,'file')
    M = load(spikes_outfile);
    spikes = M.spikes;
    params.spkExtract = M.params;
else
    fprintf('%s: spikes output does not exist. Cannot proceed.\n', str_out);
    return
end
GUI_manually_refine_spikes_imreg( spikes, tsG, dtsG, df_f, ddf_f, corr_image, masks, params, false, savedir, savename_pref )
