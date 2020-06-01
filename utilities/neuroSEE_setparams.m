% Written by Ann Go

function params = neuroSEE_setparams(mcorr_method, dofissa, default)

if nargin<2, default = true; end
    
if default
    params = load( 'default_params.mat' );
    
    % Remove irrelevant parameters 
    if strcmpi(mcorr_method,'normcorre')
        params.mcorr = rmfield(params.mcorr,'fftRigid');
    elseif strcmpi(mcorr_method,'normcorre-r')
        fields = {'normcorre_nr','fftRigid'};
        params.mcorr = rmfield(params.mcorr,fields);
    elseif strcmpi(mcorr_method,'normcorre-nr')
        fields = {'normcorre_r','fftRigid'};
        params.mcorr = rmfield(params.mcorr,fields);
    elseif strcmpi(mcorr_method,'fftRigid')
        fields = {'normcorre_r','normcorre_nr'};
        params.mcorr = rmfield(params.mcorr,fields);
    end
    
    if ~dofissa
        params = rmfield(params,'fissa');
    end

else
    params.fr = 30.9;                                   % imaging frame rate [default: 30.9]

    % motion correction
    params.mcorr.refChannel = 'green';                  % reference channel for motion correction [default: 'green']
    % Katie's method
    if strcmpi(mcorr_method,'fftRigid')
        params.mcorr.fftRigid.imscale = 1;              % image downsampling factor                                             [default: 1]
        params.mcorr.fftRigid.Nimg_ave = 10;            % no. of images to be averaged for calculating pixel shift (zippering)  [default: 10]
        params.mcorr.fftRigid.redoT = 300;              % no. of frames at start of file to redo motion correction for after 1st iteration [default: 300]
    end
    % NoRMCorre-rigid
    if or(strcmpi(mcorr_method,'normcorre'),strcmpi(mcorr_method,'normcorre-r'))
        params.mcorr.normcorre_r = NoRMCorreSetParms(...
            'd1', 512,...
            'd2', 512,...
            'max_shift',30,...          % default: 30
            'bin_width',200,...         % default: 200
            'us_fac',50,...             % default: 50
            'init_batch',200);          % default: 200
        params.mcorr.normcorre_r.print_msg = false;   % default: false
    end
    % NoRMCorre-nonrigid
    if or(strcmpi(mcorr_method,'normcorre'),strcmpi(mcorr_method,'normcorre-nr') )    
        params.mcorr.normcorre_nr = NoRMCorreSetParms(...
            'd1', 512,...
            'd2', 512,...
            'grid_size',[128,128],...   % default: [64,64]
            'overlap_pre',[64,64],...   % default: [64,64]
            'overlap_post',[64,64],...  % default: [64,64]
            'iter',1,...                % default: 1
            'use_parallel',false,...    % default: false
            'max_shift',20,...          % default: 20
            'mot_uf',4,...              % default: 4
            'bin_width',200,...         % default: 200
            'max_dev',3,...             % default: 3
            'us_fac',50,...             % default: 50
            'init_batch',200);          % default: 200
        params.mcorr.normcorre_nr.print_msg = false;    % default: false
    end

    % ROI segmentation
        params.ROIsegment.df_prctile = 5;     % percentile to be used for estimating baseline   [default: 5]
        params.ROIsegment.df_medfilt1 = 13;   % degree of smoothing for df_f                    [default: 23]
        params.ROIsegment.cellrad_FOV490 = 6;       % expected radius of a cell (pixels)    
        params.ROIsegment.maxcells_FOV490 = 300;    % estimated number of cells in FOV      
        params.ROIsegment.cellrad_FOV330 = 9;            
        params.ROIsegment.maxcells_FOV330 = 200;       

    % neuropil correction
    if dofissa
        params.fissa.ddf_prctile = 5;         % percentile to be used for estimating baseline   [default:5]
        params.fissa.ddf_medfilt1 = 17;       % degree of smoothing for ddf_f                   [default: 23]
    end

    % spike extraction
        params.spkExtract.bl_prctile = 85;    % percentile to be used for estimating baseline   [default:85]
        params.spkExtract.spk_SNR = 1;        % spike SNR for min spike value                   [default: 1]
        params.spkExtract.decay_time = 0.4;   % length of a typical transient in seconds        [default: 0.4]
        params.spkExtract.lam_pr = 0.99;      % false positive probability for determing lambda penalty   [default: 0.99]

    % PF mapping
        params.PFmap.Nepochs = 1;             % number of epochs for each 4 min video           [default: 1]
        params.PFmap.histsmoothWin = 3;       % smoothing window for histogram method           [default: 3]
        params.PFmap.Vthr = 20;               % speed threshold (mm/s) Note: David Dupret uses 20 mm/s    [default: 20]
                                              %                              Neurotar uses 8 mm/s
        params.PFmap.prctile_thr = 99;        % percentile threshold for filtering nonPCs       [default: 99]
        params.PFmap.Nlaps_thr = 0.5;         % fraction of laps that cell is required to be active to be 
                                              % considered for place cell analysis              [default: 0.5]
        params.PFmap.Nbins_1D = 50;           % no. of position bins in 103-cm linear track
        params.PFmap.Nbins_2D = [16 16];      % position bins in 325-mm diameter open field arena
end