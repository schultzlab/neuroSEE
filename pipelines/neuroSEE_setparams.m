% Written by Ann Go, adapted from NoRMCorreSetParms.m
% Struct for setting the neuroSEE pipeline parameters. Any parameter that is
% not set gets a default value

function params = neuroSEE_setparams(varargin)
Names = [
    % general
        'dofissa            ' % flag to do fissa correction (default: true)
    % motion correction (general)
        'mcorr_method       ' % motion correction method (default: 'normcorre')
        'refChannel         ' % reference channel for motion correction (default: 'green')
    % motion correction: fftRigid    
        'imscale            ' % image downsampling factor (default: 1)
        'Nimg_ave           ' % no. of images to be averaged for calculating pixel shift (zippering) (default: 10)
        'redoT              ' % no. of frames at start of file to redo motion correction for after 1st iteration (default: 300)
    % motion correction: NoRmCorre
        % dataset info
        'd1                 ' % number of rows
        'd2                 ' % number of cols
        'd3                 ' % number of planes (for 3d imaging, default: 1)
        % patches
        'grid_size_r        ' % size of non-overlapping regions - rigid correction (default: [d1,d2,d3])
        'grid_size_nr       ' % size of non-overlapping regions - nonrigid correction (default: [d1/4,d2/4,d3])
        'max_shift_r        ' % maximum rigid shift in each direction - rigid correction (default: [30,30,1])
        'max_shift_nr       ' % maximum rigid shift in each direction - nonrigid correction (default: [20,20,1])
        'overlap_pre        ' % size of overlapping region (default: [64,64,1])
        'min_patch_size     ' % minimum size of patch (default: [64,64,1])    
        'min_diff           ' % minimum difference between patches (default: [32,32,1])
        'us_fac             ' % upsampling factor for subpixel registration (default: 50)
        'mot_uf_r           ' % degree of patches upsampling - rigid correction (default: [1,1,1])
        'mot_uf_nr          ' % degree of patches upsampling - nonrigid correction (default: [4,4,1])
        
        'max_dev            ' % maximum deviation of patch shift from rigid shift (default: [3,3,1])
        'overlap_post       ' % size of overlapping region after upsampling (default: [64,64,1])
        'phase_flag         ' % flag for using phase correlation (default: false)
        'shifts_method      ' % method to apply shifts ('FFT','cubic','linear')
        % template updating
        'upd_template       ' % flag for online template updating (default: true)
        'init_batch         ' % length of initial batch (default: 100)
        'bin_width          ' % width of each bin (default: 10)
        'buffer_width       ' % number of local means to keep in memory (default: 50)
        'method             ' % method for averaging the template (default: {'median';'mean})
        'iter               ' % number of data passes (default: 1)
        'boundary           ' % method of boundary treatment 'NaN','copy','zero','template' (default: 'copy')
        % misc
        'add_value          ' % add dc value to data (default: 0)
        'use_parallel       ' % for each frame, update patches in parallel (default: false)
        'memmap             ' % flag for saving memory mapped motion corrected file (default: false)
        'mem_filename       ' % name for memory mapped file (default: 'motion_corrected.mat')
        'mem_batch_size     ' % batch size during memory mapping for speed (default: 5000)
        % plotting
        'plot_flag          ' % flag for plotting results in real time (default: false)
        'make_avi           ' % flag for making movie (default: false)
        'name               ' % name for movie (default: 'motion_corrected.avi')
        'fr                 ' % frame rate for movie (default: 30.9)
        % output type
        'output_type        ' % 'mat' (load in memory), 'memmap', 'tiff', 'hdf5', 'bin' (default:mat)
        'h5_groupname       ' % name for hdf5 dataset (default: 'mov')
        'h5_filename        ' % name for hdf5 saved file (default: 'motion_corrected.h5')
        'tiff_filename      ' % name for saved tiff stack (default: 'motion_corrected.tif')
        'output_filename    ' % name for saved file will be used if `h5_,tiff_filename` are not specified
        % use windowing
        'use_windowing      ' % flag for windowing data before fft (default: false)
        'window_length      ' % length of window on each side of the signal as a fraction of signal length
                               %    total length = length(signal)(1 + 2*window_length). (default: 0.5)
        % bitsize for reading .raw files
        'bitsize            ' % (default: 2 (uint16). other choices 1 (uint8), 4 (single), 8 (double))
        % offset from bidirectional sampling
        'correct_bidir      ' % check for offset due to bidirectional scanning (default: true)
        'nFrames            ' % number of frames to average (default: 50)
        'bidir_us           ' % upsampling factor for bidirectional sampling (default: 10)
        'col_shift          ' % known bi-directional offset provided by the user (default: [])
        'print_msg          ' % flag to print messages (default: false)
    % ROI segmentation
        'df_prctile         ' % percentile to be used for estimating baseline (default: 5)
        'df_medfilt1        ' % degree of smoothing for df_f (default: 13)
        'cellrad_FOV490     ' % expected radius of a cell in 490x490 um FOV (pixels) (default: 6)  
        'maxcells_FOV490    ' % estimated number of cells in 490x490 um FOV (default: 300) 
        'cellrad_FOV330     ' % expected radius of a cell in 330x330 um FOV (pixels) (default: 9)  
        'maxcells_FOV330    ' % estimated number of cells in 330x330 um FOV (default: 200) 
    % neuropil correction
        'ddf_prctile        ' % percentile to be used for estimating baseline (default: 5)
        'ddf_medfilt1       ' % degree of smoothing for ddf_f (default: 17)
    % spike extraction
        'bl_prctile         ' % percentile to be used for estimating baseline (default: 85)
        'spk_SNR            ' % spike SNR for min spike value (default: 1)
        'decay_time         ' % length of a typical transient in seconds (default: 0.4)
        'lam_pr             ' % false positive probability for determing lambda penalty (default: 0.99)
    % PF mapping
        'Nepochs            ' % number of epochs into which data will be divided for analysis (default: 1)
        'histsmoothWin      ' % smoothing window for histogram method (default: 5)
        'Vthr               ' % speed threshold (mm/s) Note: David Dupret uses 20 mm/s (default: 20)
                              %                              Neurotar uses 8 mm/s
        'prctile_thr        ' % percentile threshold for filtering nonPCs (default: 99)
        'pfactivet_thr      ' % fraction of dwell time in place field cell is required to be active (default: 0)
        'Nrand              ' % number of shuffles for bootstrap test
        'Nbins_1D           ' % no. of position bins in 103-cm linear track (default: 50)
        'Nbins_2D           ' % position bins in 325-mm diameter open field arena (default: [16,16])
   ]; 
   
[m,~] = size(Names);
names = lower(Names);

% Combine all leading options structures o1, o2, ... in l1Set(o1,o2,...).
options = [];
for j = 1:m
    eval(['options.' Names(j,:) '= [];']);
end
i = 1;
while i <= nargin
    arg = varargin{i};
    if ischar(arg), break; end
    if ~isempty(arg)                      % [] is a valid options argument
        if ~isa(arg,'struct')
            error(sprintf(['Expected argument %d to be a string parameter name ' ...
                'or an options structure\ncreated with OPTIMSET.'], i));
        end
        for j = 1:m
            if any(strcmp(fieldnames(arg),deblank(Names(j,:))))
                eval(['val = arg.' Names(j,:) ';']);
            else
                val = [];
            end
            if ~isempty(val)
                eval(['options.' Names(j,:) '= val;']);
            end
        end
    end
    i = i + 1;
end

% A finite state machine to parse name-value pairs.
if rem(nargin-i+1,2) ~= 0
    error('Arguments must occur in name-value pairs.');
end
expectval = 0;                          % start expecting a name, not a value
while i <= nargin
    arg = varargin{i};
    
    if ~expectval
        if ~ischar(arg)
            error(sprintf('Expected argument %d to be a string parameter name.', i));
        end
        
        lowArg = lower(arg);
        j = strmatch(lowArg,names);
        if isempty(j)                       % if no matches
            error(sprintf('Unrecognized parameter name ''%s''.', arg));
        elseif length(j) > 1                % if more than one match
            % Check for any exact matches (in case any names are subsets of others)
            k = strmatch(lowArg,names,'exact');
            if length(k) == 1
                j = k;
            else
                msg = sprintf('Ambiguous parameter name ''%s'' ', arg);
                msg = [msg '(' deblank(Names(j(1),:))];
                for k = j(2:length(j))'
                    msg = [msg ', ' deblank(Names(k,:))];
                end
                msg = sprintf('%s).', msg);
                error(msg);
            end
        end
        expectval = 1;                      % we expect a value next
        
    else
        eval(['options.' Names(j,:) '= arg;']);
        expectval = 0;
        
    end
    i = i + 1;
end

if expectval
    error(sprintf('Expected value for parameter ''%s''.', arg));
end

Values = [
    % general
        {true}                % flag to do fissa correction (default: true)
    % motion correction (general)
        {'normcorre'}         % motion correction method (default: 'normcorre')
        {'green'}             % reference channel for motion correction (default: 'green')
    % motion correction: fftRigid    
        {1}                   % image downsampling factor (default: 1)
        {10}                  % no. of images to be averaged for calculating pixel shift (zippering) (default: 10)
        {300}                 % no. of frames at start of file to redo motion correction for after 1st iteration (default: 300)
    % NoRmCorre
        % dataset info
        {512}                 % number of rows
        {512}                 % number of columns
        {1}
        % patches
        {[512,512,1]}         % size of non-overlapping regions - rigid correction (default: [d1,d2,d3])
        {[128,128,1]}         % size of non-overlapping regions - nonrigid correction (default: [d1/4,d2/4,d3])
        {[30,30,1]}           % maximum rigid shift in each direction - rigid correction (default: [30,30,1])
        {[20,20,1]}           % maximum rigid shift in each direction - nonrigid correction (default: [20,20,1])
        {[64,64,1]}           % size of overlapping region (default: [32,32,1])
        {[64,64,1]}           % minimum size of patch (default: [32,32,1])    
        {[32,32,1]}           % minimum difference between patches (default: [16,16,1])
        {50}                  % upsampling factor for subpixel registration (default: 50)
        {[1,1,1]}             % degree of patches upsampling - rigid correction (default: [1,1,1])
        {[4,4,1]}             % degree of patches upsampling - nonrigid correction (default: [4,4,1])
        {[3,3,1]}             % maximum deviation of patch shift from rigid shift (default: [3,3,1])
        {[64,64,1]}           % size of overlapping region after upsampling (default: [32,32,1])
        {false}               % use phase correlation (good for high SNR)
        {'FFT'}               % method for applying shifts ('FFT', 'linear', 'cubic')
        % template updating
        {true}                % flag for online template updating (default: true)
        {100}                 % length of initial batch (default: 100)
        {50}                  % width of each bin (default: 10)
        {50}                  % number of local means to keep in memory (default: 50)
        {{'median';'mean'}}   % method for averaging the template (default: {'median';'mean'}
        {1}                   % number of data passes (default: 1)
        {'copy'}              % method of boundary treatment (default: 'copy')
        % misc
        {0}                   % add dc value to data (default: 0)
        {false}               % for each frame, update patches in parallel (default: false)
        {false}               % flag for saving memory mapped motion corrected file (default: false)
        {'motion_corrected.mat'} % name for memory mapped file (default: 'motion_corrected.mat')
        {1000}                % batch size used during memory mapping for faster mapping
        % plotting
        {false}               % flag for plotting results in real time (default: false)
        {false}               % flag for making movie (default: false)
        {'motion_corrected.avi'} % name for movie (default: 'motion_corrected.avi')
        {30.9}                % frame rate for movie (default: 30.9)   
        % output_type    
        {'mat'}
        {'mov'}
        {'motion_corrected.h5'}
        {'motion_corrected.tif'}
        {''}
        % use_windowing
        {false}
        {0.5}
        % bitsize for reading .raw files
        {2}
        % offset from bidirectional sampling
        {true}
        {50}
        {10}
        {[]}
        {false}               % flag to print messages (default: false)
    % ROI segmentation
        {5}                   % percentile to be used for estimating baseline (default: 5)
        {13}                  % degree of smoothing for df_f (default: 13)
        {6}                   % expected radius of a cell in 490x490 um FOV (pixels) (default: 6)  
        {400}                 % estimated number of cells in 490x490 um FOV (default: 300)
        {9}                   % expected radius of a cell in 330x330 um FOV (pixels) (default: 9)  
        {300}                 % estimated number of cells in 330x330 um FOV (default: 200)
    % neuropil correction
        {5}                   % percentile to be used for estimating baseline (default:5)
        {17}                  % degree of smoothing for ddf_f (default: 17)
    % spike extraction
        {85}                  % percentile to be used for estimating baseline (default: 85)
        {1}                   % spike SNR for min spike value (default: 1)
        {1}                   % length of a typical transient in seconds (default: 0.4)
        {0.99}                % false positive probability for determing lambda penalty (default: 0.99)
    % PF mapping
        {1}                   % number of epochs into which data will be divided for analysis (default: 1)
        {5}                   % smoothing window for histogram method (default: 5)
        {20}                  % speed threshold (mm/s) Note: David Dupret uses 20 mm/s (default: 20)
                              %                              Neurotar uses 8 mm/s
        {99}                  % percentile threshold for filtering nonPCs (default: 99)
        {0}                   % fraction of dwell time in place field cell is required to be active (default: 0)
        {1000}                % number of shuffles for bootstrap test
        {50}                  % no. of position bins in 103-cm linear track (default: 50)
        {[16,16]}             % position bins in 325-mm diameter open field arena (default: [16,16])
    ];

for j = 1:m
    if eval(['isempty(options.' Names(j,:) ')'])
        eval(['options.' Names(j,:) '= Values{j};']);
    end
end

if ~isempty(options.output_filename)
    out_type = options.output_type;
    [filepath,name,~] = fileparts(options.output_filename);
    output_filename = fullfile(filepath,name);
    if strcmpi(options.h5_filename,'motion_corrected.h5') && (strcmpi(out_type,'h5') || strcmpi(out_type,'hdf5'))
        options.h5_filename = [output_filename,'.h5'];
    end
    if strcmpi(options.tiff_filename,'motion_corrected.tif') && (strcmpi(out_type,'tif') || strcmpi(out_type,'tiff'))
        options.tiff_filename = [output_filename,'.tif'];
    end    
end

if isempty(options.d1); options.d1 = input('What is the total number of rows? \n'); end
if isempty(options.d2); options.d2 = input('What is the total number of columns? \n'); end
%if options.d3 == 1; nd = 2; else nd = 3; end
if isempty(options.grid_size_r); options.grid_size_r = [options.d1,options.d2,options.d3]; end
if isempty(options.grid_size_nr); options.grid_size_nr = [options.d1/4,options.d2/4,options.d3]; end
if length(options.grid_size_r) == 1; options.grid_size_r = options.grid_size_r*ones(1,3); end
if length(options.grid_size_r) == 2; options.grid_size_r(3) = 1; end
if length(options.grid_size_nr) == 1; options.grid_size_nr = options.grid_size_nr*ones(1,3); end
if length(options.grid_size_nr) == 2; options.grid_size_nr(3) = 1; end
if length(options.min_patch_size) == 1; options.min_patch_size = options.min_patch_size*ones(1,3); end
if length(options.min_patch_size) == 2; options.min_patch_size(3) = 1; end
if length(options.min_diff) == 1; options.min_diff = options.min_diff*ones(1,3); end
if length(options.min_diff) == 2; options.min_diff(3) = 1; end
if length(options.overlap_pre) == 1; options.overlap_pre = options.overlap_pre*ones(1,3); end
if length(options.overlap_pre) == 2; options.overlap_pre(3) = 1; end
if length(options.overlap_post) == 1; options.overlap_post = options.overlap_post*ones(1,3); end
if length(options.overlap_post) == 2; options.overlap_post(3) = 1; end
if length(options.max_shift_r) == 1; options.max_shift_r = options.max_shift_r*ones(1,3); end
if length(options.max_shift_r) == 2; options.max_shift_r(3) = 1; end
if length(options.max_shift_nr) == 1; options.max_shift_nr = options.max_shift_nr*ones(1,3); end
if length(options.max_shift_nr) == 2; options.max_shift_nr(3) = 1; end
if length(options.max_dev) == 1; options.max_dev = options.max_dev*ones(1,3); end
if length(options.max_dev) == 2; options.max_dev(3) = 1; end
if length(options.mot_uf_r) == 1; options.mot_uf_r = options.mot_uf_r*ones(1,3); end
if length(options.mot_uf_r) == 2; options.mot_uf_r(3) = 1; end
if length(options.mot_uf_nr) == 1; options.mot_uf_nr = options.mot_uf_nr*ones(1,3); end
if length(options.mot_uf_nr) == 2; options.mot_uf_nr(3) = 1; end
options.mot_uf_r(options.grid_size_r >= [options.d1,options.d2,options.d3]) = 1;
options.mot_uf_nr(options.grid_size_nr >= [options.d1,options.d2,options.d3]) = 1;

% output
params.fr = options.fr;
% motion correction
params.mcorr.refChannel = options.refChannel;
if strcmpi(options.mcorr_method,'fftRigid')
    f = {'imscale';'Nimg_ave';'redoT'};
    fn = fieldnames(options);
    for i = 1:length(fn)
        for j = 1:length(f)
            if strcmpi(fn{i},f{j})
                params.fftRigid.(fn{i}) = options.(fn{i});
            end
        end
    end
end
if strcmpi(options.mcorr_method,'normcorre') || strcmpi(options.mcorr_method,'normcorre-r')
    f = {'d1'; 'd2'; 'd3'; 'grid_size_r'; 'max_shift_r'; 'mot_uf_r';...
        'overlap_pre'; 'min_patch_size'; 'min_diff'; 'us_fac'; 'max_dev';'overlap_post';...
        'phase_flag'; 'shifts_method'; 'upd_template'; 'init_batch'; 'bin_width'; 'buffer_width';...
        'method'; 'iter'; 'boundary'; 'add_value'; 'use_parallel'; 'memmap';...
        'mem_filename'; 'mem_batch_size'; 'plot_flag'; 'make_avi'; 'name'; 'fr';...
        'output_type'; 'h5_groupname'; 'h5_filename'; 'tiff_filename'; 'output_filename';...
        'use_windowing'; 'window_length'; 'bitsize'; 'correct_bidir'; 'nFrames';...
        'bidir_us'; 'col_shift'; 'print_msg '};
    fn = fieldnames(options);
    for i = 1:length(fn)
        for j = 1:length(f)
            if strcmpi(fn{i},f{j})
                params.mcorr.normcorre_r.(fn{i}) = options.(fn{i});
            end
        end
    end
    params.mcorr.normcorre_r.grid_size = params.mcorr.normcorre_r.grid_size_r;
    params.mcorr.normcorre_r.max_shift = params.mcorr.normcorre_r.max_shift_r;
    params.mcorr.normcorre_r.mot_uf = params.mcorr.normcorre_r.mot_uf_r;
    params.mcorr.normcorre_r = rmfield(params.mcorr.normcorre_r,{'grid_size_r','max_shift_r','mot_uf_r'});
end
if strcmpi(options.mcorr_method,'normcorre') || strcmpi(options.mcorr_method,'normcorre-nr')
    f = {'d1'; 'd2'; 'd3'; 'grid_size_nr'; 'max_shift_nr'; 'mot_uf_nr';...
        'overlap_pre'; 'min_patch_size'; 'min_diff'; 'us_fac'; 'max_dev';'overlap_post';...
        'phase_flag'; 'shifts_method'; 'upd_template'; 'init_batch'; 'bin_width'; 'buffer_width';...
        'method'; 'iter'; 'boundary'; 'add_value'; 'use_parallel'; 'memmap';...
        'mem_filename'; 'mem_batch_size'; 'plot_flag'; 'make_avi'; 'name'; 'fr';...
        'output_type'; 'h5_groupname'; 'h5_filename'; 'tiff_filename'; 'output_filename';...
        'use_windowing'; 'window_length'; 'bitsize'; 'correct_bidir'; 'nFrames';...
        'bidir_us'; 'col_shift'; 'print_msg '};
    fn = fieldnames(options);
    for i = 1:length(fn)
        for j = 1:length(f)
            if strcmpi(fn{i},f{j})
                params.mcorr.normcorre_nr.(fn{i}) = options.(fn{i});
            end
        end
    end
    params.mcorr.normcorre_nr.grid_size = params.mcorr.normcorre_nr.grid_size_nr;
    params.mcorr.normcorre_nr.max_shift = params.mcorr.normcorre_nr.max_shift_nr;
    params.mcorr.normcorre_nr.mot_uf = params.mcorr.normcorre_nr.mot_uf_nr;
    params.mcorr.normcorre_nr = rmfield(params.mcorr.normcorre_nr,{'grid_size_nr','max_shift_nr','mot_uf_nr'});
end
% ROI segmentation
f = {'df_prctile'; 'df_medfilt1'; 'cellrad_FOV490'; 'maxcells_FOV490';...
     'cellrad_FOV330'; 'maxcells_FOV330'};
fn = fieldnames(options);
for i = 1:length(fn)
    for j = 1:length(f)
        if strcmpi(fn{i},f{j})
            params.ROIsegment.(fn{i}) = options.(fn{i});
        end
    end
end
% neuropil correction
if options.dofissa
    f = {'ddf_prctile'; 'ddf_medfilt1'};
    fn = fieldnames(options);
    for i = 1:length(fn)
        for j = 1:length(f)
            if strcmpi(fn{i},f{j})
                params.fissa.(fn{i}) = options.(fn{i});
            end
        end
    end
end
% spike extraction
f = {'bl_prctile'; 'spk_SNR'; 'decay_time'; 'lam_pr'};
fn = fieldnames(options);
for i = 1:length(fn)
    for j = 1:length(f)
        if strcmpi(fn{i},f{j})
            params.spkExtract.(fn{i}) = options.(fn{i});
        end
    end
end
% PF mapping
f = {'Nepochs'; 'histsmoothWin'; 'Vthr'; 'prctile_thr';...
     'Nlaps_thr'; 'Nbins_1D'; 'Nbins_2D'};
fn = fieldnames(options);
for i = 1:length(fn)
    for j = 1:length(f)
        if strcmpi(fn{i},f{j})
            params.PFmap.(fn{i}) = options.(fn{i});
        end
    end
end

