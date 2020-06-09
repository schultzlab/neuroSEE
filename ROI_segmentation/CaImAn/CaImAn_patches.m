% Adapted by Ann Go from CaImAn's run_pipeline.m
%
% This function extracts ROIs and their decontaminated signals from a green
% channel image stack using CaImAn

function [C_df, masks, corr_image, F0, GUIdata] = CaImAn_patches( imG, options, display )

if nargin<3, display = false; end
if nargin<2 
    params = neuroSEE_setparams();
    options = params.ROIsegment.CaImAn;
end

fprintf('\tExtracting ROIs with CaImAn_patches\n');

sizY = size(imG);

%% Set parameters
options.d1 = sizY(1);
options.d2 = sizY(2);

patch_size = options.patch_size;          % size of each patch along each dimension (optional, default: [32,32])
overlap = options.overlap;                % amount of overlap in each dimension (optional, default: [4,4])
patches = construct_patches(sizY(1:end-1),patch_size,overlap);
npatch_fov = (sizY(1)/patch_size(1))^2;
K = round(options.maxcells/npatch_fov);   % number of components to be found in each patch
tau = options.gSig;                       % std of gaussian kernel (size of neuron) 
p = options.p;                            % order of AR model dynamics 

%% Run on patches

% [A,b,C,f,S,P,RESULTS,YrA] = run_CNMF_patches(imG,K,patches,tau,p,options);
[A,b,C,f,~,P,~,YrA] = run_CNMF_patches(imG,K,patches,tau,p,options);

%% classify components 

rval_space = classify_comp_corr(imG,A,C,b,f,options);
ind_corr = rval_space > options.space_thresh;           % components that pass the space correlation test

try  % matlab 2017b or later is needed for the CNN classifier
    % [ind_cnn,value] = cnn_classifier(A,[options.d1,options.d2],'cnn_model',options.cnn_thr);
    [ind_cnn,~] = cnn_classifier(A,[options.d1,options.d2],'cnn_model',options.cnn_thr);
catch
    ind_cnn = true(size(A,2),1);
end

fitness = compute_event_exceptionality(C+YrA,options.N_samples_exc,options.robust_std); % event exceptionality
ind_exc = (fitness < options.min_fitness);

keep = (ind_corr | ind_cnn) & ind_exc;

%% run GUI for modifying component selection (optional, close twice to save values)
corr_image = correlation_image_max(imG);  % background image for plotting
run_GUI = false;
if run_GUI
    Coor = plot_contours(A,corr_image,options,1); close;
    GUIout = ROI_GUI(A,options,corr_image,Coor,keep,ROIvars);   
    options = GUIout{2};
    keep = GUIout{3};    
end

%% re-estimate temporal components
A_throw = A(:,~keep);
C_throw = C(~keep,:);
A_keep = A(:,keep);
C_keep = C(keep,:);
options.p = 2;      % perform deconvolution
P.p = 2;
[A2,b2,C2] = update_spatial_components(imG,C_keep,f,[A_keep,b],P,options);
% [C2,f2,P2,S2,YrA2] = update_temporal_components_fast(imG,A2,b2,C2,f,P,options);
[C2,f2,~,~,YrA2] = update_temporal_components_fast(imG,A2,b2,C2,f,P,options);

%% extract DF_F
%[C_df,~] = extract_DF_F(Yr,A,C,P,options);
[C_df,F0] = detrend_df_f(A2,b2,C2,f2,YrA2,options);

%% plot results
if display
    figure;
    plot_contours(A2,corr_image,options,1);
    plot_components_GUI(imG,A2,C2,b2,f2,corr_image,options);
end

%% convert contour of spatial footprints to logical masks (added by Ann Go)
masks = zeros(sizY(1),sizY(1),numel(A2));
for i = 1:numel(A2)
    masks(:,:,i) = poly2mask(A2{i}(1,:),A2{i}(2,:),d1,d2);
end

masks = logical(masks);

GUIdata.A2 = A2; 
GUIdata.C2 = C2;
GUIdata.b2 = b2;
GUIdata.f2 = f2;

end