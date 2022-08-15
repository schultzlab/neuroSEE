% Adapted by Ann Go from CaImAn's demo_script.m
%
% This function extracts ROIs and their decontaminated signals from a green
% channel image stack using CaImAn

function [tsG, df_f, masks, corr_image, F0, A2] = CaImAn( imG, options, display )

if nargin<3, display = false; end
if nargin<2 
    params = neuroSEE_setparams();
    options = params.ROIsegment.CaImAn;
end

fprintf('\tExtracting ROIs with CaImAn\n');

% imG = imG - min(imG(:)); 
if ~isa(imG,'single');    imG = single(imG);  end     % convert to single

[d1,d2,T] = size(imG);                                % dimensions of dataset
d = d1*d2;                                            % total number of pixels

%% Set parameters
options.d1 = d1;
options.d2 = d2;
K = options.maxcells;                                 % number of components to be found
tau = options.gSig;                                   % std of gaussian kernel (half size of neuron) 
p = options.p;

%% Data pre-processing

[P,imG] = preprocess_data(imG,p);

%% fast initialization of spatial components using greedyROI and HALS

[Ain,Cin,bin,fin,center] = initialize_components(imG,K,tau,options,P);  % initialize

% display centers of found components
corr_image =  correlation_image(imG); %reshape(P.sn,d1,d2);  %max(Y,[],3); %std(Y,[],3); % image statistic (only for display purposes)
if display
    % fh1 = 
    figure; imagesc(corr_image);
        axis equal; axis tight; hold all;
        scatter(center(:,2),center(:,1),'mo');
        title('Center of ROIs found from initialization algorithm');
        drawnow;
end


%% manually refine components (optional)
refine_components = false;  % flag for manual refinement
if refine_components
    [Ain,Cin,center] = manually_refine_components(imG,Ain,Cin,center,corr_image,tau,options);
end
    
%% update spatial components
Yr = reshape(imG,d,T);
[A,b,Cin] = update_spatial_components(Yr,Cin,fin,[Ain,bin],P,options);

%% update temporal components
P.p = 0;    % set AR temporarily to zero for speed
[C,f,P,S,YrA] = update_temporal_components(Yr,A,b,Cin,fin,P,options);

%% classify components

rval_space = classify_comp_corr(imG,A,C,b,f,options);
ind_corr = rval_space > options.space_thresh;           % components that pass the correlation test
                                        % this test will keep processes
                                        
%% further classification with cnn_classifier
try  % matlab 2017b or later is needed
    % [ind_cnn,value] = cnn_classifier(A,[d1,d2],'cnn_model',options.cnn_thr);
    [ind_cnn,~] = cnn_classifier(A,[d1,d2],'cnn_model',options.cnn_thr);
catch
    ind_cnn = true(size(A,2),1);                        % components that pass the CNN classifier
end     
                            
%% event exceptionality

fitness = compute_event_exceptionality(C+YrA,options.N_samples_exc,options.robust_std);
ind_exc = (fitness < options.min_fitness);

%% select components

keep = (ind_corr | ind_cnn) & ind_exc;

%% display kept and discarded components
A_keep = A(:,keep);
C_keep = C(keep,:);
if display
    % fh2 = 
    figure;
        subplot(121); montage(extract_patch(A(:,keep),[d1,d2],[30,30]),'DisplayRange',[0,0.15]);
            title('Kept Components');
        subplot(122); montage(extract_patch(A(:,~keep),[d1,d2],[30,30]),'DisplayRange',[0,0.15])
            title('Discarded Components');
end
        
%% merge found components
% [Am,Cm,K_m,merged_ROIs,Pm,Sm] = merge_components(Yr,A_keep,b,C_keep,f,P,S,options);
[Am,Cm,K_m,merged_ROIs,Pm,~] = merge_components(Yr,A_keep,b,C_keep,f,P,S,options);

%%
display_merging = display; % flag for displaying merging example
if and(display_merging, ~isempty(merged_ROIs))
    i = 1; %randi(length(merged_ROIs));
    ln = length(merged_ROIs{i});
    % fh3 = 
    figure;
        set(gcf,'Position',[300,300,(ln+2)*300,300]);
        for j = 1:ln
            subplot(1,ln+2,j); imagesc(reshape(A_keep(:,merged_ROIs{i}(j)),d1,d2)); 
                title(sprintf('Component %i',j),'fontsize',16,'fontweight','bold'); axis equal; axis tight;
        end
        subplot(1,ln+2,ln+1); imagesc(reshape(Am(:,K_m-length(merged_ROIs)+i),d1,d2));
                title('Merged Component','fontsize',16,'fontweight','bold');axis equal; axis tight; 
        subplot(1,ln+2,ln+2);
            plot(1:T,(diag(max(C_keep(merged_ROIs{i},:),[],2))\C_keep(merged_ROIs{i},:))'); 
            hold all; plot(1:T,Cm(K_m-length(merged_ROIs)+i,:)/max(Cm(K_m-length(merged_ROIs)+i,:)),'--k')
            title('Temporal Components','fontsize',16,'fontweight','bold')
        drawnow;
end

%% refine estimates excluding rejected components

Pm.p = 0;    % restore AR value
[A2,b2,C2] = update_spatial_components(Yr,Cm,f,[Am,b],Pm,options);
[C2,f2,P2,S2,YrA2] = update_temporal_components(Yr,A2,b2,C2,f,Pm,options);


%% do some plotting

% [A_or,C_or,S_or,P_or] = order_ROIs(A2,C2,S2,P2); % order components
% K_m = size(C_or,1);
% [df_f,F0] = extract_DF_F(Yr,A_or,C_or,P_or,options); % extract DF/F values (optional)
% tsG = C_or; 
[df_f,F0] = detrend_df_f(A2,b2,C2,f2,YrA2,options);
tsG = C2; 

% fh4 = 
fig = figure;
% [Coor,json_file] = plot_contours(A_or,corr_image,options,0); % contour plot of spatial footprints
[Coor,~] = plot_contours(A2,corr_image,options,0); % contour plot of spatial footprints
% savejson('jmesh',json_file,'filename');        % optional save json file with component coordinates (requires matlab json library)

if ~display
    close(fig);
end

%% display components

% if display
%     plot_components_GUI(Yr,A_or,C_or,b2,f2,Cn,options);
% end

%% make movie
% if (0)  
%     make_patch_video(A_or,C_or,b2,f2,Yr,Coor,options)
% end

%% convert contour of spatial footprints to logical masks (added by Ann Go)
masks = zeros(d1,d2,numel(Coor));
for i = 1:numel(Coor)
    masks(:,:,i) = poly2mask(Coor{i}(1,:),Coor{i}(2,:),d1,d2);
end

masks = logical(masks);

% GUIdata.Yr = Yr; 
% GUIdata.A_or = A_or;
% GUIdata.C_or = C_or;
% GUIdata.b2 = b2;
% GUIdata.f2 = f2;

end