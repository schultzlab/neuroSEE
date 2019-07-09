%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   NeuroSEE: An automated Neuronal Source Extraction
%             and Exploration toolbox
%   
%   Author: Seigfred Prado   
%   Supervisor: Simon Schultz
%   Acknowledgment: Stephanie Reynolds, Pier Luigi Dragotti
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function [cellMasks, cellTimeSeriesG, nhbdTimeSeriesG, cellTimeSeriesR, nhbdTimeSeriesR] = ...
%                                segment_cells(phi_0, videoG, imageR, radius, options)

function [cellMasks, cellTimeSeriesG] = segment_cells(phi_0, videoG, radius, options) % edit by Ann

%%%% Fixed algorithm parameters                       
c0                 = 3;
epsilon            = 2;
timestep           = 10;
delta              = 1;
mu                 = 0.2/timestep;
noChangeNum        = 20;
bandWidth          = 3 * radius;

%%% Optional arguments  
if isfield(options, 'lambda')
    lambda = options.lambda / timestep;
else
    lambda = 100 / timestep;
end
if ~isfield(options, 'metric')
    options.metric = 'corr';
end
if isfield(options, 'minimumSize')
    minimumSize = options.minimumSize;
else
    minimumSize  = 3; 
end 
if isfield(options, 'maximumSize')
    maximumSize = options.maximumSize;
else
    maximumSize  = round(pi * radius^2 * 4); 
end 
if isfield(options, 'mergeCorr')
    mergeCorr = options.mergeCorr;
elseif isfield(options, 'snr')
    mergeCorr = 0.8/(1+(10^(-snr/10)));
else
    mergeCorr = 0.8;
end
if isfield(options, 'mergeTime')
    mergeTime = options.mergeTime;
elseif strcmp(options.metric, 'corr')
    mergeTime = 'during';
else 
    mergeTime = 'atEnd';
end 
if isfield(options, 'maxIt')
    maxIt  = options.maxIt;
else
    maxIt  = 150;
end 

% Initialise variables
it             = 1;
t_len          = size(videoG,3);
cell_num       = size(phi_0,3);
active         = true(cell_num,1);
activity       = inf*ones(cell_num,noChangeNum);
cellIndex      = 1:cell_num;
phi            = phi_0;
video_dim      = [size(videoG,1), size(videoG,2)];
timeSeriesG     = zeros(cell_num, t_len, 'single');
%%%%added codes by svp
% timeSeriesR  = zeros(cell_num, t_len, 'single');
%%%%
region         = zeros([video_dim, cell_num]);
se_narrowband  = strel('square', round(bandWidth));
se_overlap     = strel('square', ceil(radius/2));
video_reshapedG = reshape(videoG,[video_dim(1)*video_dim(2)*t_len,1,1]);
%vid_ratio = (videoG-videoR)./videoR;
%%%added code by svp
% image_reshapedR = reshape(imageR,[video_dim(1)*video_dim(2),1,1]);
%%%

mods           = (0:1:(t_len - 1))*video_dim(1)*video_dim(2);
merge_counter  = 1;
          
%%%%%%%%%   Initialise time series, phi and narrowband of each cell   %%%%%
ii = 1;
while ii <= size(phi,3)
     
   currentPhi          = phi(:,:,ii);
   currentMask         = currentPhi<0;
   % [timeSeriesG(ii,:), timeSeriesR(ii,:)] = extractTimeSeries(currentMask, mods, video_reshapedG, image_reshapedR, t_len);
   timeSeriesG(ii,:) = extractTimeSeries(currentMask, mods, video_reshapedG, t_len);
   phi(:,:,ii)         = maskToSignedDistanceFunction(currentPhi, c0);
    
   if ~any(currentMask)
        phi(:,:,ii)            = [];
        timeSeriesG(ii,:)       = [];
        % timeSeriesR(ii,:) = []; %%%added by svp
        region(:,:,ii)         = [];
        active(ii)             = [];
        activity(ii, :)        = [];
        cellIndex(ii)          = [];       
   else
       ii                   = ii + 1;
   end
    
end

nhbdTimeSeriesG = zeros(cell_num, t_len);
%nhbdTimeSeriesR = zeros(cell_num, t_len);


%%% If we're plotting progress, initialise the plots
% if isfield(options, 'plot_progress')
%     figure('units','normalized','outerposition',[0 0 1 1])
%     corrIm = options.corrIm;
%     meanIm = options.meanIm;
%     if isfield(options, 'cell_to_monitor')
%         cell_to_monitor  = options.cell_to_monitor;
%     else
%         cell_to_monitor  = round(size(phi_0,3)/2);
%     end
%     ncols = 5;
%     nrows = 4;
% end
 
%%% UPDATE: While iterations less than maximum and cells active %%%%%%%%%%%
prevstr = [];
while and(it <= maxIt, any(active))
     
    phi   = NeumannBoundCond(phi, video_dim); % Make phi satisfy neumann boundary conditions
    ii    = 1;
    
    % Are any cells touching and v. correlated? If yes, merge them
    if and(strcmp(mergeTime, 'during'), length(cellIndex) > 1)
        correlation       = corrcoef(timeSeriesG');
        correlation       = tril(correlation, -1); % get rid of duplicates (it is symmetric)
        [p, q]            = find(correlation > mergeCorr); 
        while and(~isempty(p), ~isempty(q))
            % if they are overlapping (use square shape dilation
            % ebcause its faster)
            overlap = and( imdilate(phi(:,:,p(1))<0,se_overlap), ...
                           imdilate(phi(:,:,q(1))<0,se_overlap));
            if any(overlap(:))
                ordering      = sort([p(1),q(1)]); % we merge to the first one
                merge_remove  = ordering(2); 
                merge_keep    = ordering(1);
                merge_counter = merge_counter + 1;

                [active, activity, phi, region,...
                 timeSeriesG, cellIndex] =...
                 merge(merge_keep, merge_remove, phi, c0, timeSeriesG, ...
                       video_reshapedG, t_len,...
                       mods, se_narrowband,...
                       active, activity,...
                       region, cellIndex); 

                p(1)               = [];
                q(1)               = [];
                p(p==merge_remove) = merge_keep;
                q(q==merge_remove) = merge_keep;
                p(p>merge_remove)  = p(p>merge_remove) - 1;
                q(q>merge_remove)  = q(q>merge_remove) - 1;
                unique_pq          = unique([p, q], 'rows');
                if ~isempty(unique_pq)
                    p                  = unique_pq(:,1);
                    q                  = unique_pq(:,2);
                end

            else
                p(1) = [];
                q(1) = [];
            end
        end    
    end
    
    while ii <= size(phi,3)
         
       if active(ii)  
           
           phi_removed    = 0;
           phiUpdate      = phi(:,:,ii);
           nhbd_G           = imdilate(phiUpdate<0, se_narrowband);
           region(:,:,ii) = nhbd_G;
           otherCells     = zeros(video_dim);
           ind            = find(nhbd_G);
           [ind_x, ind_y] = ind2sub(video_dim,ind) ;
           for jj = 1:length(ind_x)
               otherCells(ind_x(jj), ind_y(jj)) = any(phi(ind_x(jj),ind_y(jj),1:end ~= ii)<0,3);
           end
           nhbd_G             = and(nhbd_G, ~or(otherCells, phiUpdate<0));
           % [nG, nR]           = extractTimeSeries(nhbd_G, mods, video_reshapedG, image_reshapedR, t_len);
           nG               = extractTimeSeries(nhbd_G, mods, video_reshapedG, t_len);
           diracPhi         = Dirac(phiUpdate, epsilon); 
           evaluate_idx     = diracPhi > 0;
           onlyThisMask     = and(~otherCells, phiUpdate<0);
           % [timeSeriesG(ii,:), timeSeriesR(ii,:)] = extractTimeSeries(onlyThisMask, mods,...
           %                                           video_reshapedG, image_reshapedR, t_len);
           timeSeriesG(ii,:) = extractTimeSeries(onlyThisMask, mods, video_reshapedG, t_len);

            
           if nnz(evaluate_idx)
                [similarityVelocity] = ...
                cellSimilarityVelocity(video_dim, phi, nG, ii, evaluate_idx,...
                                       mods, video_reshapedG,  t_len, ...
                                       timeSeriesG, otherCells);
           else
                similarityVelocity = zeros(size(diracPhi));
           end
           
           %%% Calculate regularizer
           regularisationVelocity = zeros(video_dim);
           [x_loc, y_loc]         = find(region(:,:,ii));
           x_min                  = max(min(x_loc) - 2, 1);
           x_max                  = min(max(x_loc) + 2, video_dim(2));
           y_min                  = max(min(y_loc) - 2, 1);
           y_max                  = min(max(y_loc) + 2, video_dim(2));
           if and(x_max > x_min, y_max > y_min) && and(x_max<size(videoG,1), y_max < size(videoG,2))
              regularisationVelocity(x_min:x_max, y_min:y_max) =...
                distanceRegularisation(phiUpdate(x_min:x_max, y_min:y_max));
           end
          
           %%% Update phi
           phiUpdate              = phiUpdate + ...
                                    timestep * region(:,:,ii).*...
                                    (mu * regularisationVelocity - ...
                                    lambda*diracPhi.*...
                                    similarityVelocity);
                                
           remove       = ~any(phiUpdate(:)<0);
           onlyThisMask = and(~otherCells, phiUpdate<0);
           mask_size    = nnz(phi(:,:,ii)<0);                            
 
           
           if and(and(~remove, any(onlyThisMask(:))),...
                  and(mask_size < maximumSize, mask_size > minimumSize))
 
               %%% Update stored values
               changed            = nnz(or(and(phi(:,:,ii)<0, phiUpdate>0),...
                                           and(phi(:,:,ii)>0, phiUpdate<0)));
               activity(ii,:)     = [activity(ii,2:end), changed];
               phi(:,:,ii)        = phiUpdate; 
 
               %%% If converged, stop it
               if and(it > noChangeNum, activity(ii,:) < delta)
                    active(ii)       = 0;
                    % [timeSeriesG(ii,:), timeSeriesR(ii,:)] = extractTimeSeries(onlyThisMask, mods,video_reshapedG, image_reshapedR, t_len);
                    timeSeriesG(ii,:) = extractTimeSeries(onlyThisMask, mods,video_reshapedG, t_len);
                    nhbdTimeSeriesG(cellIndex(ii),:) = nG;
					% nhbdTimeSeriesR(cellIndex(ii),:) = nR;
               end
           else 
               cellIndex(ii)        = [];
               phi_removed          = 1;
               phi(:,:,ii)          = [];
               timeSeriesG(ii,:)     = [];
               % timeSeriesR(ii,:) = []; %%%added by svp
               region(:,:,ii)       = [];
               active(ii)           = [];
               activity(ii, :)      = [];
           end
%            tempPhi = phi;
%            if isfield(options, 'plot_progress')
%                
%                if cellIndex(ii) == cell_to_monitor
%                       updatePlot(ii, videoG, phi, radius, nrows, ncols, corrIm,...
%                          meanIm, onlyThisMask, nhbd_G, region, diracPhi,...
%                          similarityVelocity, regularisationVelocity,...
%                          mu, mods, video_reshapedG, t_len, timeSeriesG,...
%                          it, activity, delta, noChangeNum)   
%                       drawnow;
%                end
%            end 
            
           %%% If it's the final iteration, store the neighbouring tsG
%            if it == maxIt
%                nhbdTimeSeriesG(cellIndex(ii),:) = nG;
% 			   nhbdTimeSeriesR(cellIndex(ii),:) = nR;
%            end
           
           
           if ~phi_removed
               ii = ii + 1;
           end
            
           
       else
           %%%%% Update time series if any cell has moved into it
           otherCells     = zeros(video_dim);
           [ind_x, ind_y] = find(phi(:,:,ii)<0);
           for jj = 1:length(ind_x)
               otherCells(ind_x(jj), ind_y(jj)) = any(phi(ind_x(jj),ind_y(jj),1:end ~= ii)<0,3);
           end
           onlyThisMask     = and(~otherCells, phi(:,:,ii)<0);
           % [timeSeriesG(ii,:), timeSeriesR(ii,:)] = extractTimeSeries(onlyThisMask, mods,video_reshapedG, image_reshapedR, t_len);
           timeSeriesG(ii,:) = extractTimeSeries(onlyThisMask, mods,video_reshapedG, t_len);
           ii               = ii + 1;
       end

    end
    % disp(['Iteration ', num2str(it)]);
    currstr = sprintf( '\tIteration %g', it );
    refreshdisp( currstr, prevstr );
    prevstr = currstr; 

    it = it + 1; 
end
 
% Are any cells touching and v. correlated? If yes, merge them
if and(strcmp(mergeTime, 'atEnd'), length(cellIndex) > 1)
    correlation       = corrcoef(timeSeriesG');
    correlation       = tril(correlation, -1); % get rid of duplicates (it is symmetric)
    [p, q]            = find(correlation > mergeCorr); 
    while and(~isempty(p), ~isempty(q))
        % if they are overlapping (use square shape dilation
        % ebcause its faster)
        overlap = and( imdilate(phi(:,:,p(1))<0,se_overlap), ...
                       imdilate(phi(:,:,q(1))<0,se_overlap));
        if any(overlap(:))
            ordering      = sort([p(1),q(1)]); % we merge to the first one
            merge_remove  = ordering(2); 
            merge_keep    = ordering(1);
            merge_counter = merge_counter + 1;

%             [active, activity, phi, region,...
%              timeSeriesG, timeSeriesR, cellIndex] =...
%              merge(merge_keep, merge_remove, phi, c0, timeSeriesG, timeSeriesR,...
%                    video_reshapedG, image_reshapedR, t_len,...
%                    mods, se_narrowband,...
%                    active, activity,...
%                    region, cellIndex); 
             [active, activity, phi, region,...
             timeSeriesG, cellIndex] =...
             merge(merge_keep, merge_remove, phi, c0, timeSeriesG, ...
                   video_reshapedG, t_len,...
                   mods, se_narrowband,...
                   active, activity,...
                   region, cellIndex); 
   

            p(1)               = [];
            q(1)               = [];
            p(p==merge_remove) = merge_keep;
            q(q==merge_remove) = merge_keep;
            p(p>merge_remove)  = p(p>merge_remove) - 1;
            q(q>merge_remove)  = q(q>merge_remove) - 1;
            unique_pq          = unique([p, q], 'rows');
            if ~isempty(unique_pq)
                p                  = unique_pq(:,1);
                q                  = unique_pq(:,2);
            end

        else
            p(1) = [];
            q(1) = [];
        end
    end    
end   

% nhbdTimeSeriesG = nhbdTimeSeriesG(cellIndex,:);
% nhbdTimeSeriesR = nhbdTimeSeriesR(cellIndex,:);
cellTimeSeriesG = timeSeriesG;
% cellTimeSeriesR = timeSeriesR;
cellMasks      = phi < 0;
% disp(['\tFinished after ', num2str(it - 1), ' iterations']);
currstr = sprintf( '\tFinished after %g iterations\n', it-1 );
refreshdisp( currstr, prevstr );


function [] = updatePlot(ii, videoG, phi, radius, nrows, ncols, corrIm, meanIm,...
                         onlyThisMask, nhbd_G, region, diracPhi,...
                         similarityVelocity, regularisationVelocity,...
                         mu, mods, video_reshapedG, t_len, timeSeriesG,...
                         it, activity, delta, noChangeNum)
     
    t = 1:size(videoG,3);
    [x,y] = find(phi(:,:,ii)<0);
    d = 2*radius;
    l_x = max(1, min(x) - d);
    u_x = min(size(videoG,1), max(x) + d);
    l_y = max(1, min(y) - d);
    u_y = min(size(videoG,2), max(y) + d);
     
    %%%% Mean image plot
    h(1) = subplot(nrows, ncols, 1);
    imagesc(corrIm(l_x:u_x,l_y:u_y));
    title('Correlation image')
 
    %%%% Mean image and initial cell interior
    h(2) = subplot(nrows, ncols, 2);
    imagesc(meanIm(l_x:u_x,l_y:u_y));
    title('Mean image')
 
    %%%% Mean image and current cell interior
    h(3) = subplot(nrows, ncols, 3);   
    imshowpair(corrIm(l_x:u_x,l_y:u_y), phi(l_x:u_x,l_y:u_y,ii)<0,...
               'ColorChannels', [1, 2, 2], 'Scaling', 'none')
    title('Interior');
 
    %%%% Mean image and current cell interior
    h(4) = subplot(nrows, ncols, 4);   
    imshowpair(corrIm(l_x:u_x,l_y:u_y), onlyThisMask(l_x:u_x,l_y:u_y),...
               'ColorChannels', [1, 2, 2], 'Scaling', 'none')
    title('Interior excluding other cells');
 
     
    %%%% Mean image and narrowband around cell
    h(5) = subplot(nrows,ncols,5);   
    imshowpair(corrIm(l_x:u_x,l_y:u_y), nhbd_G(l_x:u_x,l_y:u_y),...
               'ColorChannels', [1, 2, 2], 'Scaling', 'none')
    title('Narrowband')
 
    %%%% Phi in 3D
    subplot(nrows,ncols,6);
    mesh(-phi(l_x:u_x,l_y:u_y,ii));   % for a better view, the LSF is displayed upside down
    hold on;  contour(phi(l_x:u_x,l_y:u_y,ii), [0,0], 'r','LineWidth',2);
    title('Level set function');
 
    %%%% Other cell interiors
    if size(phi,3)>1
        h(6) = subplot(nrows, ncols, 7);   
        imshowpair(corrIm(l_x:u_x,l_y:u_y), any(phi(l_x:u_x,l_y:u_y,1:end~=ii)<0,3),...
               'ColorChannels', [1, 2, 2], 'Scaling', 'none')
        title('Other cells');
    end
     
    %%%% Cell similarity velocity    
    h(7) = subplot(nrows,ncols,8);
    V    = lambda * region(l_x:u_x,l_y:u_y,ii).*...
           diracPhi(l_x:u_x,l_y:u_y).*similarityVelocity(l_x:u_x,l_y:u_y);
    imagesc(V);
    c = colorbar;
    cmap           = cmocean('balance');
    u_lim = max(max(abs(V(:))),eps);
    set(gca,'CLim',[-u_lim u_lim])
    colormap(cmap);
    title('Data-driven velocity');
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    set(c, 'ytick', 0);
    
    %%%% Regularisation velocity
    h(8) = subplot(nrows,ncols,9);
    V = -mu*region(l_x:u_x,l_y:u_y,ii).*...
           regularisationVelocity(l_x:u_x,l_y:u_y);
    imagesc(V);
    u_lim = max(max(abs(V(:))),eps);
    set(gca,'CLim',[-u_lim u_lim])
    c = colorbar;
    title('Regularisation velocity');
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    set(c, 'ytick', 0);
    
    %%%% Combined velocity
    h(9) = subplot(nrows, ncols, 10);
    V = -region(l_x:u_x,l_y:u_y,ii).*...
                                    (mu * regularisationVelocity(l_x:u_x,l_y:u_y) - ...
                                     lambda * diracPhi(l_x:u_x,l_y:u_y).*...
                                     (similarityVelocity(l_x:u_x,l_y:u_y))); 
    imagesc(V);
    u_lim = max(max(abs(V(:))),eps);
    set(gca,'CLim',[-u_lim u_lim])
    c = colorbar;
    set(c, 'ytick', 0);
    title('Combined velocity');
     set(gca,'xtick',[])
    set(gca,'ytick',[])    
 
    %%%% Time series plot: time series of cell interior and narrowband
    subplot(nrows,ncols, 2*ncols + (1:ncols));
    % [n_tsG, n_tsR]    = extractTimeSeries(nhbd_G, mods, video_reshapedG, image_reshapedR, t_len);
    n_tsG    = extractTimeSeries(nhbd_G, mods, video_reshapedG, t_len);
    ts_plot = plotProperHeight( [timeSeriesG(ii,:)./norm(timeSeriesG(ii,:));...
                                n_tsG./norm(n_tsG)] );
    plot(t, ts_plot')
    hold on; 
    legend('Interior','Narrowband', 'Orientation', 'Horizontal');
    legend boxoff
    set(gca, 'xtick', []);
    set(gca, 'ytick', []);
    xl = xlabel('Time (s)');
    set(xl, 'position', get(xl,'position')+[0,0.01,0]);
    title(['Subregion time series, iteration number: ',num2str(it)]);
    box off
    hold off;
 
    %%%% Activity plot: how many pixels joined and left contour
    subplot(nrows, ncols, ncols*3 + (1:ncols));
    plot(1:noChangeNum, activity(ii,:))
    hold on;
    line([0, noChangeNum],[delta, delta], 'LineStyle','- -','Color','r');
    xlabel( ['Previous ',num2str(noChangeNum),' iterations'] );
    ylabel( 'Number of pixels' );
    title('Num. pixels that left or entered cell per iteration');
    hold off;
    box off;
     
    linkaxes([h(1), h(2), h(3), h(4), h(5), h(6), h(7), h(8) ]);
    drawnow;
    
end
 
 
function [V] = cellSimilarityVelocity(video_dim, phi, nG,...
                                    ii, narrowband,...
                                    mods,...
                                    video_reshapedG,...
                                    t_len, timeSeriesG, otherCells) 
    

    V                     = zeros(video_dim);
    a                     = squeeze(timeSeriesG(ii,:));
    narrowbandNoOverlap   = and(narrowband,~otherCells);
    narrowbandWithOverlap = and(narrowband, otherCells);
     
    % Calculate V(x) in pixels where there is no overlap
    if any(narrowbandNoOverlap(:))
        loc          = find(narrowbandNoOverlap);
        vid_loc      = bsxfun(@plus,loc, mods);
        raw_x        = single(video_reshapedG(vid_loc(:)));
        raw_x_rs     = reshape(raw_x, length(loc), t_len);
        V_no_overlap = compare(raw_x_rs, nG, a);
        V(narrowbandNoOverlap) = V_no_overlap;  
    end 
    
    % Calculate the velocity where there are overlapping cells
    if any(narrowbandWithOverlap(:))
         
        loc               = find(narrowbandWithOverlap);
        [loc_x,loc_y]     = ind2sub(video_dim,loc) ;
        vid_loc           = bsxfun(@plus,loc, mods);
        raw_x             = single(video_reshapedG(vid_loc(:)));
        raw_x_rs          = reshape(raw_x, length(loc),t_len);
 
        for j = 1:length(loc) 
            J              = find(phi(loc_x(j),loc_y(j),:) < 0);
            J(J == ii)     = [];
            b              = sum(timeSeriesG(J,:),1);  % time series of cells that x is currently in
            x              = squeeze(raw_x_rs(j,:));     % time series of current pixel
            a_plus_b       = a + b;         
            V(loc(j))      = compare(x, b, a_plus_b);
            
        end
    end
end
 
 
function dist = compare(x, a_out, a_in)
  
   dist = zeros(size(x,1),1,'single');
   
   if strcmp(options.metric, 'corr')
       dist = corr(x', a_in') - corr(x', a_out');
   elseif strcmp(options.metric, 'euclid')
       dist = sum(bsxfun(@minus, x, a_out).^2,2) - ...
              sum(bsxfun(@minus, x, a_in).^2,2);
       dist = dist./length(a_in);
   end
   
    
end
 
 
function V = distanceRegularisation(phi)
    % This Matlab code implements an edge-based active contour model as an application of the Distance Regularized Level Set Evolution (DRLSE) formulation in Li et al's paper: C. Li, C. Xu, C. Gui, M. D. Fox, "Distance Regularized Level Set Evolution and Its Application to Image Segmentation",  IEEE Trans. Image Processing, vol. 19 (12), pp.3243-3254, 2010.
    % Author: Chunming Li, all rights reserved
    % E-mail: lchunming@gmail.com   
    %         li_chunming@hotmail.com 
    % URL:  http://www.imagecomputing.org/~cmli/    
    [phi_x,phi_y] = gradient(phi);
    s   = sqrt(phi_x.^2 + phi_y.^2);
    a   = (s>=0) & (s<=1);
    b   = (s>1);
    ps  = a.*sin(2*pi*s)/(2*pi)+b.*(s-1);  % compute first order derivative of the double-well potential p2 in eqaution (16)
    dps = ((ps~=0).*ps+(ps==0))./((s~=0).*s+(s==0));  % compute d_p(s)=p'(s)/s in equation (10). As s-->0, we have d_p(s)-->1 according to equation (18)
    V   = div(dps.*phi_x - phi_x, dps.*phi_y - phi_y) + 4*del2(phi); 
    V([1,2,end-1,end], :)   = 0;
    V(:, [1,2,end-1,end])   = 0;
end
 
function f = div(nx,ny)
    % This Matlab code implements an edge-based active contour model as an application of the Distance Regularized Level Set Evolution (DRLSE) formulation in Li et al's paper: C. Li, C. Xu, C. Gui, M. D. Fox, "Distance Regularized Level Set Evolution and Its Application to Image Segmentation",  IEEE Trans. Image Processing, vol. 19 (12), pp.3243-3254, 2010.
    % Author: Chunming Li, all rights reserved
    % E-mail: lchunming@gmail.com   
    %         li_chunming@hotmail.com 
    % URL:  http://www.imagecomputing.org/~cmli/
    [nxx,~]=gradient(nx);  
    [~,nyy]=gradient(ny);
    f=nxx+nyy;
end
 
function f = Dirac(x, sigma)
    % This Matlab code implements an edge-based active contour model as an application of the Distance Regularized Level Set Evolution (DRLSE) formulation in Li et al's paper: C. Li, C. Xu, C. Gui, M. D. Fox, "Distance Regularized Level Set Evolution and Its Application to Image Segmentation",  IEEE Trans. Image Processing, vol. 19 (12), pp.3243-3254, 2010.
    % Author: Chunming Li, all rights reserved
    % E-mail: lchunming@gmail.com   
    %         li_chunming@hotmail.com 
    % URL:  http://www.imagecomputing.org/~cmli/
    f=(1/2/sigma)*(1+cos(pi*x/sigma));
    b = (x<=sigma) & (x>=-sigma);
    f = f.*b;
end
 
 
function g = NeumannBoundCond(f, video_dim)
    % This Matlab code implements an edge-based active contour model as an application of the Distance Regularized Level Set Evolution (DRLSE) formulation in Li et al's paper: C. Li, C. Xu, C. Gui, M. D. Fox, "Distance Regularized Level Set Evolution and Its Application to Image Segmentation",  IEEE Trans. Image Processing, vol. 19 (12), pp.3243-3254, 2010.
    % Author: Chunming Li, all rights reserved
    % E-mail: lchunming@gmail.com   
    %         li_chunming@hotmail.com 
    % URL:  http://www.imagecomputing.org/~cmli/
    nrow = video_dim(1); ncol = video_dim(2);
    g = f;
    g([1 nrow],[1 ncol],:) = g([3 nrow-2],[3 ncol-2],:);  
    g([1 nrow],2:end-1,:) = g([3 nrow-2],2:end-1,:);          
    g(2:end-1,[1 ncol],:) = g(2:end-1,[3 ncol-2],:);  
end
 
function tsG = extractTimeSeries(mask, mods, video_reshapedG, t_len)
         % [tsG, tsR] = extractTimeSeries(mask, mods, video_reshapedG, video_reshapedR, t_len)
    mask([1,end], :) = 0;
    mask(:, [1,end]) = 0;
    loc              = find(mask);
    vid_loc          = bsxfun(@plus,loc, mods);
    raw_tsG          = video_reshapedG(vid_loc(:));
    raw_tsG          = reshape(raw_tsG, length(loc),t_len);
    tsG              = mean(raw_tsG,1, 'double'); %% much faster to calculate as double
    %%%added codes by svp
%     raw_tsR = video_reshapedR(loc(:));
%     tsR     = mean(raw_tsR,1, 'double');
    %%%
end



% function[active, activity, phi, region, timeSeriesG, timeSeriesR, cellIndex] =...
%              merge(p, q, phi, c0, timeSeriesG, timeSeriesR, video_reshapedG, video_reshapedR, t_len,...
%                    mods, se_narrowband, active, activity,...
%                    region,  cellIndex)
function[active, activity, phi, region, timeSeriesG, cellIndex] =...
             merge(p, q, phi, c0, timeSeriesG, video_reshapedG, t_len,...
                   mods, se_narrowband, active, activity,...
                   region,  cellIndex)
                
    % Update cell ii to be the union of cells ii and jj 
    phi(:,:,p)         = -c0 * any( phi(:,:,[p,q]) < 0, 3)+...
                          c0 *~any( phi(:,:,[p,q]) < 0, 3);
    % [timeSeriesG(p,:), timeSeriesR(p,:)]    = extractTimeSeries( phi(:,:,p) < 0, mods,...
    %                                                    video_reshapedG, video_reshapedR, t_len);
    timeSeriesG(p,:)  = extractTimeSeries( phi(:,:,p) < 0, mods, video_reshapedG, t_len);
    region(:,:,p)      = imdilate(phi(:,:,p)<0, se_narrowband);
     
    % Remove cell jj's information 
    phi(:,:,q)                  = [];
    timeSeriesG(q,:)             = [];
    % timeSeriesR(q,:) = []; %%%added by svp
    active(q)                   = [];
    region(:,:,q)               = [];
    activity(cellIndex == q, :) = [];
    cellIndex                   = setdiff(cellIndex, q);
    
    %disp(['ROIs ',num2str(p), ' and ', num2str(q),' merged.'])
    str = sprintf('\tROIs %g and %g merged.\n', p, q);
    cprintf('Text',str)
end
 
 
function[sd] = maskToSignedDistanceFunction(phi, c0)
 
sd         = double((phi > 0).*(bwdist(phi < 0)-0.5) - ...
                    (phi < 0).*(bwdist(phi > 0)-0.5));  
sd         = sd*c0;
sd(sd>c0)  = c0;
sd(sd<-c0) = -c0;
 
end
 
end