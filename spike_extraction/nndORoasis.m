%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   NeuroSEE: An automated Neuronal Source Extraction
%             and Exploration toolbox
%   
%   Author: Seigfred Prado   
%   Supervisor: Simon Schultz
%   Acknowledgment: Stephanie Reynolds, Pier Luigi Dragotti
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Non-Negative Derivative Spike Estimation
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   INPUT:
%       R: ratiometric time series whose size is [cells x samples]
%       thr: threshold value (usually 0.022)
%       mode: 1 if using basic NND, 2 if using OASIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function spikes = nndORoasis(R, mode, g1, lam1)
   if nargin<4, lam1 = 2.4; end
   if nargin<3, g1 = 0.91; end
   
   % Normalise the Ca traces first
   %nR = bsxfun(@rdivide, bsxfun(@minus, R, min(R,[],2)),(max(R,[],2)-min(R,[],2)));
   nR = R;
   nC = size(nR,1); % number of cells
   mAlength = 5;    % smoothing window?
   thr = 0.022; %0.009
   spikes = zeros(size(nR));

   if mode == 1 % Use basic NND
      for ci=1:nC
         nRs  = conv( nR(ci,:), ones(mAlength,1)/mAlength, 'same' );
         temp = diff(nRs);
         temp = [temp, temp(end)];
         temp(temp<thr) = 0;
         spikes(ci,:) = temp;
      end
   else % Use OASIS
      deconvR = zeros(size(R));
      % g1   = 0.91; %0.95
      % lam1 = 2.4; %2.4
      for ci=1:nC
         [~, deconvR(ci,:)] = oasisAR1(nR(ci,:), g1, lam1);
      end
       spikes = deconvR;
   end
end