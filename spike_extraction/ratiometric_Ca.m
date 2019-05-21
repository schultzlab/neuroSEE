% Written by Ann Go, adapted from Katie's ratiometric_Ca.m
%
% R = ratiometric_Ca( green, red, sm, scanrate )
% Inputs: 
%   green - green fluorsecence in format: neuron x time 
%     red - green fluorsecence in format: neuron x time 
%  smooth - if non-zero signals are smoothed by given amount over time
%           must be a positive integer (uses matlab's smooth function)
% scanrate - frequency of scanning (needed if smoothing is non-zero)
%
% See also:       smooth 

function R = ratiometric_Ca( green, red, sm )
   if nargin==0, help ratiometric_Ca; return; end
   nN = size(green,1);  % number of cells 
   
   if nargin==3 
      sm = 11;
   end
   
%    if sm>0
%       % for each cell smooth red & green timeseries over time
%       for c=1:nN
%           green(c,:) = medfilt1(green(c,:), sm); %smooth(green(c,:), sm);
%           red(c,:)   = medfilt1(red(c,:), sm); %smooth(red(c,:), sm);
%       end
%    end
   
   R  = green./red;

   % We want R0 to equal the R_ valid that is most often repeated, as we're
   % assuming that this represents the baseline R_ value, when the cell is
   % not activated. Use mode, but since the data is continuous we need to
   % histogram it first. 
   [bincounts, bins] = hist( R, 1e2 );
   [~, i] = max( bincounts );
   R0     = bins( i );
   R0_    = repmat( R0', [nN 1] );
   R      = (R - R0_) ./ R0_;
   
  if sm>0
      % for each cell smooth ratiometric timeseries 
      for c=1:nN
          R(c,:) = medfilt1(R(c,:), sm); 
      end
   end
   
end