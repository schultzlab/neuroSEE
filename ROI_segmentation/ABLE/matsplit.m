% varargout = matsplit( x )
% Allocate each column of a matrix to a separate output. If x is a vector,
% allocate each element of the vector to a separate output.
% E.g.  [a, b, c, d] = matsplit( [1 2 3 4] );
%
% See also:    deal 
function varargout = matsplit( x )
   if isvector( x ), x = x(:)'; end
   varargout = cell(1, nargout );
   for vi=1:min( length(varargout), size(x,2) )
      varargout{vi} = x(:,vi);
   end
end
