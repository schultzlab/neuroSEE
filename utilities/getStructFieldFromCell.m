% vals = getStructFieldFromCell(xcell,field)
% For cell arrays of structures, grab a field of each structure & return in
% a vector called vals, if it makes sense to, else return in a cell array 
% Inputs:
%	xcell		- cell array of structures
%	field    - field name
%  type		- 'array' or 'cell' for format of reponse
% Outputs:
%	vals     - values of the field for each cell
function vals = getStructFieldFromCell(xcell,field,type)
   if nargin==0, help getStructFieldFromCell; return; end
   if nargin<3
       try
         tmp = xcell{1}.(field);
       catch
         fprintf('getStructFieldFromCell: unknown field (%s)\n\n',field);
       end
		if ~(ischar(tmp) || isstruct(tmp) || iscell(tmp))
          type = 'array';
        else
		  type = 'cell';
		end
   end
   % For smaller cell arrays, convert to structs & grab required field
   % (gotta make sure they have same fields though)
%    if strcmp(type,'array') && numel(xcell)<500 && ...
   if numel(xcell)<500 && ...
      length(unique(cellfun(@(s) length(fieldnames(s)),xcell)))==1
      % to convert cells to struct array with [c{:}] the fieldnames must
      % match, else it'll throw an error. In case of error, grab only the
      % field being requested
      try
         tmp  = [xcell{:}]; 
         if strcmpi(type,'array')
            vals    = [tmp.(field)]';
         else
            vals    = cell(size(tmp));
            vals(:) = {tmp.(field)};
         end
      catch
         tmp  = cellfun(@(x) x.(field), xcell, 'uniformoutput', 0);
         if strcmpi(type,'array')
            vals = cell2mat(tmp);
         else
            vals = tmp(:);
         end
      end
   % If it's a large cell array, don't make an entire copy of it in a
   % structure, but rather simply access only the required field
   else
      tmp  = cellfun(@(x) x.(field), xcell, 'uniformoutput', 0);
      if strcmpi(type,'array')
         vals = cell2mat(tmp);
      else
         vals = tmp(:);
      end
   end
end