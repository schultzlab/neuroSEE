function [mouseid,expname] = find_mouseIDexpname( list )

%% Find mouseid and experiment name
tf = zeros(size(list));
for ii = 1:numel(list)  
    tf(ii) = strcmpi(list(ii),'_');
end
ind = find(tf);
mouseid = list(ind(1)+1:ind(2)-1);

switch numel(ind)
    case 2
        expname = list(ind(2)+1:end-4);
    case 3
        expname = list(ind(2)+1:end-4);
    case 4
        expname = list(ind(2)+1:ind(3)-1);
    otherwise
        expname = list(ind(2)+1:ind(3)-1);
end
