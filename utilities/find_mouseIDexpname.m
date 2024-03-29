function [mouseid,expname,fov] = find_mouseIDexpname( list )

%% Find mouseid and experiment name
ind = strfind(list,'_');
mouseid = list(ind(1)+1:ind(2)-1);
expname = list(ind(2)+1:end-4);

ind = strfind(list,'fov');
if ~isempty(ind)
    fov = list(ind:ind+3);
else
    fov = '';
end
