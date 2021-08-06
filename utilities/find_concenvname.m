function concenvname = find_concenvname( list )

%% Find mouseid and experiment name
ind = strfind(list,'_');
expname = list(ind(2)+1:end-4);

ind = strfind(expname,'-');
if ~isempty(ind)
    concenvname = expname(1:ind-1);
else
    concenvname = expname;
end
