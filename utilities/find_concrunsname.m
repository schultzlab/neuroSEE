function concrunsname = find_concrunsname( list )

%% Find mouseid and experiment name
ind = strfind(list,'_');
expname = list(ind(2)+1:end-4);

ind = strfind(expname,'-');
if ~isempty(ind)
    concrunsname = expname(1:ind-1);
else
    concrunsname = expname;
end
