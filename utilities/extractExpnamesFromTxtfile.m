% Written by Ann Go
% Function for reading individual experiment names for batch processing 
% INPUT
%   filename 
%
% OUTPUT
%   exps       : expnames 

function exps = extractExpnamesFromTxtfile(filename)

    fileID = fopen(filename);
    C = textscan(fileID,'%s');
    fclose(fileID);
    k = numel(C{1});
    
    for i = 1:k
        exps{i} = C{1}{i};
    end
end