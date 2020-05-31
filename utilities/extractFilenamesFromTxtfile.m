% Written by Ann Go
% Function for reading individual filenames for batch processing 
% INPUT
%   filename 
%           e.g.    20181019_12_23_34
%                       
% OUTPUT
%   files       : filenames 

function files = extractFilenamesFromTxtfile(filename)

    fileID = fopen(filename);
    C = textscan(fileID,'%s');
    fclose(fileID);
    k = numel(C{1});
    
    files = blanks(length(C{1}{2}));
    for i = 1:k
        files(i,:) = C{1}{i};
    end
end