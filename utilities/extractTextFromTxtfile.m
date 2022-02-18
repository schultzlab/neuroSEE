% Written by Ann Go
% Function for reading text (filenames or plain text) for batch processing 
% INPUT
%           e.g.    20181019_12_23_34
%                   14.56
%                       
% OUTPUT
%   files       : text 

function lines = extractTextFromTxtfile(textfile)

    fileID = fopen(textfile);
    lines = textscan(fileID,'%s','delimiter','\n');
    fclose(fileID);
    lines = lines{1}; 
end