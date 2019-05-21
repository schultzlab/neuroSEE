% Written by Ann Go
% Function for reading individual filenames for batch processing 
% INPUT
%   filename : txt file with three columns. The one-line header specifies
%               column names        
%               e.g.    file    1or2d   fov
%                       20181019_12_23_34   1   330
%                       20181019_23_34_56   2   490
% OUTPUT
%   files       : filenames (1st column)
%   filesDim    : no. of dimensions of environment (2nd column)
%   filesFOV    : size of one length of square FOV in um (3rd column)

function files = extractFilenamesFromTxtfile(filename)

    fileID = fopen(filename);
    C = textscan(fileID,'%s %s %s');
    fclose(fileID);
    k = numel(C{1})-1;
    
    files = blanks(length(C{1}{2}));
    filesDim = zeros(k,1);
    filesFOV = zeros(k,1);
    for i = 1:k
        files(i,:) = C{1}{i+1};
        filesDim(i,:) = str2double(C{2}{i+1});
        filesFOV(i,:) = str2double(C{3}{i+1});
    end
end