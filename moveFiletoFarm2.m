% Written by Ann Go
% Run in linux box after processing to move files from local hard drive to thefarm2

files = dir('*.tif');
filestomove = files(6:end);

for i = 1:length(filestomove)
    file = filestomove(i).name;
    dest_folder = ['/rds/general/user/mgo/projects/thefarm2/live/CrazyEights/AD_2PCa/Data/' file(1:8) ...
        '/Processed/' file];
    copyfile( file, dest_folder );
    % delete file;
end