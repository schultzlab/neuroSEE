% From CaImAn
% Edit by Ann Go: line 4, 25

function A = read_raw_file(filename,framenum,window,FOV,bitsize,skipN)

%% Starting from frame number: 'framenum', read 'window'
fid = fopen(filename);
imsize = FOV(1)*FOV(2)*bitsize;                                                   % Bit size of single frame

current_seek = ftell(fid);
fseek(fid, 0, 1);
file_length = ftell(fid);
fseek(fid, current_seek, -1);
frame_length = file_length/imsize;

if bitsize == 2
    bstring = 'uint16=>single';
else
    bstring = 'real*4=>double';
end
%%
window = min(window,frame_length);
A = zeros(FOV(1),FOV(2),window/skipN,'single');

for w = 0:skipN:window-1                                                          % For each frame inside the window
    fseek(fid,(framenum-1+w)*imsize,'bof');                                 % Position the file-indicator at the beginning of the frame
    temp = fread(fid,[FOV(1),FOV(2)],bstring,0,'b');                              % Read the frame
    A(:,:,w/skipN+1) = temp';                                                                    
end
fclose(fid);