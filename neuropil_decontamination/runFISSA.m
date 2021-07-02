% Written by Katie Davey, edited by Ann Go and Sihao Lu
% 
% Generates neuroSEE ROI files & call FISSA python code

function runFISSA( masks, tiffile, fissadir )

roidir = [fissadir 'roi/'];
% if roi directory doesn't exist, create it
if ~exist( roidir, 'dir' )
   mkdir( roidir ); fileattrib( roidir,'+w','g','s' );
end

outdir = [fissadir 'FISSAout/'];
% if output directory doesn't exist, create it
if ~exist( outdir, 'dir' )
   mkdir( outdir ); fileattrib( outdir,'+w','g','s' );
end

%% Generate ROI files
% roi's in masks, with format szx x szy x cell (e.g. 512 x 512 x 101)

N = size(masks,3); % number of ROIs
roizip   = fullfile( fissadir, 'rois.zip' );
roifiles = cell(N,1); % list of roi files
for i=1:N
   mask  = masks(:,:,i); 
   % grow mask by 1 pixel to include entire cell within the mask
   elmt  = strel('disk', 1); 
   mask  = imdilate( mask, elmt );
   fname = sprintf('%d.roi', i); % filename for current mask
   
   % generate ROI files by manually writing binary files 
   [~, roi_fname] = exportMask2ROI( roidir, fname, mask );   
   roifiles{i} = roi_fname;
end

% Generate zip file of all ROIs - note that doing this manually in Finder
% (assuming a mac imagecomputer) generated a zip file that failed when being
% input to the fissa code, but I've no idea why?!
zip( roizip, roifiles );

%% Run FISSA
% Now we need to run fissa using python
mydir  = pwd;
ind   = strfind(mydir,'/');
newdir = mydir(1:ind(end)-1);
folder = fullfile(newdir(1:ind(end)-1),'/neuropil_decontamination');
pyfun = [folder '/runFISSA.py' ' ' tiffile ' ' roizip ' ' outdir];
if exist('/Users/mgo/anaconda3/envs/neuroSEE/','dir')
    python_executable = '/Users/mgo/anaconda3/envs/neuroSEE/bin/python';
elseif exist('/home/mgo/anaconda3/envs/neuroSEE/','dir')
    python_executable = '/home/mgo/anaconda3/envs/neuroSEE/bin/python';
else
    python_executable = '/rds/general/user/mgo/home/anaconda3/envs/neuroSEE/bin/python';
end
pystr = [python_executable ' ' pyfun];
% [status, pyout] = 
system( pystr )



