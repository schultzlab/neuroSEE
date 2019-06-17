% Generate neuroSEE ROI files & call FISSA python code

%% Input folders & files

% mousedir  = '/Users/cedavey/Dropbox/code/hippocampus/m62/';

fissadir = [pwd '/fissa']; % test fissa directory 

% green channel tif directory - NOTE / TO DO that fissa python code seems to 
% want only relevant tif files in the directory, which will be problematic 
% with neuroSEE so the python code should be amended to allow a list of tif 
% file names as input instead of a directory
tifdir  = '/Volumes/thefarm2/live/CrazyEights/AD_2PCa/Data/20181016/'; % 20181011_15_16_17_2P_XYT_green.tif;

% name of mask generated by NeuroSEE from which ROIs will be generated
roifile = '/home/sihao/KozlovBox/Data/Fissa_test/test_files/20181016_10_11_35_2P_segment_output.mat';

% directory to store ImageJ ROI files to
roidir  = fullfile( fissadir, 'rois' ); 

% output directory to save FISSA outputs to 
outdir  = fullfile( fissadir, 'out' );

% load mask 
load( fullfile( roifile) );

% if roi directory doesn't exist, create it
if ~exist( roidir, 'dir' )
   mkdir( roidir );
end
% if output directory doesn't exist, create it
if ~exist( outdir, 'dir' )
   mkdir( outdir );
end

%% Generate ROI files
% roi's in masks, with format szx x szy x cell (e.g. 512 x 512 x 101)

[szX, szY, N] = size(masks); % number of ROIs
emptyImg = zeros( szX, szY );
roizip   = fullfile( fissadir, 'rois.zip' );
roifiles = cell(N,1); % list of roi files
for i=1:N
   mask  = masks(:,:,i); 
   % grow mask by 1 pixel to include entire cell within the mask
   elmt  = strel('disk', 1); 
   mask  = imdilate( mask, elmt );
   fname = sprintf('%d.roi', i); % filename for current mask
   
   % generate ROI files by manually writing binary files 
   [err, roi_fname] = exportMask2ROI( roidir, fname, mask );   
   roifiles{i} = roi_fname;
end

% Generate zip file of all ROIs - note that doing this manually in Finder
% (assuming a mac imagecomputer) generated a zip file that failed when being
% input to the fissa code, but I've no idea why?!
zip( roizip, roifiles );

%% Run FISSA
% Now we need to run fissa using python
pyfun = [pwd '/testFISSA.py'];
python_executable = '/Users/mgo/anaconda3/envs/neuroSEE/bin/python';
pystr =  [python_executable ' ' pyfun];
[status, pyout] = system( pystr );


%    roi_text = fullfile( fissadir, sprintf('roi%d.txt', i) );
%    T = table( x, y );
%    writetable( T, roi_text, 'WriteVariableNames', false );
% Select from 1 to 7 for Roi.LINE, Roi.RECTANGLE, Roi.POINT, Roi.OVAL, 
%                        Roi.POLYLINE, Roi.POLYGON, Roi.ANGLE 
%   OR   polygon=0, rect=1, oval=2, line=3, freeline=4, polyline=5, noRoi=6,
%        freehand=7, traced=8, angle=9, point=10;
%    MIJ.setRoi(roi, 3);
%    MIJ.run('ROI Manager...')




