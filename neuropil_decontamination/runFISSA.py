#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  7 17:13:55 2019

@author: cedavey, sihaolu

Python script to run Fissa through Python. Script takes three inputs from the command line in order:
 - tifdir : directory containing tif stack
 - roizip : zip file of ROIs generated through `generateImageJROIfiles.m
 - outdir : directory to save the Fissa output into

 Example:
 ```
 $ python runFISSA.py FISSAtest/ FISSAtest/rois.zip FISSAtest/out/
 ```

"""

import os
import sys
import fissa
import scipy.io

# from runFISSA import runFISSA
   
# tifdir   = '/Volumes/thefarm2/live/CrazyEights/AD_2PCa/Data/20181016/Processed/20181016_09_49_06/FISSAtest/'
# roizip   = '/Volumes/thefarm2/live/CrazyEights/AD_2PCa/Data/20181016/Processed/20181016_09_49_06/FISSAtest/rois.zip'
# outdir   = '/Volumes/thefarm2/live/CrazyEights/AD_2PCa/Data/20181016/Processed/20181016_09_49_06/FISSAtest/out'

tiffile = [sys.argv[1]]
roizip = sys.argv[2]
outdir = sys.argv[3]



# generate an experiment object
experiment = fissa.Experiment(tiffile, roizip, outdir)
# separate neuropil
experiment.separate(redo_sep=True, redo_prep=True)
# calculate df/f0
experiment.calc_deltaf(freq=30.9)
file_path_raw = outdir + 'deltaf_raw.mat'
file_path_result = outdir + 'deltaf_result.mat'
scipy.io.savemat(file_path_raw, {'deltaf_raw': experiment.deltaf_raw})
scipy.io.savemat(file_path_result, {'deltaf_result': experiment.deltaf_result})
# save to matlab
experiment.save_to_matlab()


