#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  7 17:13:55 2019
Edited on Thu Aug   5          2022

@author: cedavey, sihaolu, anngo

Python script to run Fissa through Python command line. Script takes three inputs from the command line in order:
 - tiff : can be directory containing tif stacks/s, can be tif filename, or list of tif filenames
 - roizip : zip file of ROIs generated through 'generateImageJROIfiles.m'
 - outdir : directory to save the Fissa output into

 Example:
 ```
 $ python runFISSA.py FISSAtest/ FISSAtest/rois.zip FISSAtest/out/
 $ python runFISSA.py folder1/tiffile1.tif,folder2/tiffile2.tif,folder3/tiffile3 FISSAtest/rois.zip FISSAtest/out/
 ```

"""

import os
import sys
import fissa

roizip = sys.argv[2]
outdir = sys.argv[3]

tiff = list(sys.argv[1].split(","))

# generate an experiment object
# Ann: I edited fissa.core to include deltaf_raw and deltaf_result in the matlab outputs as they weren't being saved
# Edit fissa.core in .../anaconda3/envs/neuroSEE/lib/python3.7/site-packages/fissa/core.py
experiment = fissa.Experiment(tiff, roizip, outdir)
# separate neuropil
experiment.separate(redo_sep=True, redo_prep=True)
# calculate df/f0
experiment.calc_deltaf(freq=30.9)
# save to matlab
experiment.save_to_matlab()


