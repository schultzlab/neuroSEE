#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  7 17:13:55 2019

@author: cedavey
"""

import os
import sys

# sys.path.append('/Users/cedavey/Dropbox/code/hippocampus/neuroSEE/neuropil_decontamination/')
# sys.path.append('/Users/cedavey/Dropbox/code/hippocampus/neuroSEE/neuropil_decontamination/fissa-master')
 
import fissa
# from runFISSA import runFISSA
   
tifdir   = '/home/sihao/KozlovBox/Data/Fissa_test/test_files'
roizip   = '/home/sihao/Downloads/testFissa/testFISSA/fissa/rois.zip'
outdir   = '/home/sihao/Downloads/testFissa/testFISSA/fissa/out'

# generate an experiment object
experiment = fissa.Experiment(tifdir, roizip, outdir)
# separate neuropil
experiment.separate(redo_sep=True, redo_prep=True)
# save to matlab
experiment.save_to_matlab()


