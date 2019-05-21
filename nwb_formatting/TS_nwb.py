# coding: utf-8
# !/usr/bin/env python3
# -*- coding: utf-8 -*-

# =============================================================================
# Script to build NWB Files
# =============================================================================

# Python Imports
import datetime
# import time
# import h5py
# import numpy as np
# import glob
import os
import pandas as pd
from openpyxl import Workbook, load_workbook, worksheet

# matlab Imports
# import matlab.engine

# NWB imports
from pynwb import NWBFile, get_manager
from pynwb import TimeSeries, behavior
from pynwb import NWBHDF5IO

# =============================================================================
# Reading metadata from Excel spreadsheet
# =============================================================================

def sbj_sheet_list(data_sheet):

    """ Get sheet names from Excel spreadsheet """

    data_sheet = pd.ExcelFile(data_sheet)
    sheets = data_sheet.sheet_names

    return sheets


def get_exp_data(data_sheet):

    """ Get sheet names Experimental Details """

    exp_details = pd.read_excel(data_sheet, sheet_name='Experiment Details', index_col=0)

    return exp_details


def create_mouse_id_folder(mouse_id, nwb_output_path):

    """ Creates folder with Mouse ID """

    mouse_folder_dir = os.path.join(nwb_output_path, mouse_id)
    if not os.path.exists(mouse_folder_dir):
        os.makedirs(mouse_folder_dir)

    return mouse_folder_dir


def create_mouse_date_folder(mouse_id, exp_date, nwb_output_path):

    """ Creates a dated folder within the Mouse folder """

    mouse_target_dir = create_mouse_id_folder(mouse_id, nwb_output_path)
    mouse_date_dir = os.path.join(mouse_target_dir, exp_date)
    if not os.path.exists(mouse_date_dir):
        os.makedirs(mouse_date_dir)

    return mouse_date_dir


def is_sbj_sheet(sheet):

    """ Returns list of sheets titled subject details """

    if 'template' in sheet.lower():
        return False
########    Bit hacky due to sheet naming #########
    if sheet.lower().startswith('subject details') or sheet.lower().startswith('m'):
        return True
    return False


def get_subject_sheet(data_sheet):
    sheets = sbj_sheet_list(data_sheet)
    sbj_details = [x for x in sheets if is_sbj_sheet(x)]

    return sbj_details


def get_times_from_sheet(subject_sheet):

    """ Extracts all times from experiment sheet, then a list of only
    the imaging times and behavioural times """

    times = subject_sheet.loc['Time':].copy()
    times.columns = times.iloc[0]
    times.drop(['Time'], inplace=True)
#   This returns a check on whether there is a TRUE inside the column
#   So contains list of imaging times that have a corresponding x
    imaging_times_list = times[times['2P (X)'].notnull()]

    return times, imaging_times_list


def mouse_folder_details(subject_sheet, data_sheet, updated_imaging_times_list, nwb_output_path):

    """ Building folder name, creating folder if it's not there using
        notes from experimental data sheet """

    mouse_id = subject_sheet.loc['Number', 'Unnamed: 1']
    exp_data = get_exp_data(data_sheet)
    exp_date = exp_data.iloc[4, 0].replace('.', '')

    mouse_date_folder = create_mouse_date_folder(mouse_id, exp_date, nwb_output_path)
    session_start = str(updated_imaging_times_list.index[0])[0:5]
    output_file_name = mouse_id + '_' + exp_date + '_' + session_start.replace(':', '_') + '.nwb'
    output_filepath = os.path.join(mouse_date_folder, output_file_name)

    return mouse_id, exp_date, session_start, mouse_date_folder, output_filepath

# =============================================================================
# Extracting the variables for the NWB File
# =============================================================================


def nwb_file_variables(imaging_folder, subject_sheet, mouse_id, exp_date, session_start):

    """ Required variables to create an NWB File,
        these have to be pre-set:
        source = path to the raw data
        lab = name of your laboratory
        institution = name of your institution
        experiment_description = description of your experiment

        Other variables are pulled from the metadata spreadsheet """

    source = imaging_folder
    session_description = mouse_id + ' ' + exp_date + ' ' + session_start
    identifier = mouse_id + ' ' + exp_date + ' ' + session_start
    session_start_time = session_start
    lab = 'Neural Coding'
    institution = 'Imperial College London'
    experiment_description = 'Two Photon GCamp Imaging of Alzheimers mouse models'
    virus = subject_sheet.index[6] + ':' + str(subject_sheet.loc['GCaMP Inj', 'Unnamed: 1']) + ', ' + \
        subject_sheet.index[7] + ':' + str(subject_sheet.loc['Methoxy Inj', 'Unnamed: 1'])

    return source, session_description, identifier, session_start_time, lab, institution, experiment_description, virus

# =============================================================================
# Managing errors in experimental data sheet
# =============================================================================

def update_2p_spreadsheet(data_sheet, i, imaging_folder, imaging_times_list, all_imaging_times):

    """ Gets a list of the 2p folder imaging times, compares with the times from the spreadsheet
        overwrites the incorrect times in spreadsheet with those from 2p folder list """

    exp_data = get_exp_data(data_sheet)
    exp_date = exp_data.iloc[4, 0].replace('.', '')

    twop_folders = os.listdir(imaging_folder)
    twop_folders = [folder_name for folder_name in twop_folders if exp_date in folder_name]
    all_times_to_string = sorted([time_to_str.strftime('%H:%M') if isinstance(time_to_str, datetime.time) else
                                  time_to_str for time_to_str in all_imaging_times.index])
    folder_times = sorted([folder[9:14].replace('_', ':') for folder in twop_folders])

    if len(folder_times) != len(all_times_to_string):
        print("Spreadsheet times don't match the imaging folders for " + i + ' ' + data_sheet, twop_folders,
              imaging_times_list.index)
    else:
        for folder, img_time in zip(twop_folders, imaging_times_list.index):
            folder_time, sheet_time = folder[9:14].replace('_', ':'), str(img_time)[0:5]
            if folder_time != sheet_time:
                wb = load_workbook(data_sheet, read_only=False, keep_vba=True)
                sheet = wb[i]

                for incorrect_time in sheet['A']:
                    if sheet_time == str(incorrect_time.value)[0:5]:
                        sheet['A' + str(incorrect_time.row)] = folder_time
                        wb.save(data_sheet[:-4]+'xlsm')

    return twop_folders

def run_matlab(imaging_folder, file):
    eng = matlab.engine.start_matlab()
    g_ts, r_ts = eng.neuroSEE_batch(imaging_folder, file, nargout=2)

    return g_ts, r_ts

# Building path to 2P.tif file, using imaging time as name for epoch
def twop_ts(data_sheet, updated_imaging_times_list, nwb_file, twop_folders, imaging_folder, exp_date, tif_paths):
    times_to_string = [time_str.strftime('%H:%M') if isinstance(time_str, datetime.time) else time_str for time_str
                       in updated_imaging_times_list.index]
    updated_imaging_times_list.index = times_to_string
#     tif_paths = get_tif_paths(twop_folders, imaging_folder)
    for folder_time in twop_folders:
# #         Need to check that tif_paths is only 1
        for tif_file in tif_paths:
            epoch_name = folder_time[9:14].replace('_', ':')

#         full_time = folder_time[9:17].replace('_', ':')
#         datetime_object = (datetime.datetime.strptime(full_time, '%H:%M:%S')).time().isoformat()
#         ftr = [3600,60,1] #?
        
##### FAKE DATA ######
            data = list(range(100, 200, 10))
#             delete these after, just using dummy data
            g_ts = data
            r_ts = data

# MATLAB SCRIPT CALLING HERE #
# Matlab script loading and running working, the composite image generation seigfred will incorporate into his script
#             eng = matlab.engine.start_matlab()
#             g_ts, r_ts = eng.neuroSEE_batch(imaging_folder,tif_file,10, 150, nargout=2)
#             eng.quit()
#             print(g_ts, r_ts)

            exp_data = get_exp_data(data_sheet)
            source = 'Licensee: ' + exp_data.at['Licensee', 'Unnamed: 1']

            nwb_file.create_epoch(epoch_name, tif_file, 2.0, 4.0, [epoch_name],
                                  updated_imaging_times_list.loc[epoch_name, 'Imaging details'] + '. ' +
                                  updated_imaging_times_list.loc[epoch_name, 'Remarks'])

            green_ts = TimeSeries(epoch_name + ' ' + 'green_timeseries',  # names the TS in the epoch
                              source,
                              g_ts,
                              'what is the SI unit?',
                              # Or timestamps = list, tuple or dataset for the timeseries
                              starting_time=0.0,  # has to be a float
####### rate needs to be pulled form somewhere #####
                              rate=1.0
                              )

            red_ts = TimeSeries(epoch_name + ' ' + 'red_timeseries',  # names the TS in the epoch
                            source,
                            r_ts,
                            'what is the SI unit?',
                            # Or timestamps = list, tuple or dataset for the timeseries
                            starting_time=1.0,  # has to be a float
####### rate needs to be pulled form somewhere #####
                            rate=1.0
                            )

            nwb_file.add_acquisition(green_ts, epoch_name)
            nwb_file.add_acquisition(red_ts, epoch_name)
            
# generic pipeline wrapper
def get_data(data_sheet, imaging_folder, nwb_output_path):

    sbj_details = get_subject_sheet(data_sheet)

    all_imaging_times = pd.DataFrame()

    for sheet in sbj_details:
        subject_sheet = pd.read_excel(data_sheet, sheet_name=sheet, index_col=0)
        times, imaging_times_list = get_times_from_sheet(subject_sheet)
        data = pd.DataFrame(imaging_times_list)
        all_imaging_times = all_imaging_times.append(data)
        
    for i in sbj_details:
        subject_sheet = pd.read_excel(data_sheet, sheet_name=i, index_col=0)
        times, imaging_times_list = get_times_from_sheet(subject_sheet)
    
        if len(imaging_times_list) > 0:
            twop_folders = update_2p_spreadsheet(data_sheet, i, imaging_folder, imaging_times_list, all_imaging_times)
            tif_paths = []
            for folder_time in twop_folders:
                
                for file in os.listdir(os.path.join(imaging_folder, folder_time)): 
                    if file.endswith('_XYT.tif'):
                        tif_paths.append(os.path.join(imaging_folder, folder_time, folder_time + '_XYT.tif'))
    
            if len(tif_paths) > 0:
                mouse_id, exp_date, session_start, mouse_date_folder, output_filepath = \
                    mouse_folder_details(subject_sheet, data_sheet, imaging_times_list, nwb_output_path)

                source, session_description, identifier, session_start_time, lab, institution, experiment_description, \
                    virus = nwb_file_variables(imaging_folder,subject_sheet, mouse_id, exp_date, session_start)

                nwb_file = NWBFile(
                            source=source,
                            session_description=session_description,
                            identifier=identifier,
                            session_start_time=session_start_time,
                            lab=lab,
                            institution=institution,
                            experiment_description=experiment_description,
                            virus=virus
                            )

                tifs = twop_ts(data_sheet, imaging_times_list, nwb_file, twop_folders, imaging_folder,
                                    exp_date, tif_paths)

                file = NWBHDF5IO(output_filepath, manager=get_manager(), mode='w')
                file.write(nwb_file)
                file.close()
        