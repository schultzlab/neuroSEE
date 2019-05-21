# NWB formatting integration

The scripts here will extract the data from the .ini file, the word document and 
TimeSeries matlab data then convert this into the NWB format

# Usage

## Automated 
To create a `NWB` file, use the `get_data()` function found in `TS_nwb.py`. You will need
to provide: 

- the path to the Excel template with metadata:
    ```python
    metadata_spreadsheet = '/path/to/xlsm/file'
    ```

- the path to a directory holding imaging data from different sessions:
    ```python
    image_directory = '/path/to/images'
    ```
    Currently, this directory is assumed to have the following structure:
    ```
    image_directory  
    │
    └───YYYYMMDD_HH_MM_SS_2P
    │   │   YYYYMMDD_HH_MM_SS_2P_XYTZ.tiff
    │   │   ...
    │   
    └───YYYYMMDD_HH_MM_SS_2P
        │   YYYYMMDD_HH_MM_SS_2P_XYTZ.tiff
        │   ...
    ```
    where each subdirectory is named as `YYYYMMDD_HH_MM_SS_2P` where `YYYY`,
    `MM`, `DD`, `HH`, `MM`, `SS` are the year, month, day, hour, minute and second of the recording
    respectively. The year, month, day, hour and minutes have to correspond with those given in `metadata_spreadsheet`
    
    Each subdirectory may or may not contain a `tiff` file with a filename ending in `_XYTZ.tiff`.

- an output directory:
    ```python
    output_directory = '/path/to/output/directory'
    ```
    This is where the `NWB` file will be saved to.
    
Please note that the actual images are not saved in the `NWB` file.