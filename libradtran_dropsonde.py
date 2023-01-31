# Author: Joshua Mueller
# Date: 2023-01-22
# Description: This script is used to test the libradtran library
# It is used to calculate the brightness temperature of a given channel at a given time, location, and altitude


wavelength_ranges = {
    0: "7700 12000",
    1: "8095 9195",
    2: "10350 11130",
    3: "10847 12000",
    4: "11500 12000"
    }


import re
import numpy as np
import os
import xarray as xr 
import glob
import subprocess
import time
from tqdm import tqdm
import concurrent.futures

def find_file_with_key(directory, key):
    try:
        # This for loop iterates through all the entries in the directory passed to the function
        for entry in os.scandir(directory):
            # If the current entry is a file and the key is present in the file name
            if entry.is_file() and key in entry.name:
                # Return the full path of the file
                return os.path.join(directory, entry.name)

    except FileNotFoundError:
    # If the directory passed to the function does not exist, return None
        return None

def datetime64_to_seconds(timestamps):
    # Subtract the datetime of the Unix epoch from the input timestamps
    # The result is an array of timedeltas
    delta = timestamps - np.datetime64('1970-01-01T00:00:00')
    
    # Divide the timedeltas by the timedelta representing one second
    # The result is an array of seconds
    seconds = delta / np.timedelta64(1, 's')
    
    # Take the modulo of the seconds with the number of seconds in a day
    # This is done to ensure that the seconds are within a 24 hour period
    return seconds % (24 * 3600)

def closest_file(seconds, file_seconds, sorted_files):
    closest_files = []
    # Iterate through the seconds passed in the function
    for i in seconds:
        # Calculate the absolute difference between the input seconds and the seconds associated with the files
        diff = np.abs(np.array(file_seconds) - i)
        # Find the index of the file with the smallest difference
        idx = np.argmin(diff)
        # Get the file at the index found above
        closest_file = sorted_files[idx]
        # Append the file to the list of closest files
        closest_files.append(closest_file)
    return closest_files

def filter_files(file_list, M, D):
    pattern = f'{M}{D}' +'_([0-9]{5})SOD.dat'
    matched_files = [(int(re.match(pattern, file).group(1)), file) for file in file_list if re.match(pattern, file)]
    matched_files.sort()
    sorted_files = [file for _, file in matched_files]
    digits = [int(re.match(pattern, file).group(1)) for _, file in matched_files]
    return (sorted_files, digits)

def start_libradtran(sod, ts, lat, lon, alt, atmos_file, cloud_fraction_file, ozone, date, channel):
    # Define the paths for the libradtran library and the home directory
    lib_file_path = '/projekt_agmwend/data/EUREC4A/11_VELOX-Tools/libradtran'
    home_path = '/home/jomueller/sim_test/20200202_CF'
    
    # Define a dictionary of wavelength ranges for different channels
    # Get the wavelength range for the specified channel
    wavelength_range = wavelength_ranges.get(channel, "")
    # Print an error message and return if the channel is invalid
    if not wavelength_range:
        print(f"Error: Invalid channel {channel}")
        return
    
    # Create the input file with the specified parameters
    input_file = f'{home_path}/VELOX_{channel}_TB_{date}.inp'
    output_file = f'{home_path}/VELOX_{channel}_TB_{date}_{sod}.out'
    with open(input_file, 'w') as f:
        f.write('data_files_path /opt/libradtran/2.0.4/share/libRadtran/data\n')
        f.write('rte_solver disort\n')
        f.write('mol_abs_param reptran medium\n')
        f.write('atmosphere_file /projekt_agmwend/data/EUREC4A/11_VELOX-Tools/add_data/afglt_evi.dat\n')
        f.write(f'radiosonde {atmos_file} H2O RH\n')
        f.write(f'cloud_fraction_file {cloud_fraction_file} \n')
        f.write(f'time {ts}\n')
        f.write(f'latitude {"S" if lat < 0 else "N"} {abs(lat)}\n')
        f.write(f'longitude {"W" if lon < 0 else "E"} {abs(lon)}\n')
        f.write(f'zout {alt}\n')
        f.write('phi 0.0\n')
        f.write('umu 1.0\n')
        f.write(f'wavelength {wavelength_range}\n')
        f.write('source thermal\n')
        f.write('output_process per_nm\n')
        f.write('output_process integrate\n')
        f.write('output_quantity brightness\n')
        f.write('output_user lambda sza uu\n')
        f.write('albedo_library IGBP\n')
        f.write('brdf_rpv_type 17\n')
        f.write('sur_temperature 300\n')
        f.write(f'mol_modify O3 {ozone}  DU\n')
        f.write('aerosol_default\n')
    
    # Start the uvspec program using the subprocess library
    process= subprocess.Popen(['/opt/libradtran/2.0.4/bin/uvspec'], stdin=open(input_file, 'r'), stdout=open(output_file, 'w'), stderr=subprocess.STDOUT)
    # Check the status of the process every 0.1 seconds
    while process.poll() is None:
        time.sleep(0.1)
    
    return output_file

# Define the flight number and the date
f = 'Flight_20200202a'
fnum = str(f[7:-1])

Y, M, D = fnum[0:4], fnum[4:6], fnum[6:]

# Define the path to the raw data and the processed data
data_raw_path = f'/projekt_agmwend2/data_raw/EUREC4A_raw_only/06_Flights/{f}'
data_path     = f'/projekt_agmwend/data/EUREC4A/06_Flights/{f}'

# Define the channels used for the simulation; the channels are defined as follows:    
# 0:  7700nm - 12000nm, 
# 1:  8095nm -  9195nm,
# 2: 10350nm - 11130nm,
# 3: 10847nm - 12000nm,
# 4: 11500nm - 12000nm, 

channels = [0, 1, 2, 3, 4]

# Define the path to the atmosphere files
atmos_path = f'/projekt_agmwend/data/EUREC4A/03_Soundings/RS_for_libradtran/Merged/{M}{D}/'
#atmos_path = f'/home/jomueller/sim_test/atmos_dir/'

# Define the path to the cloud fraction files
cloud_fraction_path = f'/home/jomueller/sim_test/CF_dir/'

# Define the path to the libradtran files
lib_file_path = '/projekt_agmwend/data/EUREC4A/11_VELOX-Tools/libradtran/'

# Find the parameter file using the find_file_with_key function
params_path = find_file_with_key(data_path+'/BAHAMAS', '.nc')

# Define the ozone layer thickness and the temporal resolution
ozone = 375
temporal_resolution = 60

# Open the parameter file using xarray
xrparams = xr.open_dataset(params_path)

channel_wavelength_string = '\n'.join([f'Channel {i} : {wavelength_ranges.get(i, "")}[nm] \n' for i in range(5)]),

    
# Create an xarray dataset for the navigation data; this is necessary for the resampling
# Resampling is used here to define the temporal resolution of the output data
xrHALO = xr.Dataset(
    data_vars=dict(
        lat=(["time"], xrparams['IRS_LAT'].values),
        lon=(["time"], xrparams['IRS_LON'].values),
        alt=(["time"], xrparams['IRS_ALT'].values),
        roll=(["time"], xrparams['IRS_PHI'].values),
        pitch=(["time"], xrparams['IRS_THE'].values),
        hdg=(["time"], xrparams['IRS_HDG'].values),
        gs=(["time"], xrparams['IRS_GS'].values),

    ),
    coords=dict(

        time=xrparams['TIME'].values,
    ),
)

# Resample the xrHALO dataset with the temporal_resolution specified in seconds.
# This will make sure that data is only taken every temporal_resolution seconds
nav_resample = xrHALO.resample(time=f'{temporal_resolution}S').mean() 

# Convert datetime64 values in `nav_resample` to seconds of day
array_sod = datetime64_to_seconds(nav_resample.time.values)

# Create an array of timestamps in the format 'YYYY MM DD HH MM SS'
ts = nav_resample.time.dt.strftime('%Y %m %d %H %M %S').values

# Extract the latitude, longitude, and altitude values from the resampled dataset
array_lat = nav_resample.lat.values
array_lon = nav_resample.lon.values
array_alt = nav_resample.alt.values

# Filter the files in the `atmos_path` directory and sort them by seconds of day
array_sorted_atmos_files, array_seconds_of_day_atmos_files = filter_files(os.listdir(atmos_path), M=M, D=D)
array_sorted_cloud_fraction_files, array_seconds_of_day_cloud_fraction_files = filter_files(os.listdir(cloud_fraction_path), M=M, D=D)


# Find the closest atmospheric file to each timestamp in `array_sod`
array_matched_atmos_files = closest_file(array_sod, array_seconds_of_day_atmos_files, array_sorted_atmos_files)

# Find the closest cloud fraction file to each timestamp in `array_sod`

array_matched_cloud_fraction_files = closest_file(array_sod, array_seconds_of_day_cloud_fraction_files, array_sorted_cloud_fraction_files)

# Create an array to store the calculated brightness temperature values
BT = np.zeros((len(array_sod), len(channels)))


#from concurrent.futures import ThreadPoolExecutor

atmos_files = [atmos_path+array_matched_atmos_files[i] for i in range(len(array_sod))]

cloud_fraction_files = [cloud_fraction_path+array_matched_cloud_fraction_files[i] for i in range(len(array_sod))]

futures = []

# Create a thread pool with a certain number of worker threads
with concurrent.futures.ProcessPoolExecutor(max_workers=20) as executor:
    for i in tqdm(range(0, len(array_sod))):
        for channel in channels:
                # Schedule the libradtran simulation for the current time and location, using the closest atmospheric file
                # and store the returned future object in the list
            future = executor.submit(start_libradtran, array_sod[i], ts[i], array_lat[i], array_lon[i], array_alt[i]/1e3, atmos_files[i], cloud_fraction_files[i], ozone, fnum, channel)
            futures.append(future)

    for future in tqdm(concurrent.futures.as_completed(futures), total=len(futures)):

        output_file = future.result()
        with open(output_file, "r") as output:

            string_BT = output.readlines()[-1].split()[-1]
            #print(string_BT)
            try:
                float_BT = float(string_BT)
            except ValueError:
                float_BT = np.nan
            channel = int(output_file[-25:-24])
            time    = int(output_file[-11:-6])
            
        match = (array_sod == time)
        BT[match, channel] = float_BT


xr.Dataset(
    data_vars=dict(
        BT =(["time", "channel"], BT)
    ),
    coords=dict(
        time=nav_resample.time.values,
        channel=channels,
    ),
    attrs=dict(
        
        wavelength_range = channel_wavelength_string,
        description = "Brightness temperature calculated with libradtran",
        units = "K",
        data_files_path = "/opt/libradtran/2.0.4/share/libRadtran/data",
        rte_solver = "disort",
        mol_abs_param = "reptran medium",
        atmosphere_file = "/projekt_agmwend/data/EUREC4A/11_VELOX-Tools/add_data/afglt_evi.dat",
        atmosphere_description = "AFGL (1986) thermally coupled troposphere and stratosphere",
        radiosonde = "JOANNE Level 3 Dropsonde data from EUREC4A campaign",
        zout = "HALO flight altitude",
        phi = 0.0,
        umu = 1.0,
        source = "thermal",
        output_process = ["per_nm", "integrate"],
        output_quantity = "brightness",
        output_user = "lambda sza uu",
        albedo_library = "IGBP",
        brdf_rpv_type = 17,
        sur_temperature = 300,
        mol_modify = f'O3 {ozone} DU',        
        ),

    ).to_netcdf(f'/home/jomueller/sim_test/VELOX_{fnum}_{temporal_resolution}S_TB.nc', mode='w')

