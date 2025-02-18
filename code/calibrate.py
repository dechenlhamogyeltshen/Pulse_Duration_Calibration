import numpy as np
import os
import astropy.units as u
import matplotlib.pyplot as plt
import pandas as pd
from datetime import datetime, timedelta
from astropy.time import Time, TimeDelta
import huxt as H
import huxt_analysis as HA

# --- Load data ---------------------------------------------------------------

# Read in Blair's pairing on DONKI and CR2003
project_dirs = H._setup_dirs_()
crpath = os.path.join(project_dirs['input'],'(I)CMEs_SingleEvents.csv')

# Load the CSV file into a DataFrame
crlist = pd.read_csv(crpath)

# Function to convert a date string to a datetime object
def convert_to_datetime(date_str):
    # List of possible date formats
    date_formats = [
        '%Y-%m-%d %H:%M:%S',
        '%Y/%m/%d %H%M',
        '%Y/%m/%d %H%M(%S)',
        '%Y/%m/%d %H%M(S)',
        '%Y/%m/%d %H:%M',
        '%Y/%m/%d %H:%M:%S',
        '%Y/%m/%d %H:%M(S)',
        '%d/%m/%Y %H:%M'
    ]
    for date_format in date_formats:
        try:
            return datetime.strptime(date_str, date_format)
        except ValueError:
            continue
    # If no format matches, return None
    return None

# Convert date columns to datetime objects using datetime library
date_columns = ['CME_Time', 'Time_21.5', 'ICME_Start_Time', 'ICME_End_Time', 'Disturbance_Time']
for column in date_columns:
    crlist[column] = crlist[column].apply(lambda x: convert_to_datetime(x) if isinstance(x, str) else None)

# Function to convert Time_Error from string to timedelta
def convert_to_timedelta(time_str):
    try:
        return pd.to_timedelta(time_str)
    except ValueError:
        return None

# Convert Time_Error column from string to timedelta
# crlist['Time_Error'] = crlist['Time_Error'].apply(lambda x: convert_to_timedelta(x) if isinstance(x, str) else None)

# Convert Time_Error column from timedelta to days
# def convert_timedelta_to_days(td):
    # return td.days + td.seconds / 86400 if isinstance(td, timedelta) else None

# crlist['Time_Error_Days'] = crlist['Time_Error'].apply(convert_timedelta_to_days)

# compute 21.5-215 transit time
crlist['tt_21'] = np.nan 
for irow in range(0, len(crlist)):
    crlist.loc[irow,'tt_21'] = crlist.loc[irow,'Disturbance_Time'] - crlist.loc[irow,'Time_21.5']
# Convert  from timedelta to days
crlist['tt_21'] = crlist['tt_21'].apply(lambda x: x.days + x.seconds / 86400 if isinstance(x, timedelta) else None)

crlist['duration_1au'] = np.nan 
for irow in range(0, len(crlist)):
    crlist.loc[irow,'duration_1au'] =  crlist.loc[irow,'ICME_End_Time'] - crlist.loc[irow,'ICME_Start_Time']
# Convert  from timedelta to days
crlist['duration_1au'] = crlist['duration_1au'].apply(lambda x: x.days + x.seconds / 86400 if isinstance(x, timedelta) else None)

# Drop rows where the absolute value of Time_Error_Days is greater than 2 hours
# crlist = crlist[abs(crlist['Time_Error_Days']) <= 0.08]
# crlist = crlist.reset_index()

# print('N = ' + str(len(crlist)))

# Drop rows where angular radius is shorter than angular distance
crlist['angdist'] = np.sqrt(crlist['lat']*crlist['lat'] + crlist['lon']*crlist['lon'])
crlist['Y'] = crlist['angdist']/crlist['Ang_rad']
#crlist = crlist[abs(crlist['Y']) <= 1.0]

print('N = ' + str(len(crlist)))

#===============================================================================

# --- Initialize Model --------------------------------------------------------

r_in = 21.5 * u.solRad
vsw = 350
durations = np.arange(5, 23, 0.1)  # CME durations in hours

# Background wind configuration
v_boundary = np.ones(128) * vsw * (u.km / u.s)

# Initialize the model
model = H.HUXt(
    v_boundary=v_boundary,
    latitude=0 * u.deg,
    r_min=r_in,
    simtime=10 * u.day,
    dt_scale=4,
    lon_out=0 * u.deg
)

# Initialize Data Storage
ang_width = []
transit_time = []
arrival_speed = []

# First row: Observed values
ang_width.append(['Angular Width'] + list(crlist['Ang_rad']))
transit_time.append(['Observed'] + list(crlist['tt_21']))
arrival_speed.append(['Observed'] + list(crlist['V_ICME']))

# Spheroidal
spheroidal_tt = []
spheroidal_as = []

# Iterate over all CMEs in crlist
for _, onecme in crlist.iterrows():
    spheroidal_cme = H.ConeCME(
        t_launch=0 * u.day,
        longitude=0.0 * u.deg,
        latitude=0.0 * u.deg,
        initial_height=r_in,
        width=onecme['Ang_rad'] * u.deg,
        v=onecme['V'] * (u.km / u.s),
        thickness=0 * u.solRad,
        cme_fixed_duration=False
    )

    model.solve([spheroidal_cme])

    # Compute the transit time
    stats = model.cmes[0].compute_arrival_at_body('EARTH')
    s_tt = stats['t_transit'].value
    arrival_time = stats['t_arrive']

    # Find the arrival speed within 1 day of arrival
    earth_series = HA.get_observer_timeseries(model, observer='Earth')
    mask = (Time(earth_series['time']) >= arrival_time) & (
            Time(earth_series['time']) <= arrival_time + 3 * u.day)
    s_v_1au = earth_series.loc[mask, 'vsw'].max()

    spheroidal_tt.append(s_tt)
    spheroidal_as.append(s_v_1au)

transit_time.append(['Spheroidal'] + list(spheroidal_tt))
arrival_speed.append(['Spheroidal'] + list(spheroidal_as))

for duration in durations:
    tt_row = [f"tt_{duration}h"]
    v_row = [f"v_{duration}h"]
    
    # Iterate over all CMEs in crlist
    for _, onecme in crlist.iterrows():
        cme = H.ConeCME(
            t_launch=0 * u.day,
            longitude=0.0 * u.deg,
            latitude=0.0 * u.deg,
            initial_height=r_in,
            width=onecme['Ang_rad'] * u.deg,
            v=onecme['V'] * (u.km / u.s),
            thickness=0 * u.solRad,
            cme_fixed_duration=True,
            fixed_duration=duration * 60 * 60 * u.s
        )

        model.solve([cme])

        # Compute the transit time
        stats = model.cmes[0].compute_arrival_at_body('EARTH')
        tt_val = stats['t_transit'].value
        arrival_time = stats['t_arrive']

        # Find the arrival speed within 1 day of arrival
        earth_series = HA.get_observer_timeseries(model, observer='Earth')
        mask = (Time(earth_series['time']) >= arrival_time) & (
                Time(earth_series['time']) <= arrival_time + 3 * u.day)
        v_1au = earth_series.loc[mask, 'vsw'].max()

        tt_row.append(tt_val)
        v_row.append(v_1au)

    transit_time.append(tt_row)
    arrival_speed.append(v_row)

# Save results
data_dir = project_dirs['output']
out_path = os.path.join(data_dir)

pd.DataFrame(ang_width).to_csv(os.path.join(out_path,"angular_width.csv"), index=False, header=False)
pd.DataFrame(transit_time).to_csv(os.path.join(out_path,"cme_transit_time.csv"), index=False, header=False)
pd.DataFrame(arrival_speed).to_csv(os.path.join(out_path,"cme_arrival_speed.csv"), index=False, header=False)
