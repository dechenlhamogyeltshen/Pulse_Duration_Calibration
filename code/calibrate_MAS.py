import numpy as np
import os
import astropy.units as u
import pandas as pd
from datetime import datetime, timedelta
from astropy.time import Time, TimeDelta
import huxt as H
import huxt_inputs as Hin
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

def get_earth_lat(dt):
    """
    A function to return Earth latitude for a given date, in radians
 
    Parameters
    ----------
    dt : datetime
 
    Returns
    -------
    E_lat: Earth latitude, with astropy units of radians
 
    """
    cr, cr_lon_init = Hin.datetime2huxtinputs(dt)
    # Use the HUXt ephemeris data to get Earth lat over the CR
    # ========================================================
    dummymodel = H.HUXt(v_boundary=np.ones(128)*400*(u.km/u.s), simtime=0.1*u.day, 
                         cr_num=cr,cr_lon_init=cr_lon_init, lon_out=0.0*u.deg)
    # retrieve a bodies position at each model timestep:
    earth = dummymodel.get_observer('earth')
    # get average Earth lat
    E_lat = np.nanmean(earth.lat_c)
    E_lat = E_lat.to(u.deg)  # Convert to degrees
    return E_lat

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

# Add a new column with carrington rotation number of CME time
crlist['cr_num'] = crlist['CME_Time'].apply(lambda dt: Hin.datetime2huxtinputs(dt)[0])
crlist['cr_lon_init'] = crlist['CME_Time'].apply(lambda dt: Hin.datetime2huxtinputs(dt)[1])
crlist['earth_lat'] = crlist['CME_Time'].apply(lambda dt: get_earth_lat(dt))

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

# --- Spheroidal --------------------------------------------------------

# Initialize Data Storage
transit_time = []
arrival_speed = []

# First three rows
transit_time.append(['Angular_Width'] + list(crlist['Ang_rad']))
arrival_speed.append(['Angular_Width'] + list(crlist['Ang_rad']))

transit_time.append(['Velocity'] + list(crlist['V']))
arrival_speed.append(['Velocity'] + list(crlist['V']))

transit_time.append(['Observed'] + list(crlist['tt_21']))
arrival_speed.append(['Observed'] + list(crlist['V_ICME']))

r_in = 21.5 * u.solRad # initial height
durations = np.arange(0.1, 20.0, 0.1)  # CME durations in hours

def spheroidal():
    ''' Solve for spheroidal CME transit times'''
    
    spheroidal_tt = []
    spheroidal_as = []
    
    # Iterate over all CMEs in crlist
    for _, onecme in crlist.iterrows():
        
        # Setup HUXt for a standard 30Rs run
        vr_in = Hin.get_MAS_long_profile(onecme['cr_num'], onecme['earth_lat'])

        #  Map the inner boundary MAS values inwards from 30 rS to 21.5 rS
        vr_21 = Hin.map_v_boundary_inwards(vr_in, 30*u.solRad, r_in)

        #  Now setup HUXt to run from 10Rs
        model21 = H.HUXt(v_boundary = vr_21, 
                         cr_num=onecme['cr_num'],
                         cr_lon_init=onecme['cr_lon_init'],
                         simtime=10*u.day, 
                         latitude=0*u.deg,
                         lon_out=0*u.deg,
                         dt_scale=4,
                         r_min=r_in)

        cme = H.ConeCME(t_launch=0 * u.day, 
                        longitude=onecme['lon'] * u.deg, 
                        latitude=onecme['lat'] * u.deg, 
                        initial_height=r_in, 
                        width=2.0 * onecme['Ang_rad'] * u.deg, 
                        v=onecme['V'] * (u.km / u.s), 
                        thickness=0.0 * u.solRad,
                        cme_fixed_duration=False)

        model21.solve([cme])

        # Compute the transit time
        stats = model21.cmes[0].compute_arrival_at_body('EARTH')
        tt = stats['t_transit'].value
        arrival_time = stats['t_arrive']

        # Find the arrival speed within 1 day of arrival
        earth_series = HA.get_observer_timeseries(model21, observer='Earth')
        mask = (Time(earth_series['time']) >= arrival_time) & (
                Time(earth_series['time']) <= arrival_time + 3 * u.day)
        v_1au = earth_series.loc[mask, 'vsw'].max()

        spheroidal_tt.append(tt)
        spheroidal_as.append(v_1au)
    
    return spheroidal_tt,spheroidal_as

spheroidal_tt, spheroidal_as = spheroidal()

transit_time.append(['Spheroidal'] + list(spheroidal_tt))
arrival_speed.append(['Spheroidal'] + list(spheroidal_as))

# --- Fixed pulse duration --------------------------------------------------------

for duration in durations:
    tt_row = [f"tt_{duration}h"]
    v_row = [f"v_{duration}h"]
    
    # Iterate over all CMEs in crlist
    for _, onecme in crlist.iterrows():
        
        # Setup HUXt for a standard 30Rs run
        vr_in = Hin.get_MAS_long_profile(onecme['cr_num'], onecme['earth_lat'])

        #  Map the inner boundary MAS values inwards from 30 rS to 21.5 rS
        vr_21 = Hin.map_v_boundary_inwards(vr_in, 30*u.solRad, r_in)

        #  Now setup HUXt to run from 10Rs
        model21 = H.HUXt(v_boundary = vr_21, 
                         cr_num=onecme['cr_num'],
                         cr_lon_init=onecme['cr_lon_init'],
                         simtime=10*u.day, 
                         latitude=0*u.deg, 
                         lon_out=0*u.deg,
                         dt_scale=4, 
                         r_min=r_in)
        
        cme = H.ConeCME(
            t_launch=0 * u.day,
            longitude=onecme['lon'] * u.deg,
            latitude=onecme['lat'] * u.deg,
            initial_height=r_in,
            width=2.0 * onecme['Ang_rad'] * u.deg,
            v=onecme['V'] * (u.km / u.s),
            thickness=0.0 * u.solRad,
            cme_fixed_duration=True,
            fixed_duration=duration * 60 * 60 * u.s
        )

        model21.solve([cme])

        # Compute the transit time
        stats = model21.cmes[0].compute_arrival_at_body('EARTH')
        tt_val = stats['t_transit'].value
        arrival_time = stats['t_arrive']

        # Find the arrival speed within 1 day of arrival
        earth_series = HA.get_observer_timeseries(model21, observer='Earth')
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

pd.DataFrame(transit_time).to_csv(os.path.join(out_path,"cme_transit_time.csv"), index=False, header=False)
pd.DataFrame(arrival_speed).to_csv(os.path.join(out_path,"cme_arrival_speed.csv"), index=False, header=False)
