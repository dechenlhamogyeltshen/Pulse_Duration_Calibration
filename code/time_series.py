import numpy as np
import os
import astropy.units as u
import pandas as pd
from datetime import datetime, timedelta
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

# Iterate over all CMEs in crlist
for _, onecme in crlist.iterrows():
    
    cr = onecme['cr_num'] # March, so Earth is around  - 7 degrees lat
    lat = onecme['earth_lat']

    vr_in = Hin.get_MAS_long_profile(cr,lat)
    vr_21 = Hin.map_v_boundary_inwards(vr_in, 30*u.solRad, 21.5*u.solRad)

    #now run HUXt
    model = H.HUXt(v_boundary=vr_21,
                   cr_num=cr,
                   cr_lon_init=onecme['cr_lon_init'],
                   latitude=lat,
                   simtime=10*u.day,
                   dt_scale=4,
                   r_min=21.5*u.solRad,
                   frame = 'sidereal')
    
    cme = H.ConeCME(t_launch=0 * u.day, 
                    longitude=onecme['lon'] * u.deg, 
                    latitude=onecme['lat'] * u.deg, 
                    initial_height=21.5*u.solRad, 
                    width=2.0 * onecme['Ang_rad'] * u.deg, 
                    v=onecme['V'] * (u.km / u.s), 
                    thickness=0.0 * u.solRad,
                    cme_fixed_duration=False)
    
    model.solve([cme]) 

    #the Earth time series can be plotted, along with OMNI data (downloaded on demand),using:
    fig, axs = HA.plot_earth_timeseries(model, plot_omni = True)
    
    data_dir = project_dirs['HUXt_figures']
    out_path = os.path.join(data_dir, "time_series")

    filename = f"{onecme['cr_num']}_{onecme['cr_lon_init']}.pdf"
    filepath = os.path.join(out_path, filename)

    fig.savefig(filepath, bbox_inches='tight')
