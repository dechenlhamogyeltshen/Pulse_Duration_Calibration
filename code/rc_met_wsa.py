import numpy as np
import os
import astropy.units as u
import pandas as pd
from datetime import datetime, timedelta
from astropy.time import Time
import re
from sunpy.coordinates import sun
import helio_coords as hcoords

# for constructing asw based on in-situ omni
import joblib
import onnxruntime as ort

# HUXt
import huxt as H
import huxt_inputs as Hin
import huxt_analysis as HA

#===============================================================================
# <codecell> Import ICME-CME data

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
    
# Function to convert Time_Error from string to timedelta
def convert_to_timedelta(time_str):
    try:
        return pd.to_timedelta(time_str)
    except ValueError:
        return None

# Function to get earth's latitude
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
    
# --- Load data ----------------------------------------------------------------
# Read in Blair's pairing on DONKI and CR2003
project_dirs = H._setup_dirs_()
crpath = os.path.join(project_dirs['input'],'rc_met.csv')

# Load the CSV file into a DataFrame
crlist = pd.read_csv(crpath)

# Convert date columns to datetime objects using datetime library
date_columns = ['CME_Time', 'Time_21.5', 'Time_met', 'Disturbance_Time', 'ICME_Start_Time']
for column in date_columns:
    crlist[column] = crlist[column].apply(lambda x: convert_to_datetime(x) if isinstance(x, str) else None)

# compute 21.5-215 transit time
crlist['tt_21'] = np.nan 
for irow in range(0, len(crlist)):
    crlist.loc[irow,'tt_21'] = crlist.loc[irow,'Disturbance_Time'] - crlist.loc[irow,'Time_21.5']
# Convert  from timedelta to days
crlist['tt_21'] = crlist['tt_21'].apply(lambda x: x.days + x.seconds / 86400 if isinstance(x, timedelta) else None)

# Add a new column with carrington rotation number of CME time
crlist['cr_num'] = crlist['Time_21.5'].apply(lambda dt: Hin.datetime2huxtinputs(dt)[0])
crlist['cr_lon_init'] = crlist['Time_21.5'].apply(lambda dt: Hin.datetime2huxtinputs(dt)[1])
crlist['earth_lat'] = crlist['Time_21.5'].apply(lambda dt: get_earth_lat(dt))

#===============================================================================

# <codecell> Functions for spheroidal and fixed duration ConeCMEs
def spheroidal_rc(onecme,model):
    '''Solve HUXt using a spheroidal cone CME and Richardson-Cane cone CME parameters'''
    
    # CME initialisation date
    t_cme = Time(onecme['Time_21.5'])
    # CME initialisation relative to model initialisation, in days
    dt_cme = (t_cme - model.time_init).jd * u.day
    
    cme = H.ConeCME(t_launch=dt_cme,
                    longitude=onecme['lon'] * u.deg,
                    latitude=onecme['lat'] * u.deg,
                    initial_height=rmin,
                    width=2.0 * onecme['ang'] * u.deg,
                    v=onecme['V']* (u.km / u.s),
                    thickness=0.0 * u.solRad,
                    cme_expansion=False,
                    cme_fixed_duration=False)

    model.solve([cme])
    cme_out = model.cmes[0]

    # Compute the transit time
    stats = cme_out.compute_arrival_at_body('EARTH')
    tt = stats['t_transit'].value
    v_1au = stats['v'].value
    
    return tt, v_1au
    
def spheroidal_met(onecme,model):
    '''Solve HUXt using a spheroidal cone CME and Met Office cone CME parameters'''
    
    # CME initialisation date
    t_cme = onecme['Time_met']
    # CME initialisation relative to model initialisation, in days
    dt_cme = (t_cme - model.time_init).jd * u.day
    
    cme = H.ConeCME(t_launch=dt_cme,
                    longitude=onecme['lon_met'] * u.deg,
                    latitude=onecme['lat_met'] * u.deg,
                    initial_height=rmin,
                    width=2.0 * onecme['ang_met'] * u.deg,
                    v=onecme['V_met']* (u.km / u.s),
                    thickness=0.0 * u.solRad,
                    cme_expansion=False,
                    cme_fixed_duration=False)

    model.solve([cme])
    cme_out = model.cmes[0]

    # Compute the transit time
    stats = cme_out.compute_arrival_at_body('EARTH')
    tt = stats['t_transit'].value
    v_1au = stats['v'].value
    
    return tt, v_1au
    
def fixed_duration_rc(onecme,model,duration,rmin=rmin):
    '''Solve HUXt using a fixed pulse duration cone CME for Richardson_Cane cone cmes'''
    
    # CME initialisation date
    t_cme = onecme['Time_21.5']
    # CME initialisation relative to model initialisation, in days
    dt_cme = (t_cme - model.time_init).jd * u.day
    
    cme = H.ConeCME(t_launch=dt_cme,
                    longitude=onecme['lon'] * u.deg,
                    latitude=onecme['lat'] * u.deg,
                    initial_height=rmin,
                    width=2.0 * onecme['ang'] * u.deg,
                    v=onecme['V'] * (u.km / u.s),
                    thickness=0.0 * u.solRad,
                    cme_expansion=False,
                    cme_fixed_duration=True,
                    fixed_duration=duration * 60 * 60 * u.s)

    model.solve([cme])
    cme_out = model.cmes[0]

    # Compute the transit time
    stats = cme_out.compute_arrival_at_body('EARTH')
    tt = stats['t_transit'].value
    v_1au = stats['v'].value
    
    return tt, v_1au
    
def fixed_duration_met(onecme,model,duration,rmin=rmin):
    '''Solve HUXt using a fixed pulse duration cone CME for Met Office cone cmes'''
    
    # CME initialisation date
    t_cme = onecme['Time_met']
    # CME initialisation relative to model initialisation, in days
    dt_cme = (t_cme - model.time_init).jd * u.day
    
    cme = H.ConeCME(t_launch=dt_cme,
                    longitude=onecme['lon_met'] * u.deg,
                    latitude=onecme['lat_met'] * u.deg,
                    initial_height=rmin,
                    width=2.0 * onecme['ang_met'] * u.deg,
                    v=onecme['V_met'] * (u.km / u.s),
                    thickness=0.0 * u.solRad,
                    cme_expansion=False,
                    cme_fixed_duration=True,
                    fixed_duration=duration * 60 * 60 * u.s)

    model.solve([cme])
    cme_out = model.cmes[0]

    # Compute the transit time
    stats = cme_out.compute_arrival_at_body('EARTH')
    tt = stats['t_transit'].value
    v_1au = stats['v'].value
    
    return tt, v_1au
#===============================================================================
# <codecell> Run model

# Initialize Data Storage
transit_time_rc = []
arrival_speed_rc = []
transit_time_met = []
arrival_speed_met = []

# First row for observation
transit_time_rc.append(['Observed'] + list(crlist['tt_21']))
transit_time_met.append(['Observed'] + list(crlist['tt_21']))
arrival_speed_rc.append(['Observed'] + list(crlist['V_ICME']))
arrival_speed_met.append(['Observed'] + list(crlist['V_ICME']))

durations = np.arange(1.0, 25.0, 0.5)  # CME durations in hours
rmin = 21.5*u.solRad
rmax = 230*u.solRad #outer boundary for HUXt runs
dt_scale = 4
simtime = 27.27 * u.day

# Pre-initialize rows
sph_tt_rc = ['Spheroidal']
sph_as_rc = ['Spheroidal']
sph_tt_met = ['Spheroidal']
sph_as_met = ['Spheroidal']
tt_rc = [[f"tt_{d:.1f}h"] for d in durations]
v_rc = [[f"v_{d:.1f}h"] for d in durations]
tt_met = [[f"tt_{d:.1f}h"] for d in durations]
v_met = [[f"v_{d:.1f}h"] for d in durations]
    
# Iterate over all CMEs in crlist
for _, onecme in crlist.iterrows():
    
    #========================================================
    # Set up the model at 21.5 rS with backmapped OMNI
    date = onecme['Time_21.5']
    date_str = date.strftime("%Y%m%d") + "00"
    
    wsafilepath = os.path.join(os.path.join("/data/WSA", "wsa_gong_{}.fits".format(date_str)))
    v_wsa = Hin.get_WSA_long_profile(wsafilepath, lat=0.0 * u.deg)
    
    cr, cr_lon_init = Hin.datetime2huxtinputs(date)
    
    model = H.HUXt(v_boundary=v_wsa, cr_num=cr, cr_lon_init = cr_lon_init,simtime = simtime, r_min=rmin, r_max=rmax, dt_scale=dt_scale, latitude=0*u.deg, frame = 'synodic', track_cmes = True, lon_out = 0*u.rad)
    
    #========================================================
    # Spheroidal cone cme
    sph_tt_val = np.nan
    sph_as_val = np.nan
    sph_tt_val, sph_as_val = spheroidal_rc(onecme, model)
    sph_tt_rc.append(sph_tt_val)
    sph_as_rc.append(sph_as_val)
    
    #========================================================
    # Fixed duration cone cmes
    for i, duration in enumerate(durations):
        tt_val = np.nan
        as_val = np.nan
        tt_val, as_val = fixed_duration_rc(onecme, model,duration)
        tt_rc[i].append(tt_val)
        v_rc[i].append(as_val)
        
# Iterate over all CMEs in crlist
for _, onecme in crlist.iterrows():
    
    #========================================================
    # Set up the model at 21.5 rS with backmapped OMNI
    date = onecme['Time_met']
    date_str = date.strftime("%Y%m%d") + "00"
    
    wsafilepath = os.path.join(os.path.join("/data/WSA", "wsa_gong_{}.fits".format(date_str)))
    v_wsa = Hin.get_WSA_long_profile(wsafilepath, lat=0.0 * u.deg)
    
    cr, cr_lon_init = Hin.datetime2huxtinputs(date)
    
    model = H.HUXt(v_boundary=v_wsa, cr_num=cr, cr_lon_init = cr_lon_init,simtime = simtime, r_min=rmin, r_max=rmax, dt_scale=dt_scale, latitude=0*u.deg, frame = 'synodic', track_cmes = True, lon_out = 0*u.rad)
    
    #========================================================
    # Spheroidal cone cme
    
    sph_tt_val = np.nan
    sph_as_val = np.nan
    sph_tt_val, sph_as_val = spheroidal_met(onecme, model)
    sph_tt_met.append(sph_tt_val)
    sph_as_met.append(sph_as_val)
    
    #========================================================
    # Fixed duration cone cmes
    for i, duration in enumerate(durations):
        
        tt_val = np.nan
        as_val = np.nan
        tt_val, as_val = fixed_duration_met(onecme, model,duration)
        tt_met[i].append(tt_val)
        v_met[i].append(as_val)
    
transit_time_rc.append(sph_tt_rc)
arrival_speed_rc.append(sph_as_rc)
transit_time_rc.extend(tt_rc)
arrival_speed_rc.extend(v_rc)

transit_time_met.append(sph_tt_met)
arrival_speed_met.append(sph_as_met)
transit_time_met.extend(tt_met)
arrival_speed_met.extend(v_met)

#==============================================================
# Save results

pd.DataFrame(transit_time_rc).to_csv(os.path.join(project_dirs['output'],"cme_transit_time_rc.csv"), index=False, header=False)
pd.DataFrame(transit_time_met).to_csv(os.path.join(project_dirs['output'],"cme_transit_time_met.csv"), index=False, header=False)
pd.DataFrame(arrival_speed_rc).to_csv(os.path.join(project_dirs['output'],"cme_arrival_speed_rc.csv"), index=False, header=False)
pd.DataFrame(arrival_speed_met).to_csv(os.path.join(project_dirs['output'],"cme_arrival_speed_met.csv"), index=False, header=False)
