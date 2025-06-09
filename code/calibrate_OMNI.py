import numpy as np
import os
import astropy.units as u
import pandas as pd
from datetime import datetime, timedelta
from astropy.time import Time
import re

import huxt as H
import huxt_inputs as Hin
import huxt_analysis as HA

#===============================================================================
# --- ICME-CME data ------------------------------------------------------------

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
crpath = os.path.join(project_dirs['input'],'(I)CMEs_SingleEvents.csv')

# Load the CSV file into a DataFrame
crlist = pd.read_csv(crpath)

# Convert date columns to datetime objects using datetime library
date_columns = ['CME_Time', 'Time_21.5', 'ICME_Start_Time', 'ICME_End_Time', 'Disturbance_Time']
for column in date_columns:
    crlist[column] = crlist[column].apply(lambda x: convert_to_datetime(x) if isinstance(x, str) else None)

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
# --- Load ICME  ---------------------------------------------------------------

icmelist_path = os.path.join(os.getcwd(),
                                'Richardson_Cane_Processed_ICME_list.csv')

#HUXt run parameters
dt_scale = 4
rmin = 21.5*u.solRad
rmax = 230*u.solRad #outer boundary for HUXt runs

recon_noicmes = True #whether to remove ICMEs
icme_list = 'CaneRichardson'#'DONKI'
icme_buffer = 0.1*u.day
sw_buffer = 2*u.day #window over which to average to fill in the solar wind during ICMEs
#fill values for ICME removal if interpolation isn't possible
vsw_fill = 450 *u.km/u.s
bx_fill = 0

def earth_R(mjd):
    #returns the heliocentric distance of Earth (in km). Based on Franz+Harper2002
    AU = 149597870.691
    
    #first up, switch to JD.
    JD=mjd+2400000.5
    d0=JD-2451545
    T0=d0/36525


    L2=100.4664568 + 35999.3728565*T0
    g2=L2-(102.9373481+0.3225654*T0)
    g2=g2*np.pi/180 #mean anomaly

    rAU=1.00014 - 0.01671*np.cos(g2)-0.00014*np.cos(2*g2);
    R=rAU*AU
    
    return R

def ICMElist(filepath):
    # -*- coding: utf-8 -*-
    """
    A script to read and process Ian Richardson's ICME list.

    Some pre-processing is required:
        Download the following webpage as a html file: 
            http://www.srl.caltech.edu/ACE/ASC/DATA/level3/icmetable2.htm
        Open in Excel, remove the year rows, delete last column (S) which is empty
        Cut out the data table only (delete header and footer)
        Save as a CSV.

    """
    
    icmes=pd.read_csv(filepath,header=None)
    #delete the first row
    icmes.drop(icmes.index[0], inplace=True)
    icmes.index = range(len(icmes))
    
    for rownum in range(0,len(icmes)):
        for colnum in range(0,3):
            #convert the three date stamps
            datestr=icmes[colnum][rownum]
            year=int(datestr[:4])
            month=int(datestr[5:7])
            day=int(datestr[8:10])
            hour=int(datestr[11:13])
            minute=int(datestr[13:15])
            #icmes.set_value(rownum,colnum,datetime(year,month, day,hour,minute,0))
            icmes.at[rownum,colnum] = datetime(year,month, day,hour,minute,0)
            
            #print(datestr)
            
        #tidy up the plasma properties
        for paramno in range(10,17):
            dv=str(icmes[paramno][rownum])
            
            #print(str(paramno)+ ' ' + dv)
            
            if dv == '...' or dv == 'dg' or dv == 'nan' or dv == '... P' or dv == '... Q':
                #icmes.set_value(rownum,paramno,np.nan)
                icmes.at[rownum,paramno] = np.nan
            else:
                #remove any remaining non-numeric characters
                dv=re.sub('[^0-9]','', dv)
                #icmes.set_value(rownum,paramno,float(dv))
                icmes.at[rownum,paramno] = float(dv)
        
    
    #chage teh column headings
    icmes=icmes.rename(columns = {0:'Shock_time',
                                  1:'ICME_start',
                                  2:'ICME_end',
                                  10:'dV',
                                  11: 'V_mean',
                                  12:'V_max',
                                  13:'Bmag',
                                  14:'MCflag',
                                  15:'Dst',
                                  16:'V_transit'})
    return icmes
  

#===============================================================================

def preprocess_omni(cme,recon_noicmes=recon_noicmes, icme_list=icme_list, icmelist_path=icmelist_path, icme_buffer=icme_buffer, sw_buffer=sw_buffer, rmin=rmin):
    '''Download OMNI data and produce boundary conditions for a given CME'''
    
    # Download and process OMNI
    icme_time = cme['CME_Time']
    
    # Compute the run start and end times so that the ICME is in the centre of the window
    run_start = icme_time - timedelta(days=13.6)
    run_stop =  icme_time + timedelta(days=13.6)
    simtime = (run_stop-run_start).days * u.day

    # Download an additional 28 days either side
    dl_starttime = run_start - timedelta(days=28*1.5)
    dl_endtime = run_stop + timedelta(days=28*1.5)
    omni = Hin.get_omni(dl_starttime, dl_endtime)

    #If the data don't span a large enough time range, repeat the last 27 days
    data_end_date = omni['datetime'][len(omni)-1]
    if data_end_date < run_stop:
        mask = (omni['datetime'] >= data_end_date - timedelta(days = 27.27))
        datachunk = omni[mask]
        datachunk.loc[:,'datetime'] = datachunk['datetime'] + timedelta(days = 27.27)
        datachunk.loc[:,'mjd'] = datachunk['mjd'] + 27.27
        #concatonate the dataframes
        omni = pd.concat([omni, datachunk], ignore_index=True)
        
    # Remove ICME from OMNI
    obs = omni.copy()
    if recon_noicmes:
    
        #create a copy of the OMNI data for ICME removal, before any end padding
        omni_noicmes = omni.copy()
        #load the DONKI ICME list
        if icme_list == 'DONKI':
            icmes = Hin.get_DONKI_ICMEs(dl_starttime, dl_endtime)
        elif icme_list == 'CaneRichardson':
            icmes = ICMElist(icmelist_path)
    
        #remove all ICMEs
        omni_noicmes =  Hin.remove_ICMEs(omni, icmes, interpolate = True,
                    icme_buffer = icme_buffer, interp_buffer = sw_buffer,
                    params = ['V', 'BX_GSE'], fill_vals = None)
      
        #if the no icme data don't span a large enough time range, repeat the last 27 days
        data_end_date = omni_noicmes['datetime'][len(omni_noicmes)-1]
        if data_end_date < run_stop:
            mask = (omni_noicmes['datetime'] >= data_end_date - timedelta(days = 27.27))
            datachunk = omni_noicmes[mask]
            datachunk.loc[:,'datetime'] = datachunk['datetime'] + timedelta(days = 27.27)
            datachunk.loc[:,'mjd'] = datachunk['mjd'] + 27.27
            #concatonate the dataframes
            omni_noicmes = pd.concat([omni_noicmes, datachunk], ignore_index=True)
            
        obs = omni_noicmes.copy()
        
    # Set up inner boundary conditions for HUXt from corotation of OMNI
    # Create  vCarr array with the omni time series at 1 AU
    time1au, v1au, b1au = Hin.generate_vCarr_from_OMNI(run_start, run_stop,
                                                            omni_input =obs,
                                                            corot_type='both')

    #now map each timestep back to the inner boundary
    vcarr_rmin = v1au.copy()
    
    for i in range(0, len(time1au)):
        #get the Earth heliocentric distance at this time
        Earth_R_km = earth_R(time1au[i]) *u.km
        #Map from 215 rto 21.5 rS
        vcarr_rmin[:,i] = Hin.map_v_boundary_inwards(v1au[:,i]*u.km/u.s,
                                    Earth_R_km.to(u.solRad), rmin)
                                    
    return simtime, vcarr_rmin, time1au, run_start
    
def spheroidal(onecme,model,rmin=rmin):
    '''Solve HUXt using a spheroidal cone CME'''
    
    cme = H.ConeCME(t_launch=13.6 * u.day,
                    longitude=onecme['lon'] * u.deg,
                    latitude=onecme['lat'] * u.deg,
                    initial_height=rmin,
                    width=2.0 * onecme['Ang_rad'] * u.deg,
                    v=onecme['V'] * (u.km / u.s),
                    thickness=0.0 * u.solRad,
                    cme_fixed_duration=False)

    model.solve([cme])
    cme_out = model.cmes[0]

    # Compute the transit time
    stats = cme_out.compute_arrival_at_body('EARTH')
    tt = stats['t_transit'].value
    arrival_time = stats['t_arrive']
    
    # Find the arrival speed within 1 day of arrival
    earth_series = HA.get_observer_timeseries(model, observer='Earth')
    mask = (Time(earth_series['time']) >= arrival_time) & (
            Time(earth_series['time']) <= arrival_time + 3 * u.day)
    v_1au = earth_series.loc[mask, 'vsw'].max()
    
    return tt, v_1au
    
def fixed_duration(onecme,model,duration,rmin=rmin):
    '''Solve HUXt using a fixed pulse duration cone CME'''
    
    cme = H.ConeCME(t_launch=13.6 * u.day,
                    longitude=onecme['lon'] * u.deg,
                    latitude=onecme['lat'] * u.deg,
                    initial_height=rmin,
                    width=2.0 * onecme['Ang_rad'] * u.deg,
                    v=onecme['V'] * (u.km / u.s),
                    thickness=0.0 * u.solRad,
                    cme_fixed_duration=True,
                    fixed_duration=duration * 60 * 60 * u.s)

    model.solve([cme])
    cme_out = model.cmes[0]

    # Compute the transit time
    stats = cme_out.compute_arrival_at_body('EARTH')
    tt = stats['t_transit'].value
    arrival_time = stats['t_arrive']
    
    # Find the arrival speed within 1 day of arrival
    earth_series = HA.get_observer_timeseries(model, observer='Earth')
    mask = (Time(earth_series['time']) >= arrival_time) & (
            Time(earth_series['time']) <= arrival_time + 3 * u.day)
    v_1au = earth_series.loc[mask, 'vsw'].max()
    
    return tt, v_1au
    
#=====Run Model=========================================================
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

durations = np.arange(0.1, 20.0, 0.1)  # CME durations in hours

# Pre-initialize rows
sph_tt_row = ['Spheroidal']
sph_as_row = ['Spheroidal']
tt_rows = [[f"tt_{d:.1f}h"] for d in durations]
v_rows = [[f"v_{d:.1f}h"] for d in durations]
    
# Iterate over all CMEs in crlist
for _, onecme in crlist.iterrows():
    
    #========================================================
    # Obtain OMNI boundary conditions at 21.5 rS
    simtime, vcarr_rmin, time1au, run_start = preprocess_omni(onecme)
    
    #========================================================
    # Set up the model at 21.5 rS with backmapped OMNI
    model = Hin.set_time_dependent_boundary(vcarr_rmin, time1au,
                        run_start, simtime = simtime, r_min=rmin, r_max=rmax,
                        dt_scale=dt_scale, latitude=0*u.deg, frame = 'sidereal',
                        track_cmes = True)
    
    #========================================================
    # Spheroidal cone cme
    tt_val, as_val = spheroidal(onecme, model)
    sph_tt_row.append(tt_val)
    sph_as_row.append(as_val)
    
    #========================================================
    # Fixed duration cone cmes
    for i, duration in enumerate(durations):
        tt_val, as_val = fixed_duration(onecme, model, duration)
        tt_rows[i].append(tt_val)
        v_rows[i].append(as_val)
    
transit_time.append(sph_tt_row)
arrival_speed.append(sph_as_row)
transit_time.extend(tt_rows)
arrival_speed.extend(v_rows)

#==============================================================
# Save results
data_dir = project_dirs['output']
out_path = os.path.join(data_dir)

pd.DataFrame(transit_time).to_csv(os.path.join(out_path,"cme_transit_time.csv"), index=False, header=False)
pd.DataFrame(arrival_speed).to_csv(os.path.join(out_path,"cme_arrival_speed.csv"), index=False, header=False)
