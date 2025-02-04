# -*- coding: utf-8 -*-
"""
A script to investigate the effect of CME perturbation shape/duration on 
model transit time and arrival speed

NOTE: Uses its own HUXt install.

@author: vy902033
"""

import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pandas as pd
from datetime import datetime, timedelta
import os
from astropy.time import Time, TimeDelta
import copy

import huxt as H
import huxt_analysis as HA

import huxt_sine_cone as Hsin

# crpath = os.path.join(os.getenv('DBOX'), 'Apps','Overleaf', 'Cone model CMEs',
#                       'RicCan_Dtags_blair_2.csv')
crpath = os.path.join('RicCan_Dtags_blair_2.csv')
r_in = 21.5*u.solRad
vsw =350
# <codecell> 2d plots of CME properties




#==============================================================================
#CME width variation
#==============================================================================
#  Form boundary conditions - background wind of 400 km/s with two fast streams.
v_boundary = np.ones(128) * 300 * (u.km/u.s)

#  Add spheroidal CMEs
cme_sphere1 = H.ConeCME(t_launch=0*u.day, longitude=360*u.deg, latitude = 0*u.deg, 
                width=120*u.deg, v=600*(u.km/u.s), thickness=0*u.solRad,
                initial_height = r_in,
                cme_expansion = False,
                cme_fixed_duration = False, 
                fixed_duration = 10*60*60*u.s)

cme_sphere2 = H.ConeCME(t_launch=0*u.day, longitude=90*u.deg, latitude = 0*u.deg, 
                width=30*u.deg, v=600*(u.km/u.s), thickness=0*u.solRad,
                initial_height = r_in,
                cme_expansion = False,
                cme_fixed_duration = False, 
                fixed_duration = 10*60*60*u.s)

cme_sphere3 = H.ConeCME(t_launch=0*u.day, longitude=180*u.deg, latitude = 0*u.deg, 
                width=90*u.deg, v=600*(u.km/u.s), thickness=0*u.solRad,
                initial_height = r_in,
                cme_expansion = False,
                cme_fixed_duration = False, 
                fixed_duration = 10*60*60*u.s)

cme_sphere4 = H.ConeCME(t_launch=0*u.day, longitude=270*u.deg, latitude = 0*u.deg, 
                width=60*u.deg, v=600*(u.km/u.s), thickness=0*u.solRad,
                initial_height = r_in,
                cme_expansion = False,
                cme_fixed_duration = False, 
                fixed_duration = 10*60*60*u.s)


cme_sphere_list = [cme_sphere1, cme_sphere2, cme_sphere3, cme_sphere4]
model_sphere = H.HUXt(v_boundary=v_boundary, latitude = 0*u.deg, 
                      r_min=r_in, simtime=6*u.day, dt_scale=4)
model_sphere.solve(cme_sphere_list)

cme_fixed_list = copy.deepcopy(cme_sphere_list)
for cme in cme_fixed_list:
    cme.cme_fixed_duration = True

#  Setup HUXt to do a 5 day simulation, with model output every 4 timesteps (roughly half and hour time step)
model_fixed = H.HUXt(v_boundary=v_boundary, latitude = 0*u.deg, 
                     r_min=r_in, simtime=6*u.day, dt_scale=4)
model_fixed.solve(cme_fixed_list, tag='cone_cme_test')

#==============================================================================
#2d plots of CME width variations
#==============================================================================

# Plot this out
t_interest = 3.35*u.day

fig = plt.figure(figsize=(8,12))
gs = fig.add_gridspec(3, 2, width_ratios=[1, 1], height_ratios = [1, 1, 0.05])


ax1 = fig.add_subplot(gs[0,0], projection='polar')
ax2 = fig.add_subplot(gs[0,1], projection='polar')

HA.plot(model_sphere, t_interest, minimalplot=True, fighandle = fig, axhandle = ax1)
HA.plot(model_fixed, t_interest, minimalplot=True, fighandle = fig, axhandle = ax2)

ax1.set_title(r'(a) Spheroidal CMEs', fontsize = 14)
ax2.set_title(r'(b) Fixed-$\Delta t$ CMEs', fontsize = 14)

for ax in [ax1, ax2]:
    r = 70

    ax.text(0, 1.1*r, r'$\alpha = $' + str(int(cme_sphere1.width.value/2)) +r'$^\circ$', fontsize = 14, 
             horizontalalignment = 'left', verticalalignment = 'center',
             bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.5'))
    ax.text(np.pi/2, r*1.5, r'$\alpha = $' + str(int(cme_sphere2.width.value/2)) +r'$^\circ$', fontsize = 14, 
             horizontalalignment = 'center', verticalalignment = 'center',
             bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.5'))
    ax.text(np.pi, r, r'$\alpha = $' + str(int(cme_sphere3.width.value/2)) +r'$^\circ$', fontsize = 14, 
             horizontalalignment = 'right', verticalalignment = 'center',
             bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.5'))
    ax.text(3*np.pi/2, r*1.5, r'$\alpha = $' + str(int(cme_sphere4.width.value/2)) +r'$^\circ$', fontsize = 14, 
             horizontalalignment = 'center', verticalalignment = 'center',
             bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.5'))
    
    #plot 1 AU
    radius = 215
    theta = np.linspace(0, 2 * np.pi, 100)
    
    # Plot the circle (constant radius, varying angle)
    ax.plot(theta, np.full_like(theta, radius), 'w--' )
    
    r_ticks = np.arange(0, 215, 215/4)  # Adjust this range as needed
    ax.set_yticks(r_ticks)

infotext = r'$V_{CME} = $' + str(int(cme_sphere1.v.value)) + ' km/s'  + '\n' + 't = ' + str(t_interest)
fig.text(0.5, 0.91, infotext, ha='center', va='center', fontsize=14, color='black')    






#==============================================================================
# CME speed variation
#==============================================================================
#  Form boundary conditions - background wind of 400 km/s with two fast streams.
v_boundary = np.ones(128) * 300 * (u.km/u.s)

width = 60*u.deg

#  Add spheroidal CMEs
cme_sphere1 = H.ConeCME(t_launch=0*u.day, longitude=360*u.deg, latitude = 0*u.deg, 
                width=width, v=300*(u.km/u.s), thickness=0*u.solRad,
                initial_height = r_in,
                cme_expansion = False,
                cme_fixed_duration = False, 
                fixed_duration = 10*60*60*u.s)

cme_sphere2 = H.ConeCME(t_launch=0*u.day, longitude=90*u.deg, latitude = 0*u.deg, 
                width=width, v=600*(u.km/u.s), thickness=0*u.solRad,
                initial_height = r_in,
                cme_expansion = False,
                cme_fixed_duration = False, 
                fixed_duration = 10*60*60*u.s)

cme_sphere3 = H.ConeCME(t_launch=0*u.day, longitude=180*u.deg, latitude = 0*u.deg, 
                width=width, v=1200*(u.km/u.s), thickness=0*u.solRad,
                initial_height = r_in,
                cme_expansion = False,
                cme_fixed_duration = False, 
                fixed_duration = 10*60*60*u.s)

cme_sphere4 = H.ConeCME(t_launch=0*u.day, longitude=270*u.deg, latitude = 0*u.deg, 
                width=width, v=2400*(u.km/u.s), thickness=0*u.solRad,
                initial_height = r_in,
                cme_expansion = False,
                cme_fixed_duration = False, 
                fixed_duration = 10*60*60*u.s)

cme_sphere_list = [cme_sphere1, cme_sphere2, cme_sphere3, cme_sphere4]
model_sphere = H.HUXt(v_boundary=v_boundary, latitude = 0*u.deg,
                      r_min=r_in, simtime=6*u.day, dt_scale=4)
model_sphere.solve(cme_sphere_list)


cme_fixed_list = copy.deepcopy(cme_sphere_list)
for cme in cme_fixed_list:
    cme.cme_fixed_duration = True

#  Setup HUXt to do a 5 day simulation, with model output every 4 timesteps (roughly half and hour time step)
model_fixed = H.HUXt(v_boundary=v_boundary, latitude = 0*u.deg, 
                     r_min=r_in, simtime=6*u.day, dt_scale=4)
model_fixed.solve(cme_fixed_list, tag='cone_cme_test')

#==============================================================================
#2d plots of CME speed variation
#==============================================================================


# Plot this out
t_interest = 3.05*u.day

# fig = plt.figure(figsize=(13,6))
# gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 0.05])

ax1 = fig.add_subplot(gs[1,0], projection='polar')
ax2 = fig.add_subplot(gs[1,1], projection='polar')

HA.plot(model_sphere, t_interest, minimalplot=True, fighandle = fig, axhandle = ax1)
HA.plot(model_fixed, t_interest, minimalplot=True, fighandle = fig, axhandle = ax2)

ax1.set_title(r'(c) Spheroidal CMEs', fontsize = 14)
ax2.set_title(r'(d) Fixed-$\Delta t$ CMEs', fontsize = 14)

for ax in [ax1, ax2]:
    r = 30

    ax.text(0, 170,  str(int(cme_sphere1.v.value)) + '\n' + ' km/s' 
            #+ '\n' + r'$t_0 = $' + str(cme_sphere1.t_launch.value)
            , fontsize = 14, horizontalalignment = 'left', verticalalignment = 'center',
             bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.5'))
    ax.text(np.pi/2, r*3,  str(int(cme_sphere2.v.value)) + ' km/s' 
            #+ '\n' + r'$t_0 = $' + str(cme_sphere2.t_launch.value)
            , fontsize = 14, horizontalalignment = 'center', verticalalignment = 'center',
             bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.5'))
    ax.text(np.pi, 80,  str(int(cme_sphere3.v.value)) + '\n' + ' km/s'
            #+ '\n' + r'$t_0 = $' + str(cme_sphere3.t_launch.value)
            , fontsize = 14, horizontalalignment = 'right', verticalalignment = 'center',
             bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.5'))
    ax.text(3*np.pi/2, r*3, str(int(cme_sphere4.v.value)) + ' km/s' 
            #+ '\n' + r'$t_0 = $' + str(cme_sphere4.t_launch.value)
            , fontsize = 14, horizontalalignment = 'center', verticalalignment = 'center',
             bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.5'))

    #plot 1 AU
    radius = 215
    theta = np.linspace(0, 2 * np.pi, 100)
    
    # Plot the circle (constant radius, varying angle)
    ax.plot(theta, np.full_like(theta, radius), 'w--' )
    
    r_ticks = np.arange(0, 215, 215/4)  # Adjust this range as needed
    ax.set_yticks(r_ticks)


infotext = r'$\alpha = $' + str(int(cme_sphere1.width.value)) + r'$^\circ$'   + '\n' + 't = ' + str(t_interest)
fig.text(0.5, 0.47, infotext, ha='center', va='center', fontsize=14, color='black')
 
    
# Create an axis for the colorbar in the GridSpec
cax = fig.add_subplot(gs[2, 0:])

# Add a single colorbar for both subplots
colorbar = fig.colorbar(ax2.collections[0], cax=cax
                        , orientation='horizontal')

# Set colorbar label
colorbar.set_label('Speed [km/s]', fontsize=14)
colorbar.ax.tick_params(labelsize=12)

plt.tight_layout()
plt.show()


# <codecell> Gopalswamy 2001 data

data_str = """
370,332,106.5
460,211,85.8
460,804,74.5
470,830,87.5
400,247,127.4
420,306,75.5
370,417,71.1
470,124,105.2
400,427,106.5
470,487,102.5
475,355,88.9
430,523,103.5
500,493,95.6
450,830,71.3
510,206,80.6
380,665,104.5
370,347,111.5
410,446,99.5
400,275,97
380,155,98.2
600,1016,68
650,1044,45.9
520,307,81.9
420,239,90.4
530,921,78.7
620,1123,60.2
435,835,85.1
460,282,88
550,362,98.4
580,1147,69.9
490,676,83.6
405,1079,54.6
500,953,98.5
410,222,81.4
480,1080,72.5
540,1009,60.5
400,550,92.5
490,447,97.5
560,551,97.5
490,641,78.9
470,396,112.4
760,1098,44.1
580,471,94.8
540,453,109.6
1100,1674,35.1
460,532,83.5
470,832,77
"""

# Convert the string to a numpy array
gopal = np.array([list(map(float, line.split(','))) for line in data_str.strip().split('\n')])


#correct for projection
gopal[:,1] = gopal[:,1]*1


# <codecell> Read in Blair's pairing on DONKI and CR2003


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
crlist['Time_Error'] = crlist['Time_Error'].apply(lambda x: convert_to_timedelta(x) if isinstance(x, str) else None)

# Convert Time_Error column from timedelta to days
def convert_timedelta_to_days(td):
    return td.days + td.seconds / 86400 if isinstance(td, timedelta) else None

crlist['Time_Error_Days'] = crlist['Time_Error'].apply(convert_timedelta_to_days)


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
crlist = crlist[abs(crlist['Time_Error_Days']) <= 0.08]
crlist = crlist.reset_index()

print('N = ' + str(len(crlist)))


#also check the CMEs are Earth directed
crlist['angdist'] = np.sqrt(crlist['lat']*crlist['lat'] + crlist['lon']*crlist['lon'])
crlist['Y'] = crlist['angdist']/crlist['Ang_rad']
#crlist = crlist[abs(crlist['Y']) <= 0.5]

print('N = ' + str(len(crlist)))

def lin_correl(x_vals, y_vals):
    #computes correlation coefficient ignoring NANs. Returns r and N
    
    # Remove NaN values from the data
    valid_indices = ~np.isnan(x_vals) & ~np.isnan(y_vals)
    x_values = x_vals[valid_indices]
    y_values = y_vals[valid_indices]
    
    #compute correlation
    rl = np.corrcoef(x_values, y_values)
    
    return rl[0,1], len(x_values)


# find the linear fit between the CME and ICME speeds
# Calculate the slope and intercept of the best-fit line
slope, intercept = np.polyfit(crlist['V_max'], crlist['V'], 1)
xt = 0.95
yt = 0.08

plt.figure(figsize=(10,9))
ax1 = plt.subplot(221)
ax2 = plt.subplot(222)
ax3 = plt.subplot(223)
ax4 = plt.subplot(224)

xvals = crlist['V']
yvals = crlist['Ang_rad']
r = lin_correl(xvals, yvals)
ax1.plot(xvals, yvals, 'ko')
ax1.set_xlabel(r'CME speed, $V_{CME}$ [km s$^{-1}$]')
ax1.set_ylabel(r'CME angular half-width, $\alpha$ [deg]')
ax1.text(xt, yt, r'(a) $r_L =$ ' + f"{r[0]:.2f}", transform=ax1.transAxes,
        fontsize=14, verticalalignment='top', horizontalalignment='right',
        bbox=dict(facecolor='white', edgecolor='none', pad=3))
ax1.set_ylim([0,90])


# xvals = crlist['duration_1au']
# yvals = crlist['Ang_rad']
# r = lin_correl(xvals, yvals)
# ax2.plot(xvals, yvals, 'ko')
# ax2.set_xlabel(r'ICME duration [days]')
# #ax2.set_ylabel(r'CME half-angular width, $\alpha$ [deg]')
# ax2.set_title(r'(b) $r_L =$ ' + f"{r[0]:.2f}", fontsize = 14)
# ax2.set_xlim([0,3])

xvals = crlist['duration_1au']* crlist['V_ICME'] *(24*60*60)/(1.5e8)
yvals = crlist['Ang_rad']
r = lin_correl(xvals, yvals)
# mask = crlist['V_max'] < 500
# ax2.plot(xvals[mask], yvals[mask] , 'bo')
# mask = (crlist['V_max'] >= 500) & (crlist['V_max'] < 600) 
# ax2.plot(xvals[mask], yvals[mask] , 'ko')
# mask = crlist['V_max'] >= 600
# ax2.plot(xvals[mask], yvals[mask] , 'ro')
ax2.plot(xvals, yvals , 'ko')
ax2.set_xlabel(r'ICME radial width, $r_{1AU}$ [AU]')
ax2.set_xlim([0,0.7])
#ax2.set_ylabel(r'CME half-angular width, $\alpha$ [deg]')
ax2.set_yticklabels([])
ax2.text(xt, yt, r'(b) $r_L =$ ' + f"{r[0]:.2f}", transform=ax2.transAxes,
        fontsize=14, verticalalignment='top', horizontalalignment='right',
        bbox=dict(facecolor='white', edgecolor='none', pad=3))
#plot the model relations
angs = np.arange(0,91,5)*np.pi/180
ax2.plot(2*0.1*np.sin(angs), angs*180/np.pi,'r', label = 'Spheroidal CME (sin)')
ax2.plot(2*0.1*np.tan(angs), angs*180/np.pi,'r--', label = 'Spheroidal CME (tan)')
ax2.plot(2*0.1*np.sin(angs) + 0.2, angs*180/np.pi,'r:', label = 'Spherocylinder')
#ax2.plot(1.5*400*8*60*60/1.5e8*np.array([1,1]),[15,75],'b:', label = 'Fixed $\Delta t, V_{CME} = 400$ '  )
#ax2.plot(1.5*2000*8*60*60/1.5e8*np.array([1,1]),[15,75],'b--', label = 'Fixed $\Delta t, V_{CME} = 2000$ '  )
ax2.legend(loc = 'center',bbox_to_anchor=(0.7, 1.02), ncol=1)
ax2.set_ylim([0,90])

xvals = crlist['V']
yvals = crlist['duration_1au']* crlist['V_ICME'] *(24*60*60)/(1.5e8)
r = lin_correl(xvals, yvals)
ax3.plot(xvals, yvals, 'ko')
ax3.set_xlabel(r'CME speed, $V_{CME}$  [km s$^{-1}$]')
ax3.set_ylabel(r'ICME radial width, $r_{1AU}$ [AU]')
ax3.set_ylim([0, 0.8])
Vs = np.arange(150,2200,50)
ax3.plot(Vs, 2*Vs*8*60*60/1.5e8,'b', label = r'Fixed-$\Delta t$ CME')
ax3.legend()
ax3.text(xt, yt, r'(c) $r_L =$ ' + f"{r[0]:.2f}", transform=ax3.transAxes,
        fontsize=14, verticalalignment='top', horizontalalignment='right',
        bbox=dict(facecolor='white', edgecolor='none', pad=3))


xvals = crlist['V_max']
yvals = crlist['duration_1au']* crlist['V_ICME'] *(24*60*60)/(1.5e8)
r = lin_correl(xvals, yvals)
ax4.plot(xvals, yvals, 'ko')
ax4.set_xlabel(r'ICME speed, $V_{1AU}$  [km s$^{-1}$]')
ax4.set_ylim([0, 0.8])
#ax4.set_ylabel(r'ICME radial width, $r_{1AU}$ [AU]')
ax4.set_yticklabels([])
Vs = np.arange(300,901,50)
ax4.plot(Vs, 2*slope*(Vs+intercept)*8*60*60/1.5e8,'b', label = r'Fixed-$\Delta t$ CME')
ax4.legend()
ax4.text(xt, yt, r'(d) $r_L =$ ' + f"{r[0]:.2f}", transform=ax4.transAxes,
        fontsize=14, verticalalignment='top', horizontalalignment='right',
        bbox=dict(facecolor='white', edgecolor='none', pad=3))

plt.tight_layout()


# plt.figure(figsize=(13,5))
# ax1 = plt.subplot(131)
# ax2 = plt.subplot(132)
# ax3 = plt.subplot(133)

# xvals = crlist['V_ICME']
# yvals = crlist['duration_1au']
# r = lin_correl(xvals, yvals)
# ax1.plot(xvals, yvals, 'ko')
# ax1.set_xlabel(r'ICME speed [km s$^{-1}$]')
# ax1.set_ylabel(r'ICME duration [days]')
# ax1.set_title(r'(a) $r_L =$ ' + f"{r[0]:.2f}", fontsize = 14)

# xvals = crlist['V_max']
# yvals = crlist['duration_1au']
# r = lin_correl(xvals, yvals)
# ax2.plot(xvals, yvals, 'ko')
# ax2.set_xlabel(r'ICME max speed [km s$^{-1}$]')
# ax2.set_ylabel(r'ICME duration [days]')
# ax2.set_title(r'(a) $r_L =$ ' + f"{r[0]:.2f}", fontsize = 14)

# <codecell> plot observations for transit properties

plt.figure(figsize=(14,5))

plotcr = True

ax1 = plt.subplot(131)
ax2 = plt.subplot(132)
ax3 = plt.subplot(133)

colours = ['b', 'k', 'r']
linewidth=2



if plotcr:
    # ax1.plot(crlist['V'], crlist['tt_21'], 'ko')
    # ax4.plot(crlist['V'], crlist['tt_21'], 'ko')
    
    # ax2.plot(crlist['V'], crlist['V_max'], 'ko')
    # ax5.plot(crlist['V'], crlist['V_max'], 'ko')
    
    # ax3.plot(crlist['V'], 1000*(crlist['V_max'] - crlist['V'])/( crlist['tt_21']*24*60*60), 'ko')
    # ax6.plot(crlist['V'], 1000*(crlist['V_max'] - crlist['V'])/( crlist['tt_21']*24*60*60), 'ko')
    alpha = 0.5
    
    
    for ax in [ax1]:
        mask = crlist['Ang_rad'] < 35
        ax.plot(crlist.loc[mask,'tt_21'], crlist.loc[mask,'V'], 'o', color=colours[0], 
                alpha = alpha, label=r'$\alpha < 35^\circ$')
        mask = (crlist['Ang_rad'] >= 35) & (crlist['Ang_rad'] <55)
        ax.plot(crlist.loc[mask,'tt_21'],  crlist.loc[mask,'V'], 'o', color=colours[1], 
                alpha = alpha, label=r'$35^\circ \geq \alpha < 55^\circ$')
        mask = crlist['Ang_rad'] >=55
        ax.plot(crlist.loc[mask,'tt_21'], crlist.loc[mask,'V'], 'o', color=colours[2], 
                alpha = alpha, label=r'$\alpha \geq 55^\circ$')
        
    for ax in [ax2]:
        mask = crlist['Ang_rad'] < 35
        ax.plot(crlist.loc[mask,'V_max'], crlist.loc[mask,'V'],  'o', color=colours[0], 
                alpha = alpha, label=r'$\alpha < 35^\circ$')
        mask = (crlist['Ang_rad'] >= 35) & (crlist['Ang_rad'] <55)
        ax.plot(crlist.loc[mask,'V_max'], crlist.loc[mask,'V'], 'o', color=colours[1], 
                alpha = alpha, label=r'$35^\circ \geq \alpha < 55^\circ$')
        mask = crlist['Ang_rad'] >55
        ax.plot(crlist.loc[mask,'V_max'], crlist.loc[mask,'V'], 'o', color=colours[2], 
                alpha = alpha, label=r'$\alpha \geq 55^\circ$')
        
    for ax in [ax3]:
        mask = crlist['Ang_rad'] < 35
        ax.plot(1000*(crlist.loc[mask,'V_max'] - crlist.loc[mask,'V'])/( crlist.loc[mask,'tt_21']*24*60*60), 
                 crlist.loc[mask,'V'], 'o', color=colours[0], alpha = alpha)
        mask = (crlist['Ang_rad'] >= 35) & (crlist['Ang_rad'] <55)
        ax.plot(1000*(crlist.loc[mask,'V_max'] 
                                             - crlist.loc[mask,'V'])/( crlist.loc[mask,'tt_21']*24*60*60), 
                 crlist.loc[mask,'V'], 'o', color=colours[1], alpha = alpha)
        mask = crlist['Ang_rad'] >55
        ax.plot(1000*(crlist.loc[mask,'V_max'] 
                                             - crlist.loc[mask,'V'])/( crlist.loc[mask,'tt_21']*24*60*60), 
                 crlist.loc[mask,'V'], 'o', color=colours[2], alpha = alpha)

#add panel labels
for ax, label in zip([ax1, ax2, ax3], 
                     ['(a)', '(b)', '(c)', '(d)', '(e)', '(f)']):
    ax.text(0.9, 0.9, label, transform=ax.transAxes,
            fontsize=14, verticalalignment='top', horizontalalignment='right',
            bbox=dict(facecolor='white', edgecolor='none', pad=3))



ax1.set_xlabel(r'Transit time, $\tau$ [days]')
ax1.set_ylabel(r'CME speed, $V_{CME}$ [km s$^{-1}$]')
ax1.set_ylim([0,2300])
ax1.set_xlim([1,5])


ax2.set_xlabel(r'ICME speed, $V_{1AU}$ [km s$^{-1}$]')
#ax2.set_ylabel(r'CME speed, $V_{CME}$ [km s$^{-1}$]')
ax2.set_ylim([0,2300])
ax2.legend(loc = 'center', bbox_to_anchor=(0.5, 1.1), ncol=3)
ax2.set_yticklabels([])


ax3.set_xlabel(r'Acceleration, $a_{IP}$ [m s$^{-2}$]')
ax3.set_ylim([0,2300])
#ax3.set_ylabel(r'CME speed, $V_{CME}$ [km s$^{-1}$]')
ax3.set_yticklabels([])

#plt.tight_layout()
plt.subplots_adjust(wspace=0.1, hspace=0.1,top=0.82, bottom=0.2) 
plt.show()
# <codecell> plot transit times for various sperhical CME implementations

#look at a range of CME speeds
v_boundary = np.ones(128) * vsw * (u.km/u.s)
cme_speeds = np.arange(100,2500,50) 
cme_widths = [50, 80, 130]

plt.figure(figsize=(14,10))

plotcr = True

ax1 = plt.subplot(331)
ax2 = plt.subplot(332)
ax3 = plt.subplot(333)
ax4 = plt.subplot(334)
ax5 = plt.subplot(335)
ax6 = plt.subplot(336)
ax7 = plt.subplot(337)
ax8 = plt.subplot(338)
ax9 = plt.subplot(339)

colours = ['b', 'k', 'r']
linewidth=2



if plotcr:
  
    alpha = 0.5

    for ax in [ax1, ax4, ax7]:
        ax.plot(crlist['tt_21'], crlist['V'], 'ko', alpha = alpha)
        
    for ax in [ax2, ax5, ax8]:
        ax.plot(crlist['V_max'], crlist['V'], 'ko',  alpha = alpha)
       
    for ax in [ax3, ax6, ax9]:
        ax.plot(1000*(crlist['V_max'] 
                                             - crlist['V'])/( crlist['tt_21']*24*60*60), 
                crlist['V'],  'ko', alpha = alpha)

#add panel labels
for ax, label in zip([ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9], 
                     ['(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '(g)', '(h)', '(i)' ]):
    ax.text(0.11, 0.11, label, transform=ax.transAxes,
            fontsize=14, verticalalignment='top', horizontalalignment='right',
            bbox=dict(facecolor='white', edgecolor='none', pad=3))


#ax1.set_xlabel(r'Transit time, $\tau$ [days]')
#ax1.set_ylabel(r'CME speed, $V_{CME}$ [km s$^{-1}$]')
ax1.set_xlim([0.9,6.1])
ax1.set_xticklabels([])
#ax1.set_xlabel('CME speed [km/s]')

#ax2.set_xlabel(r'ICME speed, $V_{1AU}$ [km s$^{-1}$]')
ax2.set_xlim([190,1110])
ax2.set_yticklabels([])
ax2.set_xticklabels([])

#ax3.set_xlabel(r'Acceleration, $a_{IP}$ [m s$^{-2}$]')
ax3.set_yticklabels([])
ax3.set_xticklabels([])
ax3.set_ylabel(r'Spheroidal CME, $\tan \alpha$', rotation=270, labelpad=20)
ax3.yaxis.set_label_position('right')

#ax4.set_xlabel(r'Transit time, $\tau$ [days]')
ax4.set_ylabel(r'CME speed, $V_{CME}$ [km s$^{-1}$]')
ax4.set_xlim([0.9,6.1])
ax4.set_xticklabels([])

#ax5.set_xlabel(r'ICME speed, $V_{1AU}$ [km s$^{-1}$]')
ax5.set_xlim([190,1110])
ax5.set_yticklabels([])
ax5.set_xticklabels([])

#ax6.set_xlabel(r'Acceleration, $a_{IP}$ [m s$^{-2}$]')
ax6.set_yticklabels([])
ax6.set_xticklabels([])
ax6.set_ylabel(r'Spheroidal CME, $\sin \alpha$', rotation=270, labelpad=20)
ax6.yaxis.set_label_position('right')


#ax7.set_ylabel(r'CME speed, $V_{CME}$ [km s$^{-1}$]')
ax7.set_xlabel(r'Transit time, $\tau$ [days]')
ax7.set_xlim([0.9,6.1])

#ax8.set_ylabel(r'CME speed, $V_{CME}$ [km s$^{-1}$]')
ax8.set_xlabel(r'ICME speed, $V_{1AU}$ [km s$^{-1}$]')
ax8.set_xlim([190,1110])
ax8.set_yticklabels([])

#ax9.set_ylabel(r'CME speed, $V_{CME}$ [km s$^{-1}$]')
ax9.set_xlabel(r'Acceleration, $a_{IP}$ [m s$^{-2}$]')
ax9.set_yticklabels([])
ax9.set_ylabel(r'Spherocyclinder CME', rotation=270, labelpad=20)
ax9.yaxis.set_label_position('right')




model = H.HUXt(v_boundary=v_boundary, latitude = 0*u.deg, 
               r_min=r_in, simtime=10*u.day, 
               dt_scale=4, lon_out = 0*u.deg)
i=0
#spheroidal - tan alpha CMEs
#===============================================================================
for width in cme_widths:
    transit_times = []
    arrival_speeds = []
    ip_acc = []

    for cme_speed in cme_speeds:
        


        cme = H.ConeCME(t_launch=1*u.day, longitude=360*u.deg, latitude = 0*u.deg, 
                        initial_height = r_in,
                        width=width*u.deg, v=cme_speed* (u.km/u.s), thickness = 0*u.solRad,
                        cme_fixed_duration = False)
        model.solve([cme])
        
        
        #get the transit time
        stats = model.cmes[0].compute_arrival_at_body('EARTH')
        tt = stats['t_transit'].value
        transit_times.append(tt)
        arrival_time = stats['t_arrive']
        
        #find the arrival speed as the max withing 1 day of arrival
        earth_series = HA.get_observer_timeseries(model, observer = 'Earth')
        mask = (Time(earth_series['time']) >= arrival_time) & (
            Time(earth_series['time']) <= arrival_time + 3*u.day)
        
        v_1au = earth_series.loc[mask, 'vsw'].max()
        arrival_speeds.append(v_1au)
        
        #compute hte interplanetary acceleration
        ip_acc.append(1000* (v_1au - cme_speed) / (tt*24*60*60) )
        

    
    ax1.plot(transit_times, cme_speeds,  linewidth = linewidth,
             label = r'$\alpha$ = ' +str(width/2) +r'$^\circ$', color = colours[i])

    ax2.plot(arrival_speeds, cme_speeds,  linewidth = linewidth,
             label = r'$\alpha$ = '  +str(width/2) +r'$^\circ$',color = colours[i])
    
    ax3.plot(ip_acc, cme_speeds,  linewidth = linewidth,
             label = r'$\alpha$ = ' + str(width/2) +r'$^\circ$', color = colours[i])
    
    i=i+1



title_handle = Line2D([0], [0], color='none')
# Get existing handles and labels
handles, labels = ax1.get_legend_handles_labels()
# Add the title handle and text at the start
handles.insert(0, title_handle)
labels.insert(0, r'Spheroidal CME (tan $\alpha$)')

ax2.legend(ncol =3, loc = 'center',bbox_to_anchor=(0.5, 1.1))




#spheriodal CMEs sin alpha
#===============================================================================

transit_times = []
arrival_speeds = []
ip_acc = []
model = Hsin.HUXt(v_boundary=v_boundary, latitude = 0*u.deg, 
               r_min=r_in, simtime=10*u.day, 
               dt_scale=4, lon_out = 0*u.deg)

i=0
for width in cme_widths:
    transit_times = []
    arrival_speeds = []
    ip_acc = []

    for cme_speed in cme_speeds:

        cme = Hsin.ConeCME(t_launch=1*u.day, longitude=360*u.deg, latitude = 0*u.deg, 
                        initial_height = r_in,
                        width=width*u.deg, v=cme_speed* (u.km/u.s), thickness = 0*u.solRad,
                        cme_fixed_duration = False, fixed_duration = 8*60*60*u.s)
        model.solve([cme])
        
        
        #get the transit time
        stats = model.cmes[0].compute_arrival_at_body('EARTH')
        tt = stats['t_transit'].value
        transit_times.append(tt)
        arrival_time = stats['t_arrive']
        
        #find the arrival speed as the max withing 1 day of arrival
        earth_series = HA.get_observer_timeseries(model, observer = 'Earth')
        mask = (Time(earth_series['time']) >= arrival_time) & (
            Time(earth_series['time']) <= arrival_time + 3*u.day)
        
        v_1au = earth_series.loc[mask, 'vsw'].max()
        arrival_speeds.append(v_1au)
        
        #compute hte interplanetary acceleration
        ip_acc.append(1000* (v_1au - cme_speed) / (tt*24*60*60) )
    

    ax4.plot(transit_times, cme_speeds,   linewidth = linewidth,
             label = r'$\alpha$ = ' +str(width/2) +r'$^\circ$', color = colours[i])

    ax5.plot(arrival_speeds, cme_speeds,  linewidth = linewidth,
             label = r'$\alpha$ = '  +str(width/2) +r'$^\circ$',color = colours[i])
    
    ax6.plot(ip_acc, cme_speeds,  linewidth = linewidth,
             label = r'$\alpha$ = ' + str(width/2) +r'$^\circ$', color = colours[i])
    
    i = i+1


title_handle = Line2D([0], [0], color='none')
# Get existing handles and labels
handles, labels = ax4.get_legend_handles_labels()
# Add the title handle and text at the start
handles.insert(0, title_handle)
labels.insert(0, r'Spheroidal CME (sin $\alpha$)')

#ax4.legend()



#spherocylcinder CMEs sin alpha
#===============================================================================

transit_times = []
arrival_speeds = []
ip_acc = []
model = Hsin.HUXt(v_boundary=v_boundary, latitude = 0*u.deg, 
               r_min=r_in, simtime=10*u.day, 
               dt_scale=4, lon_out = 0*u.deg)

i=0
for width in cme_widths:
    transit_times = []
    arrival_speeds = []
    ip_acc = []

    for cme_speed in cme_speeds:

        cme = Hsin.ConeCME(t_launch=1*u.day, longitude=360*u.deg, latitude = 0*u.deg, 
                        initial_height = r_in,
                        width=width*u.deg, v=cme_speed* (u.km/u.s), thickness = 72*u.solRad,
                        cme_fixed_duration = False, fixed_duration = 8*60*60*u.s)
        model.solve([cme])
        
        
        #get the transit time
        stats = model.cmes[0].compute_arrival_at_body('EARTH')
        tt = stats['t_transit'].value
        transit_times.append(tt)
        arrival_time = stats['t_arrive']
        
        #find the arrival speed as the max withing 1 day of arrival
        earth_series = HA.get_observer_timeseries(model, observer = 'Earth')
        mask = (Time(earth_series['time']) >= arrival_time) & (
            Time(earth_series['time']) <= arrival_time + 3*u.day)
        
        v_1au = earth_series.loc[mask, 'vsw'].max()
        arrival_speeds.append(v_1au)
        
        #compute hte interplanetary acceleration
        ip_acc.append(1000* (v_1au - cme_speed) / (tt*24*60*60) )
    

    ax7.plot(transit_times, cme_speeds,  linewidth = linewidth,
             label = r'$\alpha$ = ' +str(width/2) +r'$^\circ$', color = colours[i])

    ax8.plot(arrival_speeds, cme_speeds,  linewidth = linewidth,
             label = r'$\alpha$ = '  +str(width/2) +r'$^\circ$',color = colours[i])
    
    ax9.plot(ip_acc, cme_speeds, linewidth = linewidth,
             label = r'$\alpha$ = ' + str(width/2) +r'$^\circ$', color = colours[i])
    
    i = i+1


title_handle = Line2D([0], [0], color='none')
# Get existing handles and labels
handles, labels = ax4.get_legend_handles_labels()
# Add the title handle and text at the start
handles.insert(0, title_handle)
labels.insert(0, r'Spheroidal CME (sin $\alpha$)')

#ax4.legend()

plt.subplots_adjust(wspace=0.1, hspace=0.1,top=0.92, bottom=0.1) 
plt.show()


#plt.tight_layout()

# <codecell> Now plot fixed duration CMEs

plt.figure(figsize=(14,5))

plotcr = True

ax1 = plt.subplot(131)
ax2 = plt.subplot(132)
ax3 = plt.subplot(133)

colours = ['b', 'k', 'r']
linewidth=2




if plotcr:
  
    alpha = 0.5

    for ax in [ax1]:
        ax.plot(crlist['tt_21'], crlist['V'], 'ko', alpha = alpha)
        
    for ax in [ax2]:
        ax.plot(crlist['V_max'], crlist['V'], 'ko',  alpha = alpha)
       
    for ax in [ax3]:
        ax.plot(1000*(crlist['V_max'] 
                                             - crlist['V'])/( crlist['tt_21']*24*60*60), 
                crlist['V'],  'ko', alpha = alpha)

#add panel labels
for ax, label in zip([ax1, ax2, ax3], 
                     ['(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '(g)', '(h)', '(i)' ]):
    ax.text(0.11, 0.11, label, transform=ax.transAxes,
            fontsize=14, verticalalignment='top', horizontalalignment='right',
            bbox=dict(facecolor='white', edgecolor='none', pad=3))

ax1.set_xlabel(r'Transit time, $\tau$ [days]')
ax1.set_ylabel(r'CME speed, $V_{CME}$ [km s$^{-1}$]')
ax1.set_ylim([0,2300])
ax1.set_xlim([0.9,6.1])


ax2.set_xlabel(r'ICME speed, $V_{1AU}$ [km s$^{-1}$]')
#ax2.set_ylabel(r'CME speed, $V_{CME}$ [km s$^{-1}$]')
ax2.set_ylim([0,2300])
ax2.set_xlim([190,1110])
ax2.set_yticklabels([])


ax3.set_xlabel(r'Acceleration, $a_{IP}$ [m s$^{-2}$]')
ax3.set_ylim([0,2300])
#ax3.set_ylabel(r'CME speed, $V_{CME}$ [km s$^{-1}$]')
ax3.set_ylabel(r'Fixed $\Delta t$ CME', rotation=270, labelpad=20)
ax3.set_yticklabels([])
ax3.yaxis.set_label_position('right')


plt.subplots_adjust(wspace=0.1, hspace=0.1,top=0.82, bottom=0.2) 
plt.show()



#
#===============================================================================

transit_times = []
arrival_speeds = []
ip_acc = []
model = H.HUXt(v_boundary=v_boundary, latitude = 0*u.deg, 
               r_min=r_in, simtime=10*u.day, 
               dt_scale=4, lon_out = 0*u.deg)

i=0

transit_times = []
arrival_speeds = []
ip_acc = []

width = 60

for cme_speed in cme_speeds:

    cme = H.ConeCME(t_launch=1*u.day, longitude=360*u.deg, latitude = 0*u.deg, 
                    initial_height = r_in,
                    width=width*u.deg, v=cme_speed* (u.km/u.s), thickness = 0*u.solRad,
                    cme_fixed_duration = True, fixed_duration = 8*60*60*u.s)
    model.solve([cme])
    
    
    #get the transit time
    stats = model.cmes[0].compute_arrival_at_body('EARTH')
    tt = stats['t_transit'].value
    transit_times.append(tt)
    arrival_time = stats['t_arrive']
    
    #find the arrival speed as the max withing 1 day of arrival
    earth_series = HA.get_observer_timeseries(model, observer = 'Earth')
    mask = (Time(earth_series['time']) >= arrival_time) & (
        Time(earth_series['time']) <= arrival_time + 3*u.day)
    
    v_1au = earth_series.loc[mask, 'vsw'].max()
    arrival_speeds.append(v_1au)
    
    #compute hte interplanetary acceleration
    ip_acc.append(1000* (v_1au - cme_speed) / (tt*24*60*60) )


ax1.plot(transit_times, cme_speeds,  linewidth = linewidth,
         label = r'$\alpha$ = ' +str(width/2) +r'$^\circ$', color ='r')

ax2.plot(arrival_speeds, cme_speeds,  linewidth = linewidth,
         label = r'$\alpha$ = '  +str(width/2) +r'$^\circ$',color = 'r')

ax3.plot(ip_acc, cme_speeds, linewidth = linewidth,
         label = r'$\alpha$ = ' + str(width/2) +r'$^\circ$', color = 'r')

i = i+1
# # <codecell>sherocylinder with variable r_C

# plt.figure(figsize=(14,5))

# plotcr = True

# ax1 = plt.subplot(131)
# ax2 = plt.subplot(132)
# ax3 = plt.subplot(133)

# colours = ['b', 'k', 'r']
# linewidth=2




# if plotcr:
  
#     alpha = 0.5

#     for ax in [ax1]:
#         ax.plot(crlist['tt_21'], crlist['V'], 'ko', alpha = alpha)
        
#     for ax in [ax2]:
#         ax.plot(crlist['V_max'], crlist['V'], 'ko',  alpha = alpha)
       
#     for ax in [ax3]:
#         ax.plot(1000*(crlist['V_max'] 
#                                              - crlist['V'])/( crlist['tt_21']*24*60*60), 
#                 crlist['V'],  'ko', alpha = alpha)

# #add panel labels
# for ax, label in zip([ax1, ax2, ax3], 
#                      ['(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '(g)', '(h)', '(i)' ]):
#     ax.text(0.11, 0.11, label, transform=ax.transAxes,
#             fontsize=14, verticalalignment='top', horizontalalignment='right',
#             bbox=dict(facecolor='white', edgecolor='none', pad=3))

# ax1.set_xlabel(r'Transit time, $\tau$ [days]')
# ax1.set_ylabel(r'CME speed, $V_{CME}$ [km s$^{-1}$]')
# ax1.set_ylim([0,2300])
# ax1.set_xlim([0.9,6.1])


# ax2.set_xlabel(r'ICME speed, $V_{1AU}$ [km s$^{-1}$]')
# #ax2.set_ylabel(r'CME speed, $V_{CME}$ [km s$^{-1}$]')
# ax2.set_ylim([0,2300])
# ax2.set_xlim([190,1110])
# ax2.set_yticklabels([])


# ax3.set_xlabel(r'Acceleration, $a_{IP}$ [m s$^{-2}$]')
# ax3.set_ylim([0,2300])
# #ax3.set_ylabel(r'CME speed, $V_{CME}$ [km s$^{-1}$]')
# ax3.set_ylabel(r'Fixed $\Delta t$ CME', rotation=270, labelpad=20)
# ax3.set_yticklabels([])
# ax3.yaxis.set_label_position('right')


# plt.subplots_adjust(wspace=0.1, hspace=0.1,top=0.82, bottom=0.2) 
# plt.show()



# #
# #===============================================================================



# #spheriodal CMEs sin alpha with r_C as funciton of V_CME
# #===============================================================================

# transit_times = []
# arrival_speeds = []
# ip_acc = []
# model = Hsin.HUXt(v_boundary=v_boundary, latitude = 0*u.deg, 
#                r_min=r_in, simtime=10*u.day, 
#                dt_scale=4, lon_out = 0*u.deg)

# i=0
# for width in cme_widths:
#     transit_times = []
#     arrival_speeds = []
#     ip_acc = []

#     for cme_speed in cme_speeds:

#         cme = Hsin.ConeCME(t_launch=1*u.day, longitude=360*u.deg, latitude = 0*u.deg, 
#                         initial_height = r_in,
#                         width=width*u.deg, v=cme_speed* (u.km/u.s), 
#                         thickness = (70/2000)*(cme_speed - 400)*u.solRad,
#                         cme_fixed_duration = False, fixed_duration = 8*60*60*u.s)
#         model.solve([cme])
        
        
#         #get the transit time
#         stats = model.cmes[0].compute_arrival_at_body('EARTH')
#         tt = stats['t_transit'].value
#         transit_times.append(tt)
#         arrival_time = stats['t_arrive']
        
#         #find the arrival speed as the max withing 1 day of arrival
#         earth_series = HA.get_observer_timeseries(model, observer = 'Earth')
#         mask = (Time(earth_series['time']) >= arrival_time) & (
#             Time(earth_series['time']) <= arrival_time + 3*u.day)
        
#         v_1au = earth_series.loc[mask, 'vsw'].max()
#         arrival_speeds.append(v_1au)
        
#         #compute hte interplanetary acceleration
#         ip_acc.append(1000* (v_1au - cme_speed) / (tt*24*60*60) )
    

#     ax1.plot(transit_times, cme_speeds,   linewidth = linewidth,
#              label = r'$\alpha$ = ' +str(width/2) +r'$^\circ$', color = colours[i])

#     ax2.plot(arrival_speeds, cme_speeds,  linewidth = linewidth,
#              label = r'$\alpha$ = '  +str(width/2) +r'$^\circ$',color = colours[i])
    
#     ax3.plot(ip_acc, cme_speeds,  linewidth = linewidth,
#              label = r'$\alpha$ = ' + str(width/2) +r'$^\circ$', color = colours[i])
    
#     i = i+1

# # <codecell> Reprocess the CME list to estimate properties at 10 rS, assuming
# # they were ballistically mapped to 21.5 rS

# r_inner = 10*u.solRad

# # compute 21.5-215 transit time
# crlist['time_at_10rs'] = np.nan 
# for irow in range(0, len(crlist)):
#     dt = (21.5 - r_inner.to(u.solRad).value) * 696340/ crlist.loc[irow,'V']
#     crlist.loc[irow,'time_at_10rs'] = (crlist.loc[irow,'Time_21.5'] 
#                                        - timedelta(seconds = dt))
# # recompute transit times
# crlist['tt_10'] = np.nan 
# for irow in range(0, len(crlist)):
#     crlist.loc[irow,'tt_10'] = crlist.loc[irow,'Disturbance_Time'] - crlist.loc[irow,'time_at_10rs']
# # Convert  from timedelta to days
# crlist['tt_10'] = crlist['tt_10'].apply(lambda x: x.days + x.seconds / 86400 if isinstance(x, timedelta) else None)


# v_boundary = np.ones(128) * vsw * (u.km/u.s)
# cme_speeds = np.arange(100,2500,50) 

# plt.figure(figsize=(14,5))

# plotcr = True

# ax1 = plt.subplot(131)
# ax2 = plt.subplot(132)
# ax3 = plt.subplot(133)

# colours = ['b', 'k', 'r']
# linewidth=2




# if plotcr:
  
#     alpha = 0.5

#     for ax in [ax1]:
#         ax.plot(crlist['tt_10'], crlist['V'], 'ko', alpha = alpha)
        
#     for ax in [ax2]:
#         ax.plot(crlist['V_max'], crlist['V'], 'ko',  alpha = alpha)
       
#     for ax in [ax3]:
#         ax.plot(1000*(crlist['V_max'] 
#                                              - crlist['V'])/( crlist['tt_10']*24*60*60), 
#                 crlist['V'],  'ko', alpha = alpha)

# #add panel labels
# for ax, label in zip([ax1, ax2, ax3], 
#                      ['(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '(g)', '(h)', '(i)' ]):
#     ax.text(0.11, 0.11, label, transform=ax.transAxes,
#             fontsize=14, verticalalignment='top', horizontalalignment='right',
#             bbox=dict(facecolor='white', edgecolor='none', pad=3))

# ax1.set_xlabel(r'Transit time, $\tau$ [days]')
# ax1.set_ylabel(r'CME speed, $V_{CME}$ [km s$^{-1}$]')
# ax1.set_ylim([0,2300])
# ax1.set_xlim([0.9,6.1])


# ax2.set_xlabel(r'ICME speed, $V_{1AU}$ [km s$^{-1}$]')
# #ax2.set_ylabel(r'CME speed, $V_{CME}$ [km s$^{-1}$]')
# ax2.set_ylim([0,2300])
# ax2.set_xlim([190,1110])
# ax2.set_yticklabels([])


# ax3.set_xlabel(r'Acceleration, $a_{IP}$ [m s$^{-2}$]')
# ax3.set_ylim([0,2300])
# #ax3.set_ylabel(r'CME speed, $V_{CME}$ [km s$^{-1}$]')
# ax3.set_ylabel(r'Fixed $\Delta t$ CME', rotation=270, labelpad=20)
# ax3.set_yticklabels([])
# ax3.yaxis.set_label_position('right')


# plt.subplots_adjust(wspace=0.1, hspace=0.1,top=0.82, bottom=0.2) 
# plt.show()



# transit_times = []
# arrival_speeds = []
# ip_acc = []
# model = H.HUXt(v_boundary=v_boundary, latitude = 0*u.deg, 
#                r_min=r_inner, simtime=10*u.day, 
#                dt_scale=4, lon_out = 0*u.deg)

# i=0
# width = 60
# transit_times = []
# arrival_speeds = []
# ip_acc = []

# for cme_speed in cme_speeds:
    
#     fixed_dur = 8 #+ 2*cme_speed/600

#     cme = H.ConeCME(t_launch=1*u.day, longitude=360*u.deg, latitude = 0*u.deg, 
#                     initial_height = r_inner,
#                     width=width*u.deg, v=cme_speed* (u.km/u.s), thickness = 0*u.solRad,
#                     cme_fixed_duration = True, fixed_duration = fixed_dur*60*60*u.s)
#     model.solve([cme])
    
    
#     #get the transit time
#     stats = model.cmes[0].compute_arrival_at_body('EARTH')
#     tt = stats['t_transit'].value
#     transit_times.append(tt)
#     arrival_time = stats['t_arrive']
    
#     #find the arrival speed as the max withing 1 day of arrival
#     earth_series = HA.get_observer_timeseries(model, observer = 'Earth')
#     mask = (Time(earth_series['time']) >= arrival_time) & (
#         Time(earth_series['time']) <= arrival_time + 3*u.day)
    
#     v_1au = earth_series.loc[mask, 'vsw'].max()
#     arrival_speeds.append(v_1au)
    
#     #compute hte interplanetary acceleration
#     ip_acc.append(1000* (v_1au - cme_speed) / (tt*24*60*60) )


# ax1.plot(transit_times, cme_speeds,  linewidth = linewidth,
#          label = r'$\alpha$ = ' +str(width/2) +r'$^\circ$', color ='r')

# ax2.plot(arrival_speeds, cme_speeds,  linewidth = linewidth,
#          label = r'$\alpha$ = '  +str(width/2) +r'$^\circ$',color = 'r')

# ax3.plot(ip_acc, cme_speeds, linewidth = linewidth,
#          label = r'$\alpha$ = ' + str(width/2) +r'$^\circ$', color = 'r')

# i = i+1
