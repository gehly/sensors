import numpy as np
import sys
from datetime import datetime, timedelta
import copy
import time
import csv
import os

metis_dir = r'C:\Users\sgehly\Documents\code\metis'
sys.path.append(metis_dir)

from sensors import sensors as sens
from sensors import measurement_functions as mfunc
from utilities import astrodynamics as astro
from utilities import coordinate_systems as coord
from utilities import eop_functions as eop
from utilities import tle_functions as tle
from utilities.constants import GME, Re, wE


def build_tle_dict(obj_id_list=[], UTC_list=[]):
        
    tle_dict, tle_df = tle.get_spacetrack_tle_data(obj_id_list, UTC_list)
    
    return tle_dict


def build_geo_catalog(UTC_list, tle_dict={}):
    
    
    start = time.time()
    
    # # Generate TLE dictionary    
    # obj_id_list2 = [40146]
    # UTC_list2 = [datetime(2023, 8, 1, 0, 0, 0)]
    
    if len(tle_dict) == 0:
        tle_dict = build_tle_dict(UTC_list=UTC_list)
    
    retrieve_tle_time = time.time() - start

    # Define GEO params
    
    # # Space-track definitions
    # n_bounds = [0.99, 1.01] # rev/day
    # e_bounds = [0., 0.01]
    
    # # Flohrer definition (Dissertation Section 3.2.3)
    # n_bounds = [0.9, 1.1] # rev/day
    # e_bounds = [0., 0.2]
    # i_bounds = [0., 30.] # deg
    
    # # Siminski (Dissertation Eq 3.3)
    # a_bounds = [30000., 50000.]
    # e_bounds = [0., 0.3]
    # i_bounds = [0., 30.]
    
    # ESA Classification of Geosynchronous objects (2018) - for DISCOS
    # Extended GEO includes an inclination bound at [0, 25] deg
    # Inclined GEO then covers [25, 180]
    # a_bounds = [37948., 46380.]  # km
    # e_bounds = [0., 0.25]
    # i_bounds = [0., 180.]  # deg
    
    # # GEO slot (Holzinger 2016 - ish)
    Rgeo = 42164.2    
    a_bounds = [Rgeo - 50., Rgeo + 50.]
    e_bounds = [0., 0.01]
    i_bounds = [0., 15.]
    
    
    
    # Process and downselect to fit within bounds
    tle_dict2 = {}
    for obj_id in tle_dict:
        line2 = tle_dict[obj_id]['line2_list'][0]
        elem = tle.parse_tle_line2(line2)
        
        a = elem[0]
        e = elem[1]
        i = elem[2]
        
        if ((a >= a_bounds[0]) and (a <= a_bounds[1]) and 
            (e >= e_bounds[0]) and (e <= e_bounds[1]) and 
            (i >= i_bounds[0]) and (i <= i_bounds[1])):
            
            tle_dict2[obj_id] = tle_dict[obj_id]    
            
    start = time.time()
    
    # Generate all rotation matrices needed for time
    eop_alldata_text = eop.get_celestrak_eop_alldata()
    GCRF_TEME_list, ITRF_GCRF_list = \
        eop.batch_eop_rotation_matrices(UTC_list, eop_alldata_text)
        
    eop_time = time.time() - start
    
    start = time.time()
        
    obj_id_list = list(tle_dict2.keys())
    catalog = tle.propagate_TLE(obj_id_list, UTC_list, tle_dict=tle_dict2,
                                frame_flag=False)
    
    tle_prop_time = time.time() - start
    
    start = time.time()
    
    # Loop over times and convert to GCRF/ITRF
    for UTC in UTC_list:
        ii = UTC_list.index(UTC)
        
        # Retrieve rotation data
        GCRF_TEME = GCRF_TEME_list[ii]
        ITRF_GCRF = ITRF_GCRF_list[ii]
        
        # Loop over objects in catalog
        for obj_id in catalog:
            
            # Initialize
            if ii == 0:
                catalog[obj_id]['r_GCRF'] = []
                catalog[obj_id]['v_GCRF'] = []
                catalog[obj_id]['r_ITRF'] = []
            
            # Retrieve data and rotate coordinates to desired frames
            r_TEME = catalog[obj_id]['r_TEME'][ii]
            v_TEME = catalog[obj_id]['v_TEME'][ii]
            r_GCRF = np.dot(GCRF_TEME, r_TEME)
            v_GCRF = np.dot(GCRF_TEME, v_TEME)
            r_ITRF = np.dot(ITRF_GCRF, r_GCRF)
            
            # Note that v_ITRF transformation is more complex and can't be done
            # with simple matrix rotation, but also not needed for most cases
            
            # Store output
            catalog[obj_id]['r_GCRF'].append(r_GCRF)
            catalog[obj_id]['v_GCRF'].append(v_GCRF)
            catalog[obj_id]['r_ITRF'].append(r_ITRF)
            
            
    frame_rotate_time = time.time() - start
    
    print('retrieve_tle_time', retrieve_tle_time)
    print('eop_time', eop_time)
    print('tle_prop_time', tle_prop_time)
    print('frame_rotate_time', frame_rotate_time)
            

    return catalog


def catalog_to_csv(catalog, csv_file):
    
    
    obj_id_list = sorted(list(catalog.keys()))
    UTC_list = catalog[obj_id_list[0]]['UTC']
    
    with open(csv_file, 'w', newline='') as f:
        writer = csv.writer(f)
        
        # Header
        writer.writerow(['UTC', 'NORAD', 'x [GCRF,km]', 'y [GCRF,km]', 'z [GCRF,km]', 'vx [GCRF,km/s]', 'vy [GCRF,km/s]', 'vz [GCRF,km/s]'])
        
        for ii in range(len(UTC_list)):
            
            UTC = UTC_list[ii]
            UTC_str = UTC.strftime('%Y-%m-%d %H:%M:%S.%f')
            
            for obj_id in obj_id_list:
                
                x_GCRF = float(catalog[obj_id]['r_GCRF'][ii][0])
                y_GCRF = float(catalog[obj_id]['r_GCRF'][ii][1])
                z_GCRF = float(catalog[obj_id]['r_GCRF'][ii][2])
                dx_GCRF = float(catalog[obj_id]['v_GCRF'][ii][0])
                dy_GCRF = float(catalog[obj_id]['v_GCRF'][ii][1])
                dz_GCRF = float(catalog[obj_id]['v_GCRF'][ii][2])
                
                writer.writerow([UTC_str, str(obj_id), str(x_GCRF), str(y_GCRF), str(z_GCRF), str(dx_GCRF), str(dy_GCRF), str(dz_GCRF)])
            
    
    
    
    return



if __name__ == '__main__':
    
    # obj_id_list = [40146]
    # UTC_list = [datetime(2023, 9, 15, 0, 0, 0)]
    
    UTC_list = []
    
    tle_dict = build_tle_dict(UTC_list=UTC_list)
    
    print('num objects', len(tle_dict))
    
    datadir = '../data/2023_10_14_geo_catalog'
    
    for hr in range(20,24):
        fname = os.path.join(datadir, 'catalog_geo_slot_2023_10_14_' + str(hr).zfill(2) + 'h_' + str(hr+1).zfill(2) + 'h.csv')
        
        UTC0 = datetime(2023, 10, 14, hr, 0, 0)
        dt_vec = np.arange(0., 3600., 10.)
        UTC_list = [UTC0 + timedelta(seconds=dt) for dt in dt_vec]
        
        catalog = build_geo_catalog(UTC_list, tle_dict)
        
        catalog_to_csv(catalog, fname)
    
    
    
    
    