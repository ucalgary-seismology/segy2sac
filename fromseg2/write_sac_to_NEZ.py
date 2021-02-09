""" Read SAC files in 1,2,3 system, correct to N, E, Z and write to new SAC files"""

import pandas as pd
from obspy import read, Stream
from obspy.io.sac import SACTrace
import os
from helper_functions import rotate_to_zne
import sys
from glob import glob


def main(data_dir, output_dir):

    # Get borehole geode info
    stats = pd.read_csv("walkaround_borehole_orientation_stats.csv")

    # Loop over stations/geodes
    for sta in ['BH001', 'BH002', 'BH004', 'BH005', 'BH006', 'BH008', 'BH011', 'BH012', 'BH013', 'BH014', 'BH015', 'BH016', 'BH017', 'BH018']:  # geodes.Geode:

        print("Processing station %s" % sta)
        # Rotation matrix for correction
        inc_corr = - stats.loc[stats['geode'] == sta, 'inclination'].values[0]
        baz_corr = - stats.loc[stats['geode'] == sta, 'mean'].values[0]

        # Read data
        flist1 = glob(os.path.join(data_dir, "8O.%s..DP1*" % sta))
        flist2 = glob(os.path.join(data_dir, "8O.%s..DP2*" % sta))
        flist3 = glob(os.path.join(data_dir, "8O.%s..DP3*" % sta))
        if not flist1 or not flist2 or not flist3:
            continue
            
        # Look for suffix
        file = os.path.split(flist1[0])[1]
        if file.find("unit") != -1:
            idx = file.find("unit")
            suffix = "_" + file[idx:].split(".")[0]
        else:
            suffix = ""
        
        st = Stream()
        st += read(os.path.join(data_dir, "8O.%s..DP1*" % sta), format="SAC")[0]
        st += read(os.path.join(data_dir, "8O.%s..DP2*" % sta), format="SAC")[0]
        st += read(os.path.join(data_dir, "8O.%s..DP3*" % sta), format="SAC")[0]
        
        if len(st) == 0:
            continue
            
        st.resample(500.0)
        starttime = st[0].stats.starttime
        endtime = st[0].stats.endtime

        # Save rotated coordinates
        st3 = st.copy()
        compN, compE, compZ = rotate_to_zne(st3[0].data, st3[1].data, st3[2].data, inc_corr, baz_corr)
        trN = SACTrace.from_obspy_trace(st3[0])
        trN.data = compN
        trN.kcmpnm = "DPN"
        trN.write(os.path.join(output_dir,
                               "8O.%s..DPN.%s_%s%s.sac" % (sta, starttime.strftime("%Y%m%d%H%M%S"), endtime.strftime("%Y%m%d%H%M%S"), suffix)))
        trE = SACTrace.from_obspy_trace(st3[1])
        trE.data = compE
        trE.kcmpnm = "DPE"
        trE.write(os.path.join(output_dir,
                               "8O.%s..DPE.%s_%s%s.sac" % (sta, starttime.strftime("%Y%m%d%H%M%S"), endtime.strftime("%Y%m%d%H%M%S"), suffix)))
        trZ = SACTrace.from_obspy_trace(st3[2])
        trZ.data = compZ
        trZ.kcmpnm = "DPZ"
        trZ.write(os.path.join(output_dir,
                               "8O.%s..DPZ.%s_%s%s.sac" % (sta, starttime.strftime("%Y%m%d%H%M%S"), endtime.strftime("%Y%m%d%H%M%S"), suffix)))


if __name__ == '__main__':
    folder = sys.argv[1]
    #data_dir = "/home/gilbert_lab/cami_frs/borehole_data/sac_hourly_125Hz/%s" % folder
    data_dir = "/home/gilbert_lab/cami_frs/borehole_data/sac_daily_500Hz/%s" % folder
    output_dir = "/home/gilbert_lab/cami_frs/borehole_data/sac_daily_nez_500Hz/%s" % folder
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    main(data_dir, output_dir)
