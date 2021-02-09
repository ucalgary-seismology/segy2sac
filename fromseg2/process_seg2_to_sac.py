import numpy as np
from obspy.io.sac import SACTrace
from obspy import Stream, read, UTCDateTime
import os
from glob import glob
import sys
import warnings


SENSITIVITY = 19.69  # Sensitivity of borehole geodes GS 32-CT in V/m/s


def get_station_component(n):
    """ Return the station code, channel and depth from the channel number """
    if n < 25:
        component = "DP1"
        n2 = n
    elif 24 < n < 49:
        component = "DP2"
        n2 = n - 24
    elif n > 48:
        component = "DP3"
        n2 = n - 48
    else:
        raise ValueError("Unknown channel code %d." % n)

    ntr = abs(n2 - 24) + 1
    station = "BH%03d" % ntr

    # trace 24 = 190.65 m below surface, trace 1 = 305.65 m = 190.65 + 5 * (24-1)
    depth = 305.65 - 5 * (24 - ntr)

    return station, component, depth


def process_seg2(flist, cha_list, outdir, day):
    """ Given a list of 30-s long SEG2 files, merge into one file per channel, decimate to 100 Hz and write as SAC"""
    
    # Create output directory with name of date for start time
    day_dir = os.path.join(outdir, day)
    if not os.path.exists(day_dir):
        os.mkdir(day_dir)

    # Loop over channels
    for channel in cha_list:

        sta, comp, dep = get_station_component(channel)
        print("Assembling data for channel %d: geode %s, component %s, depth %f m" % (channel, sta, comp, dep))
#        if glob(os.path.join(day_dir, "8O.%s..%s.*.sac" % (sta, comp))):
#            print("Station %s already processed." % sta)
#            continue
                                                
        # Initialize a Stream object
        st = Stream()

        # Loop over 30-s long files and append data
        for file in sorted(flist):
            tmp = read(file, format="SEG2")
            if not tmp:
                raise ValueError("Data file read is empty: %s" % file)
            trace = []
            for tr in tmp:
                if int(tr.stats.seg2.CHANNEL_NUMBER) == channel:
                    
                    # Check for NaN
                    if np.isnan(tr.data).any():
                        print("There are NaN values in this file: %s" % file)
                        idx_nan = np.argwhere(np.isnan(tr.data))
                        # Replace by 0
                        tr.data[idx_nan] = 0
                        print("Replaced %d samples by 0." % len(idx_nan))
                        continue
                        
                    trace = tr
                    # correct start_time for SEG-2 acquisition delay
                    trace.stats.starttime = tr.stats.starttime + float(tr.stats.seg2.DELAY)
                    # correct amplitude for amplifier-gain
                    amp_gain = float(tr.stats.seg2.FIXED_GAIN.split(" ")[0])
                    instrument_correction_scalar = SENSITIVITY * (10.**(amp_gain/20))  # See 7.7 from Eaton 2018
                    trace.data = trace.data * float(tr.stats.seg2.DESCALING_FACTOR)/float(tr.stats.seg2.STACK) * 1e-3  # convert to Volt
                    trace.data = trace.data / instrument_correction_scalar
                    break
            if not trace:
                print("Could not find data for this channel in file %s" % file)
                continue
            st.append(trace)
            st.merge(method=1, interpolation_samples=-1, fill_value=0)

        # Check stream has only one merged trace
        if len(st) == 0:
            raise ValueError("Something went wrong: no data found for channel %d" % channel)
        elif len(st) > 1:
            print(st)
            raise ValueError("Something went wrong: there should be only one trace in the stream.")
        
        # Resample to 500Hz
        #st.resample(500.0)
        
        # Check data for gaps
        if st.get_gaps():
            st.print_gaps()
            warnings.warn("There are gaps in the data.")

        # Decimate data by factor of 5: from 500Hz to 100Hz
        #st.decimate(factor=5)

        # Cut day
        ts = UTCDateTime(int(day[0:4]),int(day[4:6]),int(day[6:]),0,0,0)
        te = ts + 3600*24
        st.trim(starttime=ts,endtime=te)
        print(st)

        # Add header info
        sac = SACTrace.from_obspy_trace(st[0])
        sac.kcmpnm = comp
        sac.kstnm = sta
        sac.stel = 779.60  # Ground elevation in m
        sac.stdp = dep
        sac.stla = 50.45031  # coordinates of geophysics well
        sac.stlo = -112.12087  # coordinates of geophysics well
        sac.knetwk = "8O"  # give network code BH for "borehole"

        # Define output file name
        sacname = "8O.%s..%s.%s_%s_%dHz_unitmps.sac" % (sta, comp,
                                                 st[0].stats.starttime.strftime("%Y%m%d%H%M%S"),
                                                 st[0].stats.endtime.strftime("%Y%m%d%H%M%S"), int(st[0].stats.sampling_rate))

        sac.write(os.path.join(day_dir, sacname))


def main(output_dir, parent_dir, folder_list, channel_list, day):

    for folder in folder_list:
        print("Processing data in directory %s" % folder)

        # Get list of SEG2 files to process in that folder
        file_list = glob(os.path.join(parent_dir, folder, "*.dat"))
        if not file_list:
            print("No data found in directory %s." % folder)
            continue

        # Process SEG2 to SAC
        process_seg2(file_list, channel_list, output_dir, day)


if __name__ == "__main__":

        # Define directory where to write SAC output files:
    output_dir = "/home/gilbert_lab/cami_frs/borehole_data/sac_daily_500Hz/"
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
        print("Created output directory: %s" % output_dir)

    # Define the parent directory containing the day sub-directories with SEG2 data:
    parent_dir = "/home/gilbert_lab/cami_frs/borehole_data/raw_seg2"
    if not os.path.exists(parent_dir):
        raise ValueError("Parent data directory does not exists!")

    # Get list of day sub-directories
    folder_list = [sys.argv[1]]
    #if not folder_list:
    #    raise ValueError("List of SEG-2 data subdirectories is empty!")
    
    day = sys.argv[1]

    # List of channels to extract:
    # If all 72 channels desired:
    #channel_list = np.arange(1, 73, 1).tolist()
#    channel_list = [7, 31, 55,  # G18 # Each line corresponds to 3 channels of one geode
#                    8, 32, 56,  # G17
#                    9, 33, 57,  # G16
#                    10, 34, 58,  # G15
#                    11, 35, 59,  # G14
#                    12, 36, 60,  # G13
#                    13, 37, 61,  # G12
#                    14, 38, 62,  # G11
#                    17, 41, 65,  # G8
#                    19, 43, 67,  # G6
#                    20, 44, 68,  # G5
#                    21, 45, 69,  # G4
#                    23, 47, 71,  # G2
#                    24, 48, 72]  # G1 (shallowest geode)
    channel_list = [18, 42, 66] # G7
#                    15, 39, 63, # G10
#                    4, 28, 52, # G21
#                    6, 30, 54, # G19
#                    22, 46, 70] # G3
        
    if not channel_list:
        raise ValueError("List of channels is empty!")

    
    
    main(output_dir, parent_dir, folder_list, channel_list, day)

