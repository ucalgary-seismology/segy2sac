from obspy import Stream, UTCDateTime
from obspy.io.segy.segy import iread_segy
import os
import numpy as np
import pandas as pd
from pyproj import Proj
from obspy.io.sac import SACTrace
import sys


# Constant parameters
NETWORK_CODE = "8O"
CHANNEL_CODE = "DP"  # Band and instrument code for channel naming
HAWK_SENSITIVITY = 22.5  # Geophone sensitivity V/m/s
HAWK_GAIN = 24.0  # Pre-amplifier gain in dB
INSTRUMENT_CORRECTION_SCALAR = HAWK_SENSITIVITY * (10. ** (HAWK_GAIN / 20))  # See 7.7 from Eaton 2018
SPIKE_THRESHOLD = 1e5  # Amplitude threshold for spikes. Anything above replaced with 0's
LOCAL_DECLINATION = 12  # Local magnetic declination in degrees
LOCAL_ELEVATION = 779.  # Local elevation in meters


def get_station_info(sgy_filepath, station_file):
    '''
    Get station information from CSV file for corresponding segy filepath.
    :param sgy_filepath: Full path of segy file
    :param station_file: CSV file containing station information
    :return: station name, latitude, longitude, elevation
    '''

    fname = os.path.split(sgy_filepath)[1]  # Get filename only without path

    # Get station info from CSV file:
    stainfo = pd.read_csv(station_file)
    station = str(stainfo.loc[stainfo['segy_filename'] == fname, 'StationName'].values[0])
    xutm = stainfo.loc[stainfo['segy_filename'] == fname, 'xUTM'].values[0]
    yutm = stainfo.loc[stainfo['segy_filename'] == fname, 'yUTM'].values[0]

    # Convert UTM to lat,long using WGS84 datum
    myProj = Proj("+proj=utm +zone=12U, +north +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
    stlo, stla = myProj(float(xutm), float(yutm), inverse=True)
    stel = LOCAL_ELEVATION  # Local elevation at CaMI.FRS

    return station, stla, stlo, stel


def save_to_file(st, station_file, segy_filepath, output_dir):
    '''
    Save trace to SAC file, filling header with station information
    :param st: Obspy stream with only one trace, the one to save
    :param station_file: station CSV filename (full path) with station info
    :param segy_filepath: SEGY filepath to find matching station information
    :param output_dir: output directory where to write SAC files.
    :return:
    '''

    # Get station info
    station, stla, stlo, stel = get_station_info(segy_filepath, station_file=station_file)

    # Output directory for station
    raw_out_dir = output_dir + station + '/'
    if not os.path.exists(raw_out_dir):
        os.mkdir(raw_out_dir)

    # Correct for spikes (amplitudes ~ 1e20, in general > 1)
    trace = st[0]
    data = trace.data
    data[data > SPIKE_THRESHOLD] = 0
    data[data < -SPIKE_THRESHOLD] = 0

    # Correct for instrument response with a scalar (corrects for amplitudes above corner frequency)
    data = data / INSTRUMENT_CORRECTION_SCALAR  # Convert amplitude to m/s
    trace.data = data

    # Create SAC trace
    sac = SACTrace.from_obspy_trace(trace)

    # Get channel name and orientation
    comp = trace.meta.segy.trace_header.trace_sequence_number_within_line
    if comp == 1:  # Vertical channel
        channel = CHANNEL_CODE + 'Z'
        sac.cmpaz = 0
        sac.cmpinc = 0
    elif comp == 2:  # Magnetic North-South
        channel = CHANNEL_CODE + '1'
        sac.cmpaz = LOCAL_DECLINATION
        sac.cmpinc = 90
    elif comp == 3:  # Magnetic East-West
        channel = CHANNEL_CODE + '2'
        sac.cmpaz = LOCAL_DECLINATION + 90
        sac.cmpinc = 90
    sac.kcmpnm = channel

    # Fill in header with other station info
    sac.kstnm = station
    sac.stla = stla
    sac.stlo = stlo
    sac.stel = stel
    sac.khole = ""
    sac.knetwk = NETWORK_CODE

    # Create output file name:
    tstart = trace.meta.starttime
    tend = trace.meta.endtime
    raw_fname = os.path.join(raw_out_dir,
                             "{0}.{1}..{2}.D.{3}.{4}_unitmps.sac".format(NETWORK_CODE, station, channel,
                                                                         tstart.strftime("%Y%m%d.%H%M%S"),
                                                                         tend.strftime("%H%M%S")))
    # Check if file already exists. If yes don't save
    if os.path.exists(raw_fname):
        print("%s already exists." % raw_fname)
        return

    # Now save SAC data to file
    print('Saving data to file: ' + raw_fname)
    sac.write(raw_fname)


def add_trace(trace_to_add, day_stream, previous_date, station_file, segy_filepath, output_dir):
    '''
    Add Obspy trace to stream. If the trace to add has a different date, then save the day stream, otherwise just merge.
    :param trace_to_add: trace to add to day stream
    :param day_stream: Obspy stream containing data for one day
    :param previous_date: UTCDateTime date corresponding to day stream
    :param station_file: station CSV filename (full path) with station info
    :param segy_filepath: SEGY filepath to find matching station information
    :param output_dir: output directory where to write SAC files.
    :return:
    '''

    current_date = trace_to_add.stats.endtime.date

    # Check if changing day compared to previous trace, if yes save previous day to SAC
    if current_date != previous_date:

        print('Changing date from ' + str(previous_date) + ' to ' + str(current_date))
        ts = trace_to_add.stats.starttime
        te = trace_to_add.stats.endtime
        dayend = UTCDateTime(te.year, te.month, te.day, 0, 0, 0)

        # Check if date changes within the trace itself
        if ts.day != te.day:

            # Split trace where date changes
            trace_previous_day = trace_to_add.copy()
            trace_previous_day.trim(starttime=ts, endtime=dayend)
            trace_next_day = trace_to_add.copy()
            trace_next_day.trim(starttime=dayend, endtime=te)

            # Merge the trace to complete previous day
            day_stream += trace_previous_day
            day_stream.merge(method=1, fill_value='interpolate', interpolation_samples=-1)
            # save to file
            save_to_file(day_stream, station_file, segy_filepath, output_dir)

            # Now create new stream for next day and merge with other trimmed trace
            day_stream = Stream()  # re-initialize
            day_stream += trace_next_day
            day_stream.merge(method=1, fill_value='interpolate', interpolation_samples=-1)
        else:
            save_to_file(day_stream, station_file, segy_filepath, output_dir)
            day_stream = Stream()
            day_stream += trace_to_add
            day_stream.merge(method=1, fill_value='interpolate', interpolation_samples=-1)

        # update previous date value
        previous_date = current_date

    else:
        # Continue merging trace
        day_stream += trace_to_add
        day_stream.merge(method=1, fill_value='interpolate', interpolation_samples=-1)

    return day_stream, previous_date


def main(sgy_file, station_file, output_dir_base):
    print('PROCESSING segy file ' + sgy_file)

    # Find starting date of SEGY file and initialize streams
    for trace in iread_segy(sgy_file):
        previous_date1 = trace.stats.endtime.date
        break
    previous_date2 = previous_date1
    previous_date3 = previous_date1
    comp1 = Stream()
    comp2 = Stream()
    comp3 = Stream()

    # Loop over traces and merge
    for trace in iread_segy(sgy_file):

        # Get date and channel of current trace
        current_comp = trace.stats.segy.trace_header.trace_sequence_number_within_line

        # Check if there are any gaps. If yes, fill with 0's
        if isinstance(trace.data, np.ma.masked_array):
            trace.data = trace.data.filled()

        # 1st component
        if current_comp == 1:  # DPZ

            comp1, previous_date1 = add_trace(trace, day_stream=comp1, previous_date=previous_date1,
                                              station_file=station_file,
                                              segy_filepath=sgy_file,
                                              output_dir=output_dir_base)

        # 2nd component
        elif current_comp == 2:  # DP1

            comp2, previous_date2 = add_trace(trace, day_stream=comp2, previous_date=previous_date2,
                                              station_file=station_file,
                                              segy_filepath=sgy_file,
                                              output_dir=output_dir_base)

        # 3rd component
        elif current_comp == 3:  # DP2

            comp3, previous_date3 = add_trace(trace, day_stream=comp3, previous_date=previous_date3,
                                              station_file=station_file,
                                              segy_filepath=sgy_file,
                                              output_dir=output_dir_base)

        else:  # Something's wrong
            raise NameError('Unrecognized sequence number')


if __name__ == '__main__':

    if len(sys.argv) < 4:
        raise ValueError(
            "Missing output arguments! Needs 3 arguments: \n python process_segy.py $segy_filepath $station_filepath $output_directory")

    # SEGY filename is given as input argument (full path)
    segy_filepath = sys.argv[1]

    # CSV file with station information for each segy file
    station_filepath = sys.argv[2]

    # Parent directory where to put data
    out_dir_parent = sys.argv[3]

    main(sgy_file=segy_filepath, station_file=station_filepath, output_dir_base=out_dir_parent)
