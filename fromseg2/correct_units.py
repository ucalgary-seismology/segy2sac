""" Correct previously converted sac files for seg2 header deascaling factor so unit is in m/s"""

import os
import sys
from glob import glob
from obspy import read
import logging

Logger = logging.getLogger(__name__)
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s\t%(name)s\t%(levelname)s\t%(message)s")

DESCALING_FACTOR = 0.00016985
STACK = 1.0


def main(data_dir):
    flist = glob(os.path.join(data_dir, "*.sac"))
    for file in flist:
        if file.find("unit") == -1:
            Logger.info("Correcting units for file: %s" % file)            
            try:
                tr = read(file, format="SAC")[0]
                tr.data = tr.data * DESCALING_FACTOR / STACK * 1e-3
                new_fname = file.replace(".sac", "_unitmps.sac")
                tr.write(new_fname, format="SAC")
                Logger.info("Wrote new file: %s" % new_fname)
                os.remove(file)
                Logger.info("Remove previous file: %s" % file)
            except:
                Logger.error("Something wrong happened while trying to process file %s" % file)
        else:
            Logger.info("File already in units m/s: %s" % file)


if __name__ == '__main__':
    folder = sys.argv[1]
    data_dir_root = sys.argv[2]
    data_dir = os.path.join(data_dir_root, folder)
    Logger.info("Processing files in folder: %s" % data_dir)

    main(data_dir)
