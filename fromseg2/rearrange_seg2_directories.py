import os
from glob import glob
from obspy import read, UTCDateTime
from shutil import copyfile
import sys

indir = "/home/gilbert_lab/cami_frs/borehole_data/raw_segy/"
outdir = "/home/gilbert_lab/cami_frs/borehole_data/raw_seg2/"
folder = sys.argv[1]

flist = glob(os.path.join(indir,folder, "*.dat"))
flist.sort()
for i, file in enumerate(flist):
    try:
        tr = read(file, headonly=True, format="SEG2")[0]
    except:
        print("Problem reading file %s" % file)
        continue
    starttime = tr.stats.starttime + float(tr.stats.seg2.DELAY)
    endtime = tr.stats.endtime + float(tr.stats.seg2.DELAY)
    print("folder: %s, file %d/%d: %s, start time: %s" % (os.path.split(folder)[1], i, len(flist), os.path.split(file)[1], starttime.strftime("%Y-%m-%d %H:%M:%S")))
    startdate = starttime.strftime("%Y%m%d")
    enddate = endtime.strftime("%Y%m%d")
    if not os.path.exists(os.path.join(outdir, startdate)):
        os.makedirs(os.path.join(outdir, startdate))
    fname = "%s_%s" % (starttime.strftime("%Y%m%d%H%M%S"), os.path.split(file)[1])
    outfile1 = os.path.join(outdir, startdate, fname)
    os.rename(file, outfile1)
    print("Moving %s to \n\t %s" % (file, outfile1))
    if endtime.date != starttime.date:
        if not os.path.exists(os.path.join(outdir, enddate)):
            os.makedirs(os.path.join(outdir, enddate))
        outfile2 = os.path.join(outdir, enddate, fname)
        print("Moving %s to \n\t %s" % (file, outfile2))
        copyfile(outfile1, outfile2)
