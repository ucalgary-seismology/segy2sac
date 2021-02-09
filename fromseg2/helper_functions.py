""" Collection of helper functions """

import numpy as np
from scipy.linalg import lstsq
import matplotlib.pyplot as plt
import math
from scipy.spatial.transform import Rotation as R
from obspy.core.inventory import Inventory, Network, Station, Channel, Site
from obspy.realtime.signal import kurtosis
from obspy import UTCDateTime
from statsmodels import robust

def mccc(seis, dt, twin, ccmin, comp='Z'):
    """ FUNCTION [TDEL,RMEAN,SIGR] = MCCC(SEIS,DT,TWIN);
    Function MCCC determines optimum relative delay times for a set of seismograms based on the
    VanDecar & Crosson multi-channel cross-correlation algorithm. SEIS is the set of seismograms.
    It is assumed that this set includes the window of interest and nothing more since we calculate the
    correlation functions in the Fourier domain. DT is the sample interval and TWIN is the window about
    zero in which the maximum search is performed (if TWIN is not specified, the search is performed over
    the entire correlation interval).
    APP added the ccmin, such that only signals that meet some threshold similarity contribute to the delay times. """

    # Set nt to twice length of seismogram section to avoid
    # spectral contamination/overlap. Note we assume that
    # columns enumerate time samples, and rows enumerate stations.
    # Note in typical application ns is not number of stations...its really number of events
    # all data is from one station
    nt = np.shape(seis)[1] * 2
    ns = np.shape(seis)[0]
    tcc = np.zeros([ns, ns])

    # Copy seis for normalization correction
    seis2 = np.copy(seis)

    # Set width of window around 0 time to search for maximum
    # mask = np.ones([1,nt])
    # if nargin == 3:
    itw = int(np.fix(twin / (2 * dt)))
    mask = np.zeros([1, nt])[0]
    mask[0:itw + 1] = 1.0
    mask[nt - itw:nt] = 1.0

    # Zero array for sigt and list on non-zero channels
    sigt = np.zeros(ns)

    # First remove means, compute autocorrelations, and find non-zeroed stations.
    for iss in range(0, ns):
        seis[iss, :] = seis[iss, :] - np.mean(seis[iss, :])
        ffiss = np.fft.fft(seis[iss, :], nt)
        acf = np.real(np.fft.ifft(ffiss * np.conj(ffiss), nt))
        sigt[iss] = np.sqrt(max(acf))

    # Determine relative delay times between all pairs of traces.
    r = np.zeros([ns, ns])
    tcc = np.zeros([ns, ns])

    # Two-Channel normalization ---------------------------------------------------------

    # This loop gets a correct r by checking how many channels are actually being compared
    if comp == 'NE':

        # First find the zero-channels (the np.any tool will fill in zeroNE)
        # zeroNE ends up with [1,0], [0,1], or [1,1] for each channel, 1 meaning there IS data
        zeroNE = np.zeros([ns, 2])
        dum = np.any(seis2[:, 0:nt / 4], 1, zeroNE[:, 0])
        dum = np.any(seis2[:, nt / 4:nt / 2], 1, zeroNE[:, 1])

        # Now start main (outer) loop
        for iss in range(0, ns - 1):
            ffiss = np.conj(np.fft.fft(seis[iss, :], nt))

            for jss in range(iss + 1, ns):

                ffjss = np.fft.fft(seis[jss, :], nt)
                # ccf  = np.real(np.fft.ifft(ffiss*ffjss,nt))*mask
                ccf = np.fft.fftshift(np.real(np.fft.ifft(ffiss * ffjss, nt)) * mask)
                cmax = np.max(ccf)

                # chcor for channel correction sqrt[ abs( diff[jss] - diff[iss]) + 1]
                # This would be perfect correction if N,E channels always had equal power, but for now is approximate
                chcor = np.sqrt(abs(zeroNE[iss, 0] - zeroNE[jss, 0] - zeroNE[iss, 1] + zeroNE[jss, 1]) + 1)

                # OLD, INCORRECT chcor
                # chcor = np.sqrt( np.sum(zeroNE[iss,:])+np.sum(zeroNE[jss,:]) - (zeroNE[iss,0]*zeroNE[jss,0]+zeroNE[iss,1]*zeroNE[jss,1]) )

                rtemp = cmax * chcor / (sigt[iss] * sigt[jss])

                # Quadratic interpolation for optimal time (only if CC found > ccmin)
                if rtemp > ccmin:

                    ttemp = np.argmax(ccf)

                    x = np.array(ccf[ttemp - 1:ttemp + 2])
                    A = np.array([[1, -1, 1], [0, 0, 1], [1, 1, 1]])

                    [a, b, c] = lstsq(A, x)[0]

                    # Solve dy/dx = 2ax + b = 0 for time (x)
                    tcc[iss, jss] = -b / (2 * a) + ttemp

                    # Estimate cross-correlation coefficient
                    # r[iss,jss] = cmax/(sigt[iss]*sigt[jss])
                    r[iss, jss] = rtemp
                else:
                    tcc[iss, jss] = nt / 2

                    # Reguar Normalization Version -------------------------------------------------------
    elif comp != 'NE':
        for iss in range(0, ns - 1):
            ffiss = np.conj(np.fft.fft(seis[iss, :], nt))
            for jss in range(iss + 1, ns):

                ffjss = np.fft.fft(seis[jss, :], nt)
                # ccf  = np.real(np.fft.ifft(ffiss*ffjss,nt))*mask
                ccf = np.fft.fftshift(np.real(np.fft.ifft(ffiss * ffjss, nt)) * mask)
                cmax = np.max(ccf)

                rtemp = cmax / (sigt[iss] * sigt[jss])

                # Quadratic interpolation for optimal time (only if CC found > ccmin)
                if rtemp > ccmin:

                    ttemp = np.argmax(ccf)

                    x = np.array(ccf[ttemp - 1:ttemp + 2])
                    A = np.array([[1, -1, 1], [0, 0, 1], [1, 1, 1]])

                    [a, b, c] = lstsq(A, x)[0]

                    # Solve dy/dx = 2ax + b = 0 for time (x)
                    tcc[iss, jss] = -b / (2 * a) + ttemp

                    # Estimate cross-correlation coefficient
                    # r[iss,jss] = cmax/(sigt[iss]*sigt[jss])
                    r[iss, jss] = rtemp
                else:
                    tcc[iss, jss] = nt / 2

                    #######################################################

    # Some r could have been made > 1 due to approximation, fix this
    r[r >= 1] = 0.99

    # Fisher's transform of cross-correlation coefficients to produce
    # normally distributed quantity on which Gaussian statistics
    # may be computed and then inverse transformed
    z = 0.5 * np.log((1 + r) / (1 - r))
    zmean = np.zeros(ns)
    for iss in range(0, ns):
        zmean[iss] = (np.sum(z[iss, :]) + np.sum(z[:, iss])) / (ns - 1)
    rmean = (np.exp(2 * zmean) - 1) / (np.exp(2 * zmean) + 1)

    # Correct negative delays (for fftshifted times)
    # ix = np.where( tcc>nt/2);  tcc[ix] = tcc[ix]-nt
    tcc = tcc - nt / 2

    # Subtract 1 to account for sample 1 at 0 lag (Not in python)
    # tcc = tcc-1

    # Multiply by sample rate
    tcc = tcc * dt

    # Use sum rule to assemble optimal delay times with zero mean
    tdel = np.zeros(ns)

    # I changed the tdel calculation to not include zeroed-out waveform pairs in normalization
    for iss in range(0, ns):
        ttemp = np.append(tcc[iss, iss + 1:ns], -tcc[0:iss, iss])
        tdel[iss] = np.sum(ttemp) / (np.count_nonzero(ttemp) + 1)
        # tdel[iss] = ( np.sum(tcc[iss,iss+1:ns])-np.sum(tcc[0:iss,iss]) )/ns

    # Compute associated residuals
    res = np.zeros([ns, ns])
    sigr = np.zeros(ns)
    for iss in range(0, ns - 1):
        for jss in range(iss + 1, ns):
            res[iss, jss] = tcc[iss, jss] - (tdel[iss] - tdel[jss])

    for iss in range(0, ns):
        sigr[iss] = np.sqrt((np.sum(res[iss, iss + 1:ns] ** 2) + np.sum(res[0:iss, iss] ** 2)) / (ns - 2))

    return tdel, rmean, sigr, r, tcc


def shift(seis, dt, t0):
    """ Produces shifted time series using FFT (sinc(x)) interpolation
    seis: input time series as a numpy array or a matrix with one row per trace
    dt: sampling interval
    t0: desired time shit (scalar or array of length equal to number of rows in seis.
    """

    # Definitions
    n1 = max(seis.shape)
    m1 = min(seis.shape)
    omega = np.arange(0, n1 - 1) * 2 * np.pi / (n1 * dt)  # fft frequencies
    seis_shifted = np.empty(seis.shape)

    # Loop through traces. Check for even/odd length since this affects definition of Nyquist frequency
    if n1 % 2 == 0:  # even length case

        for i1 in range(0, m1):
            ffn = np.fft.fft(seis[i1, :], n1)
            ffn[0:int(n1 / 2)] = ffn[0:int(n1 / 2)] * np.exp(-1j * omega[0:int(n1 / 2)] * t0[i1])
            ffn[int(n1 / 2)] = ffn[int(n1 / 2)] * np.cos(np.pi * t0[i1] / dt)

            # We can use negative frequencies only to get back in time domain
            # (faster just to use conjugate relation for real signal
            ffn[int(n1 / 2 + 1):] = np.conj(ffn[int(n1 / 2):1:-1])
            dum = np.real(np.fft.ifft(ffn, n1))
            seis_shifted[i1, :] = dum[0:n1]

    else:  # odd length case
        # NEEDS TO BE FIXED
        for i1 in range(0, m1):
            ffn = np.fft.fft(seis[i1, :], n1)
            ffn[0:int((n1 + 1) / 2)] = ffn[0:int((n1 + 1) / 2)] * np.exp(-1j * omega[0:int((n1 + 1) / 2)] * t0[i1])
            ffn[int((n1 + 1) / 2 + 1):] = np.conj(ffn[int((n1 + 1) / 2):1:-1])
            dum = np.real(np.fft.ifft(ffn, n1))
            seis_shifted[i1, :] = dum[0:n1]

    return seis_shifted


def section_plot(seis, times, aflag=-1, labels="", title=""):
    import numpy.matlib
    import matplotlib.dates as mdates

    nt = max(seis.shape)
    ns = min(seis.shape)
    seismat = np.empty(seis.shape)
    tmat = np.matlib.repmat(times, ns, 1)

    # Normalization
    if aflag < 0:
        for iy in range(0, ns):
            seismat[iy, :] = iy + 1 - 0.7 * seis[iy, :] / np.max(np.abs(seis) + 0.0000001)
    elif aflag == 0:
        for iy in range(0, ns):
            seismat[iy, :] = iy + 1 - 8.0 * seis[iy, :] / np.max(np.abs(seis[iy, :]) + 0.0000001)
    else:
        for iy in range(0, ns):
            seismat[iy, :] = iy + 1 - 1.0 * seis[iy, :] / aflag

    plt.rcParams['figure.figsize'] = [10, 10]
    fig, ax = plt.subplots(1, 1)
    for iy in range(0, ns):
        ax.plot_date(tmat[iy, :], seismat[iy, :], Marker=None, linestyle='-', color='black')
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
    ax.set_yticks(np.arange(1, ns + 1))
    ax.set_yticklabels(labels)
    ax.set_title(title)
    plt.show()


def rotate_to_zne(comp1, comp2, comp3, inc, baz):
    rotmat_inc = R.from_euler('y', inc, degrees=True).as_matrix()  # inclination correction
    rotmat_baz = R.from_euler('z', baz, degrees=True).as_matrix()  # backazimuth correction
    rotmat = np.dot(rotmat_baz, rotmat_inc)  # combined rotation matrix

    out = np.dot(rotmat, np.array([comp1, comp2, comp3]))
    compN = out[0, :]
    compE = out[1, :]
    compZ = out[2, :]

    return compN, compE, compZ


def hodogram(comp1, comp2, compz, title="", ndt=0.001, azimuth=None, incidence=None):
    fig, axs = plt.subplots(2, 2)
    # fig.set_size_inches(8, 8)
    axs[0][0].plot(comp1, comp2)
    axs[0][0].set_xlabel('horizontal 1')
    axs[0][0].set_ylabel('horizontal 2')

    axs[0][1].plot(comp1, compz)
    axs[0][1].set_xlabel('horizontal 1')
    axs[0][1].set_ylabel('vertical')

    axs[1][0].plot(comp2, compz)
    axs[1][0].set_xlabel('horizontal 2')
    axs[1][0].set_ylabel('vertical')

    # Set axes limits
    mmax = np.max(np.abs(np.array([comp1, comp2, compz])))
    axs[0][0].set_xlim(-mmax, mmax)
    axs[0][0].set_ylim(-mmax, mmax)
    axs[0][1].set_xlim(-mmax, mmax)
    axs[0][1].set_ylim(-mmax, mmax)
    axs[1][0].set_xlim(-mmax, mmax)
    axs[1][0].set_ylim(-mmax, mmax)

    # Add azimuth and incidence angle to hodogram
    if azimuth is not None:
        azi_rad = np.deg2rad(azimuth)
        x = mmax * np.sin(azi_rad)
        y = mmax * np.cos(azi_rad)
        axs[0][0].plot([x, -x], [y, -y], marker=None, linestyle="-", color="red", linewidth=2)
    if incidence is not None:
        inc_rad = np.deg2rad(incidence)
        x = mmax * np.sin(inc_rad)
        y = -mmax * np.cos(inc_rad)
        axs[0][1].plot([x, -x], [y, -y], marker=None, linestyle="-", color="red", linewidth=2)
        axs[1][0].plot([x, -x], [y, -y], marker=None, linestyle="-", color="red", linewidth=2)

    t = np.arange(0, len(comp1) * ndt, ndt)
    axs[1][1].plot(t, comp1, marker=None, color='blue')
    axs[1][1].plot(t, comp2, marker=None, color='green')
    axs[1][1].plot(t, compz, marker=None, color='red')
    if azimuth is not None and incidence is not None:
        plt.suptitle("Azimuth = %f, Incidence = %f" % (azimuth, incidence))
    else:
        plt.suptitle(title)


def polarization_svd(datax, datay, dataz):
    covmat = np.zeros([3, 3])
    covmat[0][0] = np.cov(datax, rowvar=False)
    covmat[0][1] = covmat[1][0] = np.cov(datax, datay, rowvar=False)[0, 1]
    covmat[0][2] = covmat[2][0] = np.cov(datax, dataz, rowvar=False)[0, 1]
    covmat[1][1] = np.cov(datay, rowvar=False)
    covmat[1][2] = covmat[2][1] = np.cov(dataz, datay, rowvar=False)[0, 1]
    covmat[2][2] = np.cov(dataz, rowvar=False)
    eigenvec, eigenval, unit_vec = (np.linalg.svd(covmat))
    azimuth = math.degrees(np.arctan2(eigenvec[0, 0] * np.sign(eigenvec[2, 0]), eigenvec[1, 0] * np.sign(eigenvec[2, 0])))
    incidence = math.degrees(np.arccos(eigenvec[2, 0]))
    return azimuth, incidence


def dip_azimuth2nez_base_vector(dip, azimuth):
    """
    Helper function converting a vector described with azimuth and dip of unit
    length to a vector in the NEZ (North, East, Vertical) base.

    The definition of azimuth and dip is according to the SEED reference
    manual.
    """
    dip = np.deg2rad(dip)
    azimuth = np.deg2rad(azimuth)

    return np.array([np.cos(azimuth) * np.cos(dip),
                     np.sin(azimuth) * np.cos(dip),
                     -np.sin(dip)])


def get_inventory(stations, depths, lat=50.45031, long=-112.12087, elevation=779.0, dip1=0, azi1=0, dip2=0, azi2=90, dip3=90, azi3=0):
    inv = Inventory(networks=[], source="Genevieve")
    net = Network(code="BH", stations=[], description=" ", start_date=UTCDateTime(2019, 1, 1))
    for i, station in enumerate(stations):
        dep = depths[i]
        sta = Station(code=station, latitude=lat, longitude=long, elevation=elevation,
                      creation_date=UTCDateTime(2019, 1, 1), site=Site(name="borehole"))
        chaz = Channel(code="DPZ", location_code="", latitude=lat, longitude=long, elevation=elevation,
                       azimuth=azi3, dip=dip3, depth=dep, sample_rate=500)
        cha1 = Channel(code="DPN", location_code="", latitude=lat, longitude=long, elevation=elevation,
                       azimuth=azi1, dip=dip1, depth=dep, sample_rate=500)
        cha2 = Channel(code="DPE", location_code="", latitude=lat, longitude=long, elevation=elevation,
                       azimuth=azi2, dip=dip2, depth=dep, sample_rate=500)
        sta.channels.append(chaz)
        sta.channels.append(cha1)
        sta.channels.append(cha2)
        net.stations.append(sta)
    inv.networks.append(net)
    return inv


def get_arrival_from_kurtosis(st):
    tarr = []
    for tr in st.traces:
        k = kurtosis(st[0], win=0.1)
        dk = np.gradient(k)
        ix = np.argmax(np.abs(dk))
        t = tr.times("utcdatetime")
        tarr.append(t[ix])
    return min(tarr)


def semblance_weighted_stack(seis1, seis2, seis3):
    """ seis: one trace per row"""

    # Calculate semblance for each component
    semblance1 = np.sum(np.sum(np.power(seis1, 2), axis=0)) / np.power(np.sum(np.sum(seis1, axis=0)), 2)
    semblance2 = np.sum(np.sum(np.power(seis2, 2), axis=0)) / np.power(np.sum(np.sum(seis2, axis=0)), 2)
    semblance3 = np.sum(np.sum(np.power(seis3, 2), axis=0)) / np.power(np.sum(np.sum(seis3, axis=0)), 2)

    # semblance weighted stack
    stack = semblance1 * np.sum(np.sum(seis1, axis=0)) + \
            semblance2 * np.sum(np.sum(seis2, axis=0)) + \
            semblance3 * np.sum(np.sum(seis3, axis=0))

    return stack


def pad_with_white_noise(tr, buffer=1):
    """Pad and fill in with uncorrelated white Gaussian noise with a buffer at each end of buffer seconds"""
    starttime = tr.stats.starttime - buffer
    endtime = tr.stats.starttime + buffer
    noise_ref_data = np.hstack(tr.data[0:50])
    median_gap = np.median(noise_ref_data)
    mad_gap = robust.mad(noise_ref_data)
    tr.trim(starttime=starttime, endtime=endtime, pad=True, nearest_sample=True, fill_value=999)
    ixnone = np.argwhere(tr.data == 999).flatten()
    tr.data[ixnone] = mad_gap * np.random.randn(len(ixnone)) + median_gap
    return tr
