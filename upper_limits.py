import numpy as np
import optparse
from stochpy_radiometer.utils import *
from stochpy_radiometer.doppler_shift import segment_list_shifts
from stochpy_radiometer.doppler_shift import segment_list_shifts_binary
from scipy.interpolate import interp1d
from scipy.io import loadmat, savemat
from multiprocessing import (Process, Queue as ProcessQueue)
from stochpy_radiometer.getUpperLimit_corrected import getUpperLimit_marg
from stochpy_radiometer.getUpperLimit import getUpperLimit
import pickle
import h5py


def parse_command_line():
    """
    parse command line
    """
    parser = optparse.OptionParser()
    parser.add_option(
            "--ini-file", "-i", help="params file",
            dest="ini", type=str)
    params, args = parser.parse_args()
    return params
    
def load_h5file(h5file):
    """
    Properly load hdf5 file from the output of Pystoch
    Specifically, we fill in the gaps of missing frequencies with zeros.

    Parameters
    ----------
    h5file : `str`
        hdf5 file from Pystoch

    Returns
    -------
    ptEst_long : `numpy.ndarray`
        point estimate with zeros for notched frequencies
    sigma_long : `numpy.ndarray`
        sigma with zeros for notched frequencies
    f : `numpy.ndarray`
        frequency with no skips
    df : `float`
        frequency bin width
    """
    nbr = h5py.File(h5file, "r")
    print(nbr.keys())
    ptEst = nbr['ptEst'][:]
    sigma = nbr['sig'][:]
    f = nbr['f_all'][:]
    nbr.close()    

    df = np.min(f[1:] - f[:-1]) # bin size
    # fix our frequency array, point estimate, and sigma
    # to fill in missing values
    fnews = np.arange(np.min(f), np.max(f)+df, df)


    fnews = np.arange(np.min(f), np.max(f)+df, df)
    ptEst_long = np.zeros(fnews.size)
    sigma_long = np.zeros(fnews.size)
    ct = 0 # counter for array that was cut up
    for ii,fnew in enumerate(fnews):
        if np.any(f == fnew):
            ptEst_long[ii] = ptEst[ct]
            sigma_long[ii] = sigma[ct]
            ct += 1 # iterate
        else:
            pass
    f = fnews
    return ptEst_long, sigma_long, f, df
    
def load_pkl(path):
    """
    WARNING LOADING PKL INSTEAD, NOT EDITING FUNCTION NAME FOR THE MOMENT
    Properly load hdf5 file from the output of Pystoch
    Specifically, we fill in the gaps of missing frequencies with zeros.

    Parameters
    ----------
    h5file : `str`
        hdf5 file from Pystoch

    Returns
    -------
    ptEst_long : `numpy.ndarray`
        point estimate with zeros for notched frequencies
    sigma_long : `numpy.ndarray`
        sigma with zeros for notched frequencies
    f : `numpy.ndarray`
        frequency with no skips
    df : `float`
        frequency bin width
    """

    with open(path, 'rb') as f:
        struct = pickle.load(f)
    #print(struct.keys())
    
    #ptEst = struct['ptEst'][:]
    #sigma = struct['sig'][:]
    #f = struct['f_all'][:]
    ptEst = np.array(struct['ptEst'])  
    sigma = np.array(struct['sig'])    
    f = np.array(struct['f_all'])      

    df = np.min(f[1:] - f[:-1]) # bin size
    # fix our frequency array, point estimate, and sigma
    # to fill in missing values
    fnews = np.arange(np.min(f), np.max(f)+df, df)


    fnews = np.arange(np.min(f), np.max(f)+df, df)
    ptEst_long = np.zeros(fnews.size)
    sigma_long = np.zeros(fnews.size)
    ct = 0 # counter for array that was cut up
    for ii,fnew in enumerate(fnews):
        if np.any(f == fnew):
            ptEst_long[ii] = ptEst[ct]
            sigma_long[ii] = sigma[ct]
            ct += 1 # iterate
        else:
            pass
    f = fnews
    return ptEst_long, sigma_long, f, df



def run_simulation(q, sigma, f, fl, fh, nsims):
    """
    run a simulation

    Parameters
    ----------
    q : TODO
    sigma : TODO
    f : TODO
    fl : TODO
    fh : TODO
    nsims : TODO

    Returns
    -------
    TODO

    """
    np.random.seed()
    maxSNRs = np.zeros(int(nsims))
    df = f[1] - f[0]
    for ii, sim in enumerate(range(int(nsims))):
        y = sigma * np.random.randn(sigma.size)
        newy, new_sigma, nc, nc_a, nc_b, r = combine_bins(
                y, sigma, f, fl, fh, df)
        snr = newy / new_sigma
        snr[np.isnan(snr)] = 0
        maxSNRs[ii] = np.max(snr)
    q.put(maxSNRs)
    

def do_pproc(analysis_params, binary_params):
    """
    Do post post processing for a source!
    Saves ``final_results.mat`` file in
    the output directory specified in
    analysis_params.

    Parameters
    ----------
    analysis_params : `dict`
        Dict of parameters for running analysis
    binary_params : `dict`
        Dict of parameters for binary params. For unmodeled pulsar these could
        be set to zero. If analysis_params['unmodeled']=True, this is done
        autmoatically.

    Returns
    -------

    """
 
    # if the print functions show only empty arrays, the new function can be used
    ptEst, sigma, f, df = load_pkl(analysis_params['h5file'])
    #ptEst, sigma, f, df = load_h5file(analysis_params['h5file'])
    # get segment start times
    segstarttimes =\
    np.arange(analysis_params['start_time'],analysis_params['start_time']+analysis_params['ndays']*86400, 7200)
    if analysis_params['unmodeled']:
        # hack up binary params to make sure they'll
        # give no binary doppler modulation
        binary_params = {}
        binary_params['P_orb'] = 1
        binary_params['asini'] = 0
        binary_params['T_asc'] = 0
    # get doppler shifts
    segs = segment_list_shifts_binary(segstarttimes, binary_params,
            analysis_params)
    dur = segstarttimes[-1] - segstarttimes[0]
    # get lowest and highest frequencies
    # for our two different cases
    if analysis_params['unmodeled']:
        # fix the number of bins you want
        nbins = analysis_params['nbins']
        highest_frequencies = f + nbins * df + df / 2.
        # lowest possible frequency (we want 4 bins below)
        lowest_frequencies = f - nbins * df - df / 2.
    else:
        total_spindown = dur * analysis_params['max_fdot']
        # highest possible frequency for each bin
        highest_frequencies = np.max(segs.get_shifts()) * (f+df/2.)
        # lowest possible frequency
        lowest_frequencies = np.min(segs.get_shifts()) * (f-df/2.) + total_spindown
    # combine frequencies
    # combine bins
    ptEst_new, sigma_new, number_combined, nc_above, nc_below, ratios =\
            combine_bins(ptEst, sigma, f, lowest_frequencies,
                    highest_frequencies, df)
    # do spinwandering calculation
    fbin_low = f - (nc_below * df) - (df / 2.)
    fbin_high = f + (nc_above * df) + (df / 2.)
    spinwander_f_total_below = lowest_frequencies - fbin_low
    spinwander_f_total_above = fbin_high - highest_frequencies
    # get observation time
    T_observation = segstarttimes[-1] - segstarttimes[0]
    # if ratios are bigger than cutoff then set them to zero
    # basically to make sure we have at least a certain fraction
    # of the bin present if we're going to set an upper limit
    ratios[np.where(ratios >= analysis_params['ratios_cutoff'])] = 0
    ratio_cut = np.where(ratios >= analysis_params['ratios_cutoff'])
    # set zeros where we fail the ratio cut
    ptEst_new[ratio_cut] = 0
    sigma_new[ratio_cut] = 0
    snr_data_nonan = ptEst_new / sigma_new
    snr_data_nonan[np.isnan(snr_data_nonan)] = 0
    # get maximum snr
    data_maximum_snr = np.max(snr_data_nonan)
    f_max_snr = f[np.where(snr_data_nonan == data_maximum_snr)]

    chunksize = np.ceil(analysis_params['n_sims'] /
            float(analysis_params['n_proc']))

    # speed things up by multithreading these guys
    print('Done combining bins...')
    print('Max SNR is %4.6f at %4.6f' % (data_maximum_snr, f_max_snr))
    print('Getting FAP...')
    processlist = []
    if analysis_params['n_proc'] > 50:
        raise ValueError('nproc should not be larger than 50')
    queue = ProcessQueue(analysis_params['n_proc'])
    # random number generator seeds
    for i in range(analysis_params['n_proc']):
        # set up the process we want to run
        # which is "run_simulation"
        # we pass it the number of simulations
        # and the info needed to combine gaussian bins
        process = Process(target=run_simulation,
                args=(queue, sigma, f, highest_frequencies,
                    lowest_frequencies, chunksize))
        # set up process list
        process.daemon=True
        processlist.append(process)
        # start the process
        process.start()
    maxSNRs = []
    for ii,process in enumerate(processlist):
        # get list of max SNRS for the simulations
        # we run on one of the threads
        maxSNRs.extend(queue.get())
        # block
        process.join()
    maxSNRs = np.asarray(maxSNRs)
    for process in processlist:
        process.terminate()
    queue.close()
    # pval spacing
    dp = 1. / maxSNRs.size
    # sort ascending
    maxSNRs = np.sort(maxSNRs)
    # ascending pvals
    pvals = np.arange(dp, 1+dp, dp)[::-1]
    pvals_interp = interp1d(maxSNRs, pvals)
    # get FAP of our
    try:
        max_snr_pval = pvals_interp(data_maximum_snr)
    except ValueError:
        print('max snr p-value is < %f' % np.min(pvals))
        max_snr_pval = np.min(pvals)

    print('Done getting FAP...')
    print('Getting upper limits...')
    # Correct for bin edges
    ptEst_corrected = np.copy(ptEst_new)
    sigma_corrected = np.copy(sigma_new)

    # Correct noise for notched frequencies
    sigma_corrected *= np.sqrt(ratios)
#    ptEst_corrected *= ratios
    ptEst_corrected[ptEst_corrected==0] = np.nan
    sigma_corrected[sigma_corrected==0] = np.nan

    # load marginalization correction
    marg_mat = loadmat(analysis_params['ul_ratio_file'])
    ul_ratio_snrs = marg_mat['snrs'].squeeze()
    ul_ratio = marg_mat['ul_ratio'].squeeze()
    ul_ratio_one_sigma = marg_mat['ul_ratio_one_sigma'].squeeze()
    conf_matfile = marg_mat['conf'].squeeze()
    # error checking may want to put this in `fix_params` at some point.
    if not analysis_params['conf'] == conf_matfile:
        raise ValueError('You are using the wrong files for scaling upper\
        limits! You have requested %f confidence upper limits and your\
        scaling is for %f confidence' % (analysis_params['conf'], conf_matfile))
    # get upper limits from marginalization
    uls = getUpperLimit_marg(ptEst_corrected, sigma_corrected, analysis_params['conf'],
        analysis_params['calibration_error'], ul_ratio, ul_ratio_snrs)
    # get circular upper limits
    ul_circ = getUpperLimit(ptEst_corrected, sigma_corrected,
        analysis_params['conf'], analysis_params['calibration_error'])

    # 1 sigma on h0 assuming no detection
    onesig = getUpperLimit_marg(np.zeros(sigma_corrected.size), sigma_corrected, 0.68,
        analysis_params['calibration_error'], ul_ratio_one_sigma, ul_ratio_snrs)
#    ul_circ = getUpperLimit_marg(np.zeros(sigma_corrected.size), sigma_corrected, 0.68,
#        analysis_params['calibration_error'])
    ########## END OF ANALYSIS!!!! #############
    # print(some stuff to screen about results and what's in matfile)

    print('***********************************')
    print('*********** RESULTS ***************')
    print('***********************************')
    print('The maxmimum SNR is %4.4f' % data_maximum_snr)
    print('The pvalue of this SNR is: %1.6f' % max_snr_pval)
    print('Saving things in final_results.mat')
    print('***********************************')
    print('***********************************')
    # put together final matfile
    results_mat = {
        'calibration_uncertainty' : analysis_params['calibration_error'],
        'conf' : analysis_params['conf'],
        'f' : f,
        'pte' : ptEst_new,
        'sigma' : sigma_new,
        'ul' : uls,
        'ul_circ' : ul_circ,
        'maxSNR' : data_maximum_snr,
        'max_snr_pval' : max_snr_pval,
        'f_max_snr' : f_max_snr,
        'snrs': snr_data_nonan,
        'simMAXSNRs' : maxSNRs,
        'pvals' : pvals,
        'ratio_for_rescaling_bins' : ratios,
        'sigma_corrected' : sigma_corrected,
        'pte_corrected' : ptEst_corrected,
        'one_sigma': onesig,
        'number_combined' : number_combined,
        'segstarttimes' : segstarttimes,
        'doppler_shifts' : segs.get_shifts(),
        'fdot' : analysis_params['max_fdot'],
        'raHr' : analysis_params['raHr'],
        'decDeg' : analysis_params['decDeg'],
        'spinwander_f_total_below': spinwander_f_total_below,
        'spinwander_f_total_above': spinwander_f_total_above,
        'T_observation':T_observation,
        'unmodeled': analysis_params['unmodeled']}

    savemat('%s/final_results.mat'%analysis_params['outputdir'],results_mat)








#!/usr/bin/env python
#
# Copyright (C) 2012 Santosh Sundaresan
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see
# <http://www.gnu.org/licenses/>.