# -*- coding: utf-8 -*-
"""
@author: Michael Yeh
"""


from __future__ import print_function
import time
import argparse
import multiprocessing
import numpy as np
import scipy.io as sio
import pyfftw


def SEC_C(data, template, k, moveout, weight):
    """ Super Efficient Cross-correlation Python code for seismic data

    SEC_C is designed to handle multiple stations, multiple components and
    multiple templates matched filtering of time series data specifically for
    seismic applications in an efficient time with an efficient memory usage.
    We have adopted the input parameter style as Beauce et al., 2017 matched
    filter code (https://github.com/beridel/fast_matched_filter). For testing
    SEC_C user also can use the Beauce et al., 2017 test code.

    Normalization part inspired by Mass algorithm
    (http://www.cs.unm.edu/~mueen/FastestSimilaritySearch.html)

    Required libraries: numpy, scipy, and pyfftw
    Tested with Python 3.5.2, numpy 1.13.3, scipy 0.18.1, and pyfftw 0.10.4

    MATLAB implemented by: Nader Shakibay Senobari, Spring 2018
    Python implemented by: Chin-Chia Michael Yeh, Summer 2018

    The Python version is implemented based on the MATLAB version

    Citation: Shakibay Senobari et al., 2018, submitted to SRL.

    Parameters
    ----------
    data : numpy array, shape (len_data, n_component, n_station)
        input data
    template : numpy array, shape (len_template, n_component, n_station,
                                   n_template)
        input tempate
    k : int
        A tunning parameter, we recommend assigning a power of two for k (e.g.
        2^12 for one day of 20 Hz data or 2^13 for one day of 50-100 Hz data)
        for efficient performance. However, in all cases we advise running some
        test cases to find values of k to optimize the run time. For more
        information we refer you to the citation above.
    moveout : numpy array, shape (n_station, n_template)
        input moveouts
    weight : numpy array, shape (n_station, n_template)
        input weights

    Returns
    -------
    ccc_sum : numpy array, shape (len_data - len_template + 1, n_template)
        output ccc_sum

    Notes
    -----
    moveouts and weights can be extended to include components (i.e. become a
    3D array). In this case user can easily change line 178 to account for
    that.

    User also can save individual CCC for each channel-template, see below:
    """

    # let's gather some information about data and templates such as the length
    # of templates, number of components, number of stations, number of
    # tempplates and length of data.
    data = np.array(data, dtype=float, copy=True)
    template = np.array(template, dtype=float, copy=True, order='F')
    moveout = np.array(moveout, dtype=int, copy=True)
    weight = np.array(weight, dtype=float, copy=True)
    len_template, n_component, n_station, n_template = template.shape
    len_data = data.shape[0]
    k = int(k)

    # initialize fftw
    n_core = multiprocessing.cpu_count()
    buffer_shape = _get_buffer_size(len_data, k, len_template - 1)
    query_shape = [k, n_template]
    buffer_freq_shape = buffer_shape[:]
    buffer_freq_shape[0] = buffer_freq_shape[0] // 2 + 1
    query_freq_shape = query_shape[:]
    query_freq_shape[0] = query_freq_shape[0] // 2 + 1

    data_segment_fftw = pyfftw.empty_aligned(
        buffer_shape, dtype='float64', order='F')
    data_segment_freq_fftw = pyfftw.empty_aligned(
        buffer_freq_shape, dtype='complex128', order='F')
    query_fftw = pyfftw.empty_aligned(
        query_shape, dtype='float64', order='F')
    query_freq_fftw = pyfftw.empty_aligned(
        query_freq_shape, dtype='complex128', order='F')
    dot_prod_fftw = pyfftw.empty_aligned(
        buffer_shape, dtype='float64', order='F')
    dot_prod_freq_fftw = pyfftw.empty_aligned(
        buffer_freq_shape, dtype='complex128', order='F')
    data_segment_fft = pyfftw.FFTW(
        data_segment_fftw, data_segment_freq_fftw,
        axes=(0, ), threads=n_core)
    query_fft = pyfftw.FFTW(
        query_fftw, query_freq_fftw,
        axes=(0, ), threads=n_core)
    dot_prod_ifft = pyfftw.FFTW(
        dot_prod_freq_fftw, dot_prod_fftw,
        axes=(0, ), threads=n_core,
        direction='FFTW_BACKWARD')

    # initialize ccha and ccc_sum
    ccha = np.zeros((buffer_shape[0] - len_template + 1, buffer_shape[1]),
                    order='F')
    ccc_sum = np.zeros((len_data - len_template + 1, n_template))

    # we want to find a station with largest moveout (i.e. the last station
    # that detect the signal) and later on we pad zeros to the other stations
    # at the begining of the data to align them.
    moveout = np.max(moveout, axis=0) - moveout

    # uncomment below if you want to scale ccc_sum to one.
    # weight /= n_component

    # loop over stations
    for i in range(n_station):
        # loop over componets
        for j in range(n_component):
            # get the template and data for each channel
            query = template[:, j, i, :]
            # do some preprocessing for normalization, we need this later on
            query_rms = np.sqrt(np.sum(query ** 2, axis=0))

            # buffer is a function that divides the data into pieces with the
            # length of k and with overlaps of len_template - 1 samples and
            # make a 2D array from the data vector.
            data_segment = _buffer(data[:, j, i], k, len_template - 1)

            # another preprocessing step for normalization
            data_rms = np.cumsum(data_segment ** 2, axis=0)
            data_rms[len_template:, :] = (
                data_rms[len_template:, :] -
                data_rms[:-(len_template), :])
            data_rms = np.sqrt(data_rms)

            # calculating cross-correlation (CC) in frequency domain
            # reverse the templates
            query = np.flipud(query)
            # padding with zeros
            query = np.concatenate((query,
                                    np.zeros((k - len_template,
                                              n_template),
                                             order='F')))
            # transfering to the frequency domain for the
            data_segment_fftw[:] = data_segment
            data_freq = data_segment_fft()
            # transfering to the frequency domain for templates
            query_fftw[:] = query
            query_freq = query_fft()

            # loop over templates
            for l in range(n_template):
                # do the dot product
                dot_prod_freq_fftw[:] = (data_freq *
                                         np.tile(query_freq[:, l],
                                                 (data_freq.shape[1], 1)).T)
                # going back to the time domain
                dot_prod = dot_prod_ifft()

                # devide by the normalization factor
                ccha[:] = (dot_prod[len_template - 1:, :] /
                           (data_rms[len_template - 1:, :] * query_rms[l]))

                # calculate the CCC and sum over stations and components
                ccc_sum[:, l] = (
                    weight[i, l] *
                    (np.concatenate((
                        np.zeros(moveout[i, l]),
                        ccha.flatten('F')[
                            len_template - 1:len_data - moveout[i, l]]))) +
                    ccc_sum[:, l])

                # if you are interested in having individual CCCs for each
                # channel and templates comment out the above statement and
                # uncomment below; you should give a directory
                #
                # ("your_directory_") for saving the CCCs.
                # import os
                # dir_path = os.path.join('.', 'yourdirectory_')
                # if not os.path.isdir(dir_path):
                #     os.makedirs(dir_path)
                # out_fname = 'CCC_{0:d}_{1:d}_{2:d}'.format(i, j, l)
                # out_path = os.path.join(dir_path, out_fname)
                # sio.savemat(
                #     out_path,
                #     {'ccha': ccha.flatten('F')[len_template - 1:len_data]})
    return ccc_sum


def _get_buffer_size(len_data, len_segment, len_overlap):
    """ Get the size of buffer

    Parameters
    ----------
    len_data : int
        length of data
    len_segment : int
        length of segment
    len_overlap : int
        length of overlap

    Returns
    -------
    buffer_shape : list, shape (2)
        output shape of buffer
    """
    len_data += len_overlap
    n_segment = np.ceil(
        (len_data - len_overlap) / (len_segment - len_overlap))
    n_segment = int(n_segment)
    return [len_segment, n_segment]


def _buffer(data, len_segment, len_overlap):
    """ Segment the data into chunks with sliding window of size len_segment
    and overlap of size len_overlap

    Parameters
    ----------
    data : numpy array, shape (len_data)
        input data
    len_segment : int
        length of segment
    len_overlap : int
        length of overlap

    Returns
    -------
    data_segment : numpy array, shape (len_segment, n_segment)
        output segmented data
    """

    data = np.concatenate((np.zeros(len_overlap), data))
    n_segment = np.ceil(
        (data.shape[0] - len_overlap) / (len_segment - len_overlap))
    n_segment = int(n_segment)
    data_segment = np.zeros((len_segment, n_segment), order='F')
    segment_st = 0
    for i in range(n_segment):
        segment = data[segment_st:min(segment_st + len_segment,
                                      data.shape[0])]
        data_segment[:segment.shape[0], i] = segment
        segment_st += (len_segment - len_overlap)
    return data_segment


def main():
    args = parser.parse_args()
    input_mat = args.input_mat
    output_mat = args.output_mat

    parameter = sio.loadmat(input_mat)
    t_tic = time.time()
    ccc_sum = SEC_C(parameter['data'],
                    parameter['template'],
                    parameter['k'][0][0],
                    parameter['moveout'],
                    parameter['weight'])
    t_toc = time.time() - t_tic
    print('Elapsed time is {0:.6f} seconds.\n'.format(t_toc))
    sio.savemat(output_mat, {'ccc_sum': ccc_sum, 't_toc': t_toc})


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-i', '--input_mat', help='path to the input mat file')
    parser.add_argument(
        '-o', '--output_mat', help='path to the output mat file')
    main()
