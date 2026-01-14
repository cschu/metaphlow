#!/usr/bin/env python
# coding: utf-8

import xarray as xr
import numpy as np
import os
import sys
import gzip
import argparse


def read_params():
    """
    Parse command line arguments
    """
    parser = argparse.ArgumentParser(
        description='Calculate per-position, per-nucleotide statistics of variants',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
            Examples:
            %(prog)s "$input" "$output"
        '''
    )
    
    parser.add_argument(
        '-i', '--inpath',
        type=str,
        required=True,
        help='Absolute path to numpy array containing coverage info'
    )
    
    parser.add_argument(
        '-o', '--outpath',
        type=str,
        required=True,
        help='Output path of calculated statistics (NetCDF format)'
    )

    return parser.parse_args()


def load_numpy_file(input_file):
    """
        Load numpy file
    """
    if input_file.endswith('.npz'):
        return np.load(input_file, allow_pickle=True)['arr_0']
    elif input_file.endswith('.npy.gz'):
        with gzip.open(input_file, 'rb') as f:
            return np.load(f, allow_pickle=True)
    else:
        return np.load(input_file, allow_pickle=True)


def main(inpath, outpath):
    # Load file and convert NaN to 0
    np_arrs = np.nan_to_num(load_numpy_file(inpath), nan=0)
    np_shapes = np_arrs.shape

    # Convert numpy array to xarray
    coords = {
        'sample':np.arange(0, np_shapes[0]).tolist(),
        'position':np.arange(0, np_shapes[1]).tolist(),
        'nucleotide':['A', 'C', 'G', 'T']
    }
    raw_xr = xr.DataArray(np_arrs,
                          dims=('sample', 'position', 'nucleotide'),
                          coords=coords)
    freq_xr = (raw_xr > 0).astype(int)


    # Basic statistics (within-study, per-position, per-nucleotide)
    mean_cov = raw_xr.mean(dim="sample")
    median_cov = raw_xr.median(dim="sample")
    sum_cov = raw_xr.sum(dim="sample")
    mean_freq = freq_xr.mean(dim="sample")
    median_freq = freq_xr.median(dim="sample")
    sum_freq = freq_xr.sum(dim="sample")
    nvar = (sum_freq > 0).astype(int).sum(dim='nucleotide')


    # Collect statistics DataArrays (keep only covered positions)
    summary_stats = xr.Dataset({
        "mean_coverage": mean_cov,
        "median_coverage": median_cov,
        "sum_coverage": sum_cov,
        "mean_freq": mean_freq,
        "median_freq": median_freq,
        "sum_freq": sum_freq,
        "nvar": nvar
    }).sel(position=nvar.position[nvar > 0])

    # Export statistics to NetCDF
    summary_stats.to_netcdf(outpath)


if __name__ == "__main__":
    args = read_params()

    main(args.inpath, args.outpath)
