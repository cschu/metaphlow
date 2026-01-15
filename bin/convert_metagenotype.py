#!/usr/bin/env python
# coding: utf-8


import numpy as np
import pandas as pd
import xarray as xr
import sfacts as sf
import os
import sys
import gc
import gzip
import tarfile
import tempfile
import lzma
import argparse


def read_params():
    """
    Parse command line arguments
    """
    parser = argparse.ArgumentParser(
        description='Convert SameStr npz into StrainFacts Metagenotype netCDF',
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
        help='Path to SameStr array (npy/npz)'
    )
    
    parser.add_argument(
        '-o', '--outpath',
        type=str,
        required=True,
        help='Output path of StrainFacts Metagenotype (netCDF)'
    )

    parser.add_argument(
        '-s', '--study_stats',
        type=str,
        required=True,
        help='Path to within-population (within-study) statistics (tarball/netCDF)'
    )
        
    parser.add_argument(
        '-m', '--member',
        help="Path to the NetCDF file inside the tarball (required if input is a tarball)",
    )

    parser.add_argument(
        '-p', '--pos_vec',
        type=str,
        default=None,
        help='Fixed positions table (optional, TSV)'
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


def load_statistics(fpath, member=None):
    if member is None:
        pos_stats = xr.open_dataset(fpath)
    else:
        with tarfile.open(fpath, "r:*") as tar:
            with tar.extractfile(member) as src:
                with tempfile.NamedTemporaryFile(suffix=".nc") as tmp:
                    tmp.write(src.read())
                    tmp.flush()

                    pos_stats = xr.open_dataset(tmp.name)
    return pos_stats


def aggregate_nt2allele(xr_ds):
    # Mask nucleotide coverage based on allele category, then sum per category
    result_dict = {}
    for allele_type in ['ref', 'alt']:
        mask = xr_ds.coords['allele'] == allele_type
        result_dict[allele_type] = xr_ds.where(mask, 0).sum(dim='nucleotide')

    # Reconstruct xarray Dataarray
    result = xr.concat([result_dict['alt'], result_dict['ref']],
                       dim=xr.DataArray(['alt', 'ref'], dims='allele', name='allele'))
    result = result.transpose('sample', 'position', 'allele')

    return result


def main(inpath, outpath, pos_stats, positions=None):
    # Load SameStr 3D array and convert to xarray
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


    # Remove intermediates to free memory
    del np_arrs, coords
    gc.collect()


    # Load study statistics and extract major alleles (based on max allele frequency)
    study_major = pos_stats.mean_freq.idxmax('nucleotide')


    # Construct xarray Dataset (also subset to global position if provided)
    if positions is not None:
        # fixed_pos[fixed_pos['position'].isin([70304])]
        fixed_pos = pd.read_csv(positions, sep='\t', index_col=None)
        fixed_pos = fixed_pos[fixed_pos["position"].isin(raw_xr["position"].to_numpy())]

        mgen_ds = (
            raw_xr
            .to_dataset(name='metagenotype')
            .sel(position=fixed_pos['position'].values)
            .assign(major_nt=('position', fixed_pos['nucleotide'].values))
        )
    else:
        mgen_ds = (
            raw_xr
            .to_dataset(name='metagenotype')
            .assign(major_nt = study_major)
        )
    mgen_ds['metagenotype'] = mgen_ds.metagenotype.astype(int)


    # Construct allele group
    allele_grp = (
        xr.where(mgen_ds.nucleotide == mgen_ds.major_nt, 'ref', 'alt')
        .transpose('position', 'nucleotide')
    )
    mgen_ds = mgen_ds.assign_coords(allele=allele_grp)

    # Convert nucleotide to allele
    mgen_converted = aggregate_nt2allele(mgen_ds)

    # Sort by positions
    mgen_converted = mgen_converted.sortby('position')

    # Convert xarray to sfacts.data.Metagenotype
    sf_mgen = sf.data.Metagenotype(mgen_converted.metagenotype)

    # Export converted Metagenotype to netCDF file
    sf_mgen.dump(outpath)


if __name__ == "__main__":
    args = read_params()

    # Check format of study statistics
    if tarfile.is_tarfile(args.study_stats):
        if not args.member:
            parser.error("--member is required when --input is a tarball")
        else:
            pos_stats = load_statistics(args.study_stats, args.member)
    else:
        if args.member:
            parser.error("--member should not be used when --input is a NetCDF file")
        else:
            pos_stats = load_statistics(args.study_stats)

    main(args.inpath, args.outpath, pos_stats, args.pos_vec)
