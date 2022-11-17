from concurrent.futures import ProcessPoolExecutor
import numpy as np
from dateutil.parser import parse
from datetime import *; from dateutil.relativedelta import relativedelta

from pathlib import Path

import ccdproc
from astropy.nddata import CCDData
from astropy.io import fits


def average_comb(imgs, saveAs):
    toCombine = [CCDData(img, unit='adu') for img in imgs]
    med = ccdproc.combine(toCombine,
                            method='average',
                            sigma_clip=True, sigma_clip_low_thresh=3, sigma_clip_high_thresh=3,
                            sigma_clip_func=np.ma.average)
    med.meta['combined'] = True
    med.write(saveAs, overwrite=True)
    print('saved to: ' + saveAs)
    return med

def increment_date(strdate, tincrement):
    parsed = parse(strdate)
    later = parsed + relativedelta(seconds=tincrement)
    incremented_dateobj = later.strftime('%Y-%m-%dT%X.0000')
    return incremented_dateobj

def marckie_slice_combine_reduce(filepath, output_dir, chunk_size, exp_time, darkbias, flat, do_offset=False):

    extra_offset = int(filepath[-6]) * 2000 * do_offset # because the cubes are split into 2

    filepath = Path(filepath)
    print(f'working on file{filepath}')
    with fits.open(filepath) as f:
        data = f[0].data
        header = f[0].header
        start_time = header['DATE-OBS']


    for i in range(0, data.shape[0], chunk_size):
        images = data[i:i+chunk_size, ...]
        combined = images.mean(axis=0)
        new_header = header.copy()
        new_header['DATE-OBS'] = date_n_time = increment_date(start_time, exp_time * i + extra_offset)
        
        new_name = f'combined_{date_n_time}__{filepath.name}'
        print(new_name)
        # reduce it
        print(f'reducing image {i} from file {filepath.name} with timestamp offset {exp_time * i + extra_offset}')

        corrected = (combined - darkbias) / flat

        out_file = Path(output_dir).joinpath(new_name)
        
        fits.writeto(
            out_file,
            corrected,
            header=new_header,
            overwrite=True
        )
        print(f'Successfully sliced and combined {new_name}')

def fast_process(image_cube_list, output_dir, chunk_size, exp_time, darkbias, flat, do_offset=False):
    '''Combine images without slicing'''
    length = len(image_cube_list)
    with ProcessPoolExecutor() as pool:
        print('starting process')
        pool.map(marckie_slice_combine_reduce, image_cube_list, [output_dir]*length, [chunk_size]*length, [exp_time]*length, [darkbias]*length, [flat]*length, [do_offset]*length)