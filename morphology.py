#!/usr/bin/env python -tt

import numpy as np
import photutils
import statmorph
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from photutils import make_source_mask
from astropy.table import Table
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
import sys
from astropy.io import ascii
import re
from astropy.table import vstack
import astropy.units as u
#import warnings
#warnings.filterwarnings("always")

def make_segmap(imgdata):

  mask = make_source_mask(imgdata, snr = 3, npixels = 5, dilate_size = 11)
  mean, median, rms = sigma_clipped_stats(imgdata, sigma = 3, mask = mask)

  threshold = photutils.detect_threshold(imgdata, snr = 3,
  background = median, error = rms)

  segmap = photutils.detect_sources(imgdata, threshold, npixels = 5)

#  hdu = fits.PrimaryHDU(header = img[0].header, data = segmap.data)
#  hdu.writeto('segtest.fits', overwrite = True)

  return segmap, median

def measure_morphology(imgdata, segmap, median):

  imgdata -= median
  source_morphs = statmorph.source_morphology(imgdata, segmap,
  border_size = 0, cutout_extent = 10)

  return source_morphs
  
def write_gini_mask(imgdata, one_morph, filename):

  writemask = np.zeros(imgdata.shape)
  writemask[one_morph._slice_stamp][one_morph._segmap_gini] = 1.

  hdu = fits.PrimaryHDU(header = img[0].header, data = writemask)
  hdu.writeto(filename, overwrite = True)
  
def build_table(all_morphs, header):

  columndata = np.zeros((len(all_morphs), 11))

  for i, morph in enumerate(all_morphs):

    wcs = WCS(header)

    posasym = wcs.all_pix2world(morph.xc_asymmetry, 
                                morph.yc_asymmetry, 0, 
                                ra_dec_order = True)
                                
    scale = proj_plane_pixel_scales(wcs.celestial)
    
    radiusconv = (scale[0]*u.deg).to('arcsec').value

    columndata[i, :] = [morph.asymmetry, posasym[0], 
             posasym[1], morph.concentration, 
             morph.smoothness, morph.rhalf_ellip*radiusconv,
             morph.rpetro_ellip*radiusconv, morph.gini, 
             morph.m20, morph.flag,
             morph.sn_per_pixel]

  names = ('asymmetry', 'RA center for asymmetry', 'Dec center for asymmetry',
           'concentration', 'smoothness', 'half light elliptical semimajor axis length',
           'petrosian elliptical semimajor axis length', 'gini', 'm20', 
           'flag', 'S/N per pixel')

  t = Table(columndata, names = names)
  t['RA center for asymmetry'].unit = u.deg
  t['Dec center for asymmetry'].unit = u.deg
  t['half light elliptical semimajor axis length'].unit = u.arcsec
  t['petrosian elliptical semimajor axis length'].unit = u.arcsec

  return t

def main():

  args = sys.argv[1:]
  
  table = Table()
  
  for arg in args:
    img = fits.open(arg)
    fileroot = re.split('fits\Z', arg)[0]
    segm, bkg_median = make_segmap(img[0].data)
    obj_morphs = measure_morphology(img[0].data, segm, bkg_median)
    t = build_table(obj_morphs, img[0].header)
    t['filename'] = arg
    t.write(fileroot+'cat', format = 'ascii.fixed_width', overwrite = True)
    table = vstack([table, t])

  table.write('allmorphs.cat', format = 'ascii.fixed_width', overwrite = True)

if __name__ == '__main__':
  main()
