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
from scipy.ndimage.interpolation import shift
import itertools
#import warnings
#warnings.filterwarnings("always")

def make_segmap(imgdata):

  mask = make_source_mask(imgdata, snr = 3, npixels = 5, dilate_size = 11)
  mean, median, rms = sigma_clipped_stats(imgdata, sigma = 3, mask = mask)

  threshold = photutils.detect_threshold(imgdata, snr = 3,
  background = median, error = rms)

  segmap = photutils.detect_sources(imgdata, threshold, npixels = 5)

  rms_array = np.ones(imgdata.shape)*rms

  return segmap, median, rms_array
  
def return_bkg(imgdata, segmap):

  mask = segmap
  mean, median, rms = sigma_clipped_stats(imgdata, sigma = 3, mask = mask)
  
  rms_array = np.ones(imgdata.shape)*rms
  
  return median, rms_array

def measure_morphology(imgdata, segmap, median, rmsarr, extent):

  imgdata -= median
  
  source_of_interest = segmap[(segmap.shape[0]-1)/2, (segmap.shape[1]-1)/2]

  segmap = photutils.SegmentationImage(segmap)
  segmap.keep_labels(labels = source_of_interest)

  if source_of_interest is not 0:
    try:
      source_morph = statmorph.source_morphology(imgdata, segmap,
      weightmap = rmsarr, cutout_extent = extent)
    except AssertionError:
      source_morph = None
  elif source_of_interest is 0:
    source_morph = None

  return source_morph

def mod_segmap_1pix(segmap):

  source_of_interest = segmap[(segmap.shape[0]-1)/2, (segmap.shape[1]-1)/2]

  segmap = photutils.SegmentationImage(segmap)
  segmap.keep_labels(labels = source_of_interest)

  segmap = segmap.data.astype(bool)

  print(segmap[90:110, 90:110].astype(int))
  
  segmap_expanded = None

  for shift_arr in itertools.product([-1, 0, 1], repeat = 2):
    if segmap_expanded is not None:
      segmap_expanded = segmap_expanded + shift(segmap, shift_arr)
      segmap_shrunk = segmap_shrunk & shift(segmap, shift_arr)
    else:
      segmap_expanded = segmap + shift(segmap, shift_arr)
      segmap_shrunk = segmap & shift(segmap, shift_arr)
      
  segmap_expanded = segmap_expanded.astype(int)
  segmap_shrunk = segmap_shrunk.astype(int)
    
  print(segmap_expanded[90:110, 90:110])
  print(segmap_shrunk[90:110, 90:110])
    
  return segmap_expanded, segmap_shrunk
  
def write_gini_masks(imgdata, header, all_morphs, filename):

  writemask = np.zeros(imgdata.shape)

  for morph in all_morphs:

    writemask[morph._slice_stamp][morph._segmap_gini] = 1.

  hdu = fits.PrimaryHDU(header = header, data = writemask)
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
    #segm, bkg_median, bkg_rms_arr = make_segmap(img[0].data)

    #hdu = fits.PrimaryHDU(header = img[0].header, data = segm.data)
    #hdu.writeto(fileroot+'seg.fits', overwrite = True)

    segm_fname = fileroot+'seg.fits'

    segm = (fits.open(segm_fname))[0].data
    
    segm += 1
    
    bkg_median, bkg_rms_arr = return_bkg(img[0].data, segm.astype(bool))

    obj_morph = measure_morphology(img[0].data, segm, bkg_median,
    bkg_rms_arr, 2)
    
    if obj_morph:

      mod_segmap_1pix(segm)
      
      write_gini_masks(img[0].data, img[0].header, obj_morph,
      fileroot+'gini.fits')

      t = build_table(obj_morph, img[0].header)

      t['filename'] = arg

      t.write(fileroot+'cat', format = 'ascii.fixed_width', overwrite = True)
      table = vstack([table, t])

  if table:
    table.write('allmorphs.cat', format = 'ascii.fixed_width', overwrite = True)

if __name__ == '__main__':
  main()
