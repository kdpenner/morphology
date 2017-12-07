import numpy as np
import photutils
import statmorph
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from photutils import make_source_mask
from astropy.table import Table
from astropy.wcs import WCS
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

  return segmap

def measure_morphology(imgdata, segmap):

  imgdata -= median
  source_morphs = statmorph.source_morphology(imgdata, segmap,
  border_size = 0, cutout_extent = 2)

  return source_morphs
  
def write_gini_mask(imgdata, one_morph, filename)

  writemask = np.zeros(imgdata.shape)
  writemask[one_morph._slice_stamp][one_morph._segmap_gini] = 1.

  hdu = fits.PrimaryHDU(header = img[0].header, data = writemask)
  hdu.writeto(filename, overwrite = True)
  
def build_table(all_morphs, header):

  rowdata = np.zeros(11, len(all_morphs))

  for i, morph in xenumerate(all_morphs):

    wcs = WCS(header)

    posasym = wcs.all_pix2world(morph.xc_asymmetry, 
                                morph.yc_asymmetry, 1, 
                                ra_dec_order = True)


    rowdata[:, i] = [morph.asymmetry, posasym[0], 
             posasym[1], morph.concentration, 
             morph.smoothness, morph.rhalf_ellip,
             morph.rpetro_ellip, morph.gini, 
             morph.m20, morph.flag,
             morph.sn_per_pixel]

  names = ('asymmetry', 'RA center for asymmetry', 'Dec center for asymmetry',
           'concentration', 'smoothness', 'half light elliptical semimajor axis length',
           'petrosian elliptical semimajor axis length', 'gini', 'm20', 
           'flag', 'S/N per pixel')

  t = Table(rows = rowdata, names = names)

  return t
