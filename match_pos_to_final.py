# Commands to select directory based on hostname
from socket import gethostname

path_start = '/Users/itcadmin/documents/highz_elais/make_cutouts-master/'

##################################################

import numpy as np
from matplotlib import pyplot as plt
import glob

# For catalog matching
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.coordinates import match_coordinates_sky
from astropy.coordinates import search_around_sky

from astropy.table import Table
##################################################

"""
Script to match a list of positions to the radio catalogue
	- And, outptus a workflow.txt file with Source_Name for visual checks
	- Also, creates a subset of radio catalogue for these sources
"""

# Read in radio catalogue
final_fname = glob.glob("../data/final_cross_match_catalogue-v0.3.fits")[-1]
mlfin_srl = Table.read(final_fname, character_as_bytes=False)

# Read in file with a list of positions
tab_to_match = Table.read("workflow.txt", format='ascii')

radio_coords = SkyCoord(mlfin_srl["RA"], mlfin_srl["DEC"], unit='deg', frame='icrs')
match_coords = SkyCoord(tab_to_match["RA"], tab_to_match["DEC"], unit='deg', frame='icrs')

# ind_r, ind_m, sep2d, _ = search_around_sky(radio_coords, match_coords, seplimit=5*u.arcsec)
ind_r, sep2d, _ = match_coordinates_sky(match_coords, radio_coords, nthneighbor=1)

mlfin_srl["Notes"] = -1
mlfin_srl["Sep"] = np.nan
mlfin_srl["Sep"][ind_r] = sep2d.arcsec
mlfin_srl["Source_Name", "RA", "DEC", "Sep", "z1_median", "Notes"][ind_r].write("workflow.txt", format='ascii.commented_header')
