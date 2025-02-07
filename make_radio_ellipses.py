import numpy as np
from matplotlib import pyplot as plt
import glob

# For catalog matching
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.coordinates import match_coordinates_sky
from astropy.coordinates import search_around_sky

from astropy.table import Table
import os
##################################################

"""
Script to make radio source ellipses from the final catalogue based on a list of positions
"""


def make_ds9_reg(ra_positions, dec_positions, major_rad, minor_rad, pa, output_regfile, marker_colours, lgz_size, lgz_width, lgz_pa):
    first_regline = "# Region file format: DS9 version 4.1"
    global_reg_def = 'global color={0} dashlist=8 3 width=2 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1'.format(marker_colours)

    with open(output_regfile, "w") as fout:
        fout.write(first_regline + "\n")
        fout.write(global_reg_def + "\n")
        fout.write("fk5\n")

        for ii, ra in enumerate(ra_positions):
            if np.isnan(lgz_size[ii]) or np.isnan(lgz_width[ii]) or np.isnan(lgz_pa[ii]):
                fout.write('ellipse({0},{1},{2}",{3}",{4})\n'.format(ra, dec_positions[ii],
                                                                     major_rad[ii], minor_rad[ii], pa[ii]))
            else:
                fout.write('ellipse({0},{1},{2}",{3}",{4})\n'.format(ra, dec_positions[ii],
                                                                     lgz_size[ii], lgz_width[ii], lgz_pa[ii]))

    return


# RA, DEC, MAJ rad, MIN RAD, angle(check this!)
# t['RA'],t['DEC'],t['Maj']*2/overlay_scale,t['Min']*2/overlay_scale, angle=90+t['PA']

# Read in radio catalogue
final_fname = glob.glob("../data/final_cross_match_catalogue-v0.3.fits")[-1]
mlfin_srl = Table.read(final_fname, character_as_bytes=False)

# Workflow
cata = Table.read("workflow.txt", format='ascii')

wflow_coords = SkyCoord(cata["RA"], cata["DEC"], unit='deg', frame='icrs')
final_coords = SkyCoord(mlfin_srl["RA"], mlfin_srl["DEC"], unit='deg', frame='icrs')

srad = 300

OUTDIR_POS = "cata_pos/"
if not os.path.exists(OUTDIR_POS):
        os.makedirs(OUTDIR_POS)

for k in range(len(cata)):
    fin_fname = "{0}ellipse_{1}.reg".format(OUTDIR_POS, cata["Source_Name"][k])
    if not os.path.exists(fin_fname):
        print(fin_fname)
        cent_coord = SkyCoord(cata["RA"][k], cata["DEC"][k], unit='deg', frame='icrs')
        fin_near_cent = cent_coord.separation(final_coords).arcsec < srad
        print(np.sum(fin_near_cent))
        # Now write ellipses for all sources in fin_near_cent
        make_ds9_reg(mlfin_srl["RA"][fin_near_cent].tolist(), mlfin_srl["DEC"][fin_near_cent].tolist(), (mlfin_srl["Maj"][fin_near_cent]*3600).tolist(),
                     (mlfin_srl["Min"][fin_near_cent]*3600).tolist(), (mlfin_srl["PA"][fin_near_cent]+90).tolist(), fin_fname, "cyan",
                     mlfin_srl["LGZ_Size"][fin_near_cent].tolist(), mlfin_srl["LGZ_Width"][fin_near_cent].tolist(), (mlfin_srl["LGZ_PA"][fin_near_cent]+90).tolist())
