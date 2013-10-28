
# 
#  This file is part of Herschel Common Science System (HCSS).
#  Copyright 2001-2013 Herschel Science Ground Segment Consortium
# 
#  HCSS is free software: you can redistribute it and/or modify
#  it under the terms of the GNU Lesser General Public License as
#  published by the Free Software Foundation, either version 3 of
#  the License, or (at your option) any later version.
# 
#  HCSS is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#  GNU Lesser General Public License for more details.
# 
#  You should have received a copy of the GNU Lesser General
#  Public License along with HCSS.
#  If not, see <http://www.gnu.org/licenses/>.
# 
###########################################################################
###                    SPIRE HiRes Mapper Script                        ###
###########################################################################
# Title: Photometer_HiResMapping.py
#
# Purpose:
# This script will show how to set up beam profiles and to run the
# HiRes superresolution mapper on SPIRE Level 1 data
#
# Usage:
# 1) Load multiple observations (require obs ID and Pool name for each observation),
# 2) Collect the the level 1 products (scan lines) together from each observation, 
# 3) Perform baseline subtraction on level 1 timelines using a median subtraction,  
#    (for alternative baseline removal algorithms see the SPIRE data Redction Guide),
# 4) Perform Map Making on scan lines from all observations  for each SPIRE array.
#
# Required Inputs:
# a) Array of Observation IDs
# b  band to process, e.g. 'PSW', 'PMW' or 'PLW'
# c) working directory containing the beam profiles
#
# Assumptions: Data is obtained from the HSA and has been already destriped.
# 
# Last Updated: 25 September 2013
#
#
#   Written for HIPE v12.0
#    
###########################################################################
#
###########################################################################
###                     User Selectable Options                         ###
###########################################################################
# (A) List of OBSIDs in the form of an integer or hexadecimal (0x) number:
#     (Example data for M33) 
# (B) image size in arcminutes
# (C) image central coordinates in decimal degrees
# (D) band to process, e.g. 'PSW', 'PMW' or 'PLW'
#
# Maybe do the Tadpole? Or stick with M33?
# Need to specify: the array
# the beam (show how to retrieve it)
# the beam directory
# the obsids
obsids  = [1342189079, 1342189080]
imagesize = [65, 65] # y-, x- dimensions in arcminutes
imagecenter = [23.466769672650443, 30.663356134471154] # RA, Dec in decimal degrees
band = 'PLW'  # PSW, PLW, or PMW
###########################################################################
# hiresMapper demo 

# Prepare 1" beam - start with http://herschel.esac.esa.int/twiki/bin/view/Public/SpirePhotometerBeamProfile
# downloaded to the var.hcss.workdir directory
import urllib, os
workDir = Configuration.getProperty('var.hcss.workdir')
beamName = "0x5000241aL_%s_pmcorr_1arcsec_norm_beam.fits"%band
urllib.urlretrieve ("https://nhscsci.ipac.caltech.edu/spire/data/beam_profiles/"+beamName,\
    os.path.join(workDir,beamName))
bcenter = {'PSW':(700,699), 'PMW':(700,700), 'PLW':(698,700)} # Positions of peak pixel
beamfull = fitsReader(file = os.path.join(workDir,beamName))
beam = crop(beamfull, int(bcenter[band][0] - 100) , int(bcenter[band][1] - 100),\
            int(bcenter[band][0] + 101), int(bcenter[band][1] + 101))

assert beam.dimensions[0] % 2 == 1, "Beam dimensions are not odd"
assert beam.dimensions[1] % 2 == 1, "Beam dimensions are not odd"
assert beam.getIntensity(beam.dimensions[0] / 2, beam.dimensions[1] / 2) \
    == MAX(NAN_FILTER(beam.image)), "Beam is not centred on central pixel"

level1Corrected = Level1Context()
# Retrieve timeline data
for obsid in obsids:
    obs = getObservation(obsid, useHsa=True, instrument='SPIRE')
    for ref in obs.level1.refs:
        level1Corrected.addRef(ref)

level1Corrected = destriper(level1Corrected, array=band, useSink=True)[0]

# Prepare Wcs with half the pixel size of standard map
wcs = obs.level2.getProduct("psrc"+band).wcs.copy()
wcs.crval1 = imagecenter[0]
wcs.crval2 = imagecenter[1]
wcs.cdelt1 /= 2.0
wcs.cdelt2 /= 2.0
wcs.naxis1 = int(imagesize[1]/60./abs(wcs.cdelt1) + 0.5)
wcs.naxis2 = int(imagesize[0]/60./abs(wcs.cdelt2) + 0.5)
wcs.crpix1 = (wcs.naxis1 + 1) / 2.
wcs.crpix2 = (wcs.naxis2 + 1) / 2.

# FOUR EXAMPLES FOLLOW:

# (1) 20 iterations (default for maxIter parameter), output final iteration only
hiresImage, hiresBeam = hiresMapper(level1Corrected, beam=beam, wcs=wcs)

# (2) Output first iteration only
hiresImageIter1, hiresBeamIter1 = hiresMapper(level1Corrected, beam=beam, wcs=wcs, maxIter=1)

# (3) 20 iterations, output multiple iterations
hiresImages, hiresBeams = hiresMapper(level1Corrected, beam=beam, wcs=wcs, \
    storeIter=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,20])
hiresImage1 = hiresImages[1]
hiresImage20 = hiresImages[20]
hiresBeam1 = hiresBeams[1]
hiresBeam20 = hiresBeams[20]

# (4) Stop after N-iterations, adjust the image, then continue
n = 2
hiresImageIterN, hiresBeamIterN = hiresMapper(level1Corrected, beam=beam, wcs=wcs, maxIter=n)
# Do something to hiresImageIterN
# ...
# Continue processing
hiresImageFinal, hiresBeamFinal = hiresMapper(level1Corrected, beam=beam, \
    startImage=hiresImageIterN)

