# 
#  This file is part of Herschel Common Science System (HCSS).
#  Copyright 2001-2011 Herschel Science Ground Segment Consortium
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
###                     SPIRE Astrometry Fix Script                     ###
###########################################################################
# Title: Photometer_AstrometryCorrection.py
#
# Purpose:
# This script willto fix astrometry in SPIRE observations 
# Large map, Small Map or Parallel mode can be corrected. 
# BOTH Level 1 AND Level 2 products are updated in the Observation Context.
#
# Usage:
# Select the appropriate section of the script, and RUN the SELECTION
# Examples exist for:
# 1) Aligning 2 SPIRE observations
# 2) Align multiple SPIRE observations with a reference SPIRE map
# 3) Align a SPIRE Observation with an ancillary observation
# 4) Align a SPIRE Observation with a HIPE source list Product
# Before and after difference maps are produced for comparison
#
# Assumptions: 
# Observations exist in data pools in a Local Store
# 
# Limitations:
#   The current version of astrometryFix does not handle reference images where 
#   the rotation of the horizontal and vertical axes of the reference image 
#   (crot2a keword in image.wcs) are not zero
#
# Authors: Anthony Smith, Ali Dariush, Chris Pearson
# History:  
# 2011/11/4 v1.0: First User Version 
#                    
#
#   Written for HIPE v.8.0
###########################################################################



###########################################################################
###                (1) Align two SPIRE observations                     ###
###########################################################################
#
# ******** LOAD YOUR DATA ************************************************
# REQUIRED:
# obs_ref : first SPIRE observation (reference)
# obs1    : second SPIRE observation (astrometry to be corrected)

# Example below is for a pair of parallel mode observations
# first observation (reference)
myDataPool = "default"         # Name of the data pool for ref image
myObsid    = 0x5000241a                         # Observation ID. e.g.: 1342188781
obs_ref        = getObservation(myObsid, poolName=myDataPool)

# second observation (astrometry to be corrected)
myDataPool1 = "default"       # Name of the data pool for image
myObsid1    = 0x5000241b                        # Observation ID. e.g.: 1342188782
obs1        = getObservation(myObsid1, poolName=myDataPool1)


# ******** ASTROMETRY CORRECTION AND MAP COMPARISON  **********************
# Create a difference map (PSW) before correcting astrometry
differenceMapBefore = imageSubtract(image1=obs_ref.level2.getProduct("PSW"), \
                                 image2=obs1.level2.getProduct("PSW"), ref=1)

# Change the astrometry in obs1 to be consistent with obs_ref
obs1 = astrometryFix(data=obs1, reference=obs_ref)

# Create a difference map (PSW) after correcting astrometry
differenceMapAfter = imageSubtract(image1=obs_ref.level2.getProduct("PSW"), \
                                image2=obs1.level2.getProduct("PSW"), ref=1)

# Display the difference maps before and after the astrometry correction
diffMapBefore01= Display(differenceMapBefore,True,title="(1) Map Difference (before correction)")
diffMapAfter01 = Display(differenceMapAfter ,True,title="(1) Map Difference (after correction)")




###########################################################################
###            (2) Align multiple SPIRE observations                    ###
###########################################################################
#
# ******** LOAD YOUR DATA ************************************************
# REQUIRED:
# obs_ref: first SPIRE observation (reference)
# mapPsw : reference map (PSW) }
# mapPmw : reference map (PMW) } 
# mapPlw : reference map (PLW) }
# obsids : [obsid1, obsid2, ... obsidN] list of obsIDs (astrometry to be corrected)
# pools  : [pool1, pool2, ... poolN] associated list of Pools (astrometry to be corrected)

# Reference observation
myDataPool = "default"       # Name of the data pool for ref image
myObsid    = 0x5000241a                         # Observation ID. e.g.: 1342188781
obs_ref    = getObservation(myObsid, poolName=myDataPool)

# Reference maps (using reference observation products)
mapPsw=obs_ref.level2.getProduct("PSW")
mapPmw=obs_ref.level2.getProduct("PMW")
mapPlw=obs_ref.level2.getProduct("PLW")

# Combine reference maps into a MapContext
mapContext = MapContext()
mapContext.setProduct("PSW", mapPsw)
mapContext.setProduct("PMW", mapPmw)
mapContext.setProduct("PLW", mapPlw)

# List of all Obs IDs with Pools that will be corrected using the reference image
#obsids  = [0x5000241b, \
#           0x5000241c, \
#           0x5000241d]
#pools   = [ "default",
#            "default", \
#            "default"]
obsids  = [0x5000241b]
pools   = [ "default"]
if len(obsids)!=len(pools):
    raise Exception("Warning!! List of OBSIDs and Pools must contain same number of entries")

# Change the astrometry in each observation
differenceMapsBefore = []
differenceMapsAfter = []
for i in range(len(obsids)):
    #l2temp2=MapContext()
    #l2temp2.setProduct("PSW",obs.level2.getProduct("PSW"))
    #l2temp2.setProduct("PMW",obs.level2.getProduct("PMW"))
    #l2temp2.setProduct("PLW",obs.level2.getProduct("PLW"))
    #obs.setLevel2(l2temp2)

    print "Processing ObsID =", obsids[i],"("+hex(obsids[i])+") from Pool ",pools[i]
    # Load observation into variable obs
    obs = getObservation(obsids[i], poolName=pools[i])
    obs.level2.getRefs().remove("PSWRg")
    obs.level2.getRefs().remove("PMWRg")
    obs.level2.getRefs().remove("PLWRg")
    # Create a difference map (PSW) before correcting astrometry
    differenceMapsBefore.append(imageSubtract( \
        image1=mapPsw, image2=obs.level2.getProduct("PSW"), ref=1))
    #
    # Change the astrometry in obs to be consistent with reference maps
    obs = astrometryFix(data=obs, reference=mapContext)
    #
    # Create a difference map (PSW) after correcting astrometry
    differenceMapsAfter.append(imageSubtract( \
        image1=mapPsw, image2=obs.level2.getProduct("PSW"), ref=1))
pass        
    



###########################################################################
###    (3) Align SPIRE observation to an ancillary reference image      ###
###########################################################################
#
# ******** LOAD YOUR DATA ************************************************
# REQUIRED:
# image : ancillary reference image
# obs   : SPIRE observation (astrometry to be corrected)

# ******** GET REFERENCE IMAGE ********************************************                                      
# Example from the NVSS (for SPIRE obsid 1342188781):
# http://www.cv.nrao.edu/nvss/postage.shtml 
# e.g. ra=12:01:10.602 ; dec=+14:06:16.98 (desired image size ~30.0x30.0 arcmin^2)
# Download the NVSS image in fits format (for example nvss.fits)
#
# CAVEAT: The current version of astrometryFix does not handle reference images 
# where the rotation of the horizontal and vertical axes of the reference 
# image (crot2a keword in image.wcs) are not zero !
fits=FitsArchive()
image_ref = fitsReader(file = 'NVSS.fits')
#

# SPIRE observation (astrometry to be corrected)
myDataPool = "ImageDataPoolName"               # Name of the data pool for image
myObsid    = 0000000000                        # Observation ID. e.g.: 1342188781
obs        = getObservation(myObsid, poolName=myDataPool)
image=obs.level2.getProduct("PSW")

# Change the astrometry in obs to be consistent with the reference image
obs = astrometryFix(data=obs, reference=image_ref) 

# Create a difference map (PSW) after correcting astrometry
differenceMapAfter = imageSubtract(image1=image, \
                                   image2=obs.level2.getProduct("PSW"), ref=1)




###########################################################################
###        (4) Align SPIRE observation to known source positions        ###
###########################################################################
#
# ******** LOAD YOUR DATA ************************************************
# REQUIRED:
# sourceList : SourceListProduct containing positions of known sources
# obs        : observation (astrometry to be corrected)

# observation
myDataPool = "ImageDataPoolName"             # Name of the data pool for ref image
myObsid    = 0000000000                      # Observation ID. e.g.: 1342188781
obs        = getObservation(myObsid, poolName=myDataPool)

# sourceList
# The sourceList is produced by using DAOPHOT or Sussextractor tasks.
# The two tasks are listed in the Applicable folder of the Tasks view whenever
# an image is selected in the Variables view.

#Change the astrometry in obs to be consistent with the positions in sourceList
obs = astrometryFix(data=obs, reference=sourceList)


###########################################################################
###                            END OF SCRIPT                            ###
###########################################################################
