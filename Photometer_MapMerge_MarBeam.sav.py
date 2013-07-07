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
###                        SPIRE Map Merging Script                     ###
###########################################################################
# Title: Photometer_MapMerge.py
#
# Purpose:
# This script will merge two or more observations performed in SPIRE Large Map or 
# Small Map or Parallel mode to produce a single map. 
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
# b) Corresponding array of Pool names
# c) output directory for writing the maps FITS files
#
# Assumptions: A data pool has already been created on disk. 
#               The data has already been processed to Level 1
# 
# Original:  script for OT1 Workshop , ESAC, Madrid, March 2011
# History:  
#          2011/03/10 Luca Conversi:
#                     Original version 
#          2011/03/28 Chris Pearson:
#                     User version for HIPE v.7.0 
#          2011/11/07 Chris Pearson:
#                     Added instrument="SPIRE" parameter to getObservation 
#                     User version for HIPE v.8.0 
#          2011/11/11 Chris Pearson:
#                     Updated script to use new Baseline Removal Task 
#
#   Written for HIPE v.8.0
#   Compatitible with HIPE v.6.0,v.7.0
###########################################################################
#
###########################################################################
###                     User Selectable Options                         ###
###########################################################################
# (A) List of OBSIDs in the form of an integer or hexadecimal (0x) number:
#     (Example data taken from OT1 User Workshop, ESAC, March 2011) 
# (B) List of the corresponding data Pools in your Local Store:
# (C) Specify the output directory for writing the maps FITS files:
#
obsids  = [0x50004E4B,0x5000532C]
pool = HttpClientPool ("http://wakefield.bnsc.rl.ac.uk/hcss/pal", "ops")
pools   = [pool,pool]
outDir     = "/home/astrog82/spxcen/Herschel/Calibration/RelGains"
if len(obsids)!=len(pools):
    raise Exception("Warning!! List of OBSIDs and Pools must contain same number of entries")
###########################################################################
 

######################################################################
###        find source position                     ###
##################################################################

print "Reading in ObsID =", obsids[0],"("+hex(obsids[0])+") from Pool ",pools[0]

###get location of objects from first obsid
query=MetaQuery(ObservationContext,"p","p.meta['obsid'].value==%iL"%obsids[0])
store = ProductStorage (pool)    # Just the Data
store.authenticate()
refs=store.select(query)
obs_ref=refs[0].product
#obs_ref = getObservation(obsids[0], poolName=pools[0],instrument="SPIRE")
level1=obs_ref.level1

raEst=148.407835 #RA from JPL Horizons
decEst=14.19442 #Dec from JPL Horizon

##make initial map
##basic baseline removal
scans=baselineRemovalMedian(level1)
mapPsw1Arcsec=naiveScanMapper(scans, array="PSW",resolution=1.)


##Find source position in map
(yEst,xEst)=mapPsw1Arcsec.wcs.getPixelCoordinates(raEst,decEst)
srcPars = sourceFitting(image=mapPsw1Arcsec,\
			minX=xEst-300, minY=yEst-300, width=600, height=600)
xSrc=srcPars["Column1"].data[1]
ySrc=srcPars["Column1"].data[2]
(raSrcMap,decSrcMap)=mapPsw1Arcsec.wcs.getWorldCoordinates(ySrc,xSrc)

##Find source position from timelines (more accurate)
scansList=[]
for i in range(level1.count):
	scansList.append(level1.getProduct(i))
scansContext=ScanContext(scansList)

fitter=TimelineSourceFitterTask()
output=fitter(input=scansContext,array='PSW',\
	sourcePositionEstimate=Double1d([raSrcMap,decSrcMap]),rPeak=22.,\
	fitEllipticalGaussian=False,useBackInFit=False,allowVaryBackground=False)

raSrc=output.meta["fitRa"].value
decSrc=output.meta["fitDec"].value

wcsIn=mapPsw1Arcsec.wcs
#wcsOut=Wcs(crpix1=1000.,crpix2=1000.,crval1=raSrc,crval2=decSrc,\
#    	cdelt1=wcsIn.cdelt1,cdelt2=wcsIn.cdelt2,ctype1=wcsIn.ctype1,ctype2=wcsIn.ctype2,\
#    	equinox=wcsIn.equinox,naxis2=2000,naxis1=2000,crota2=wcsIn.crota2)
wcsOut=Wcs(crpix1=1000.,crpix2=1000.,crval1=raSrc,crval2=decSrc,\
    	cdelt1=wcsIn.cdelt1,cdelt2=wcsIn.cdelt2,ctype1=wcsIn.ctype1,ctype2=wcsIn.ctype2,\
    	equinox=wcsIn.equinox,naxis2=2000,naxis1=2000,crota2=wcsIn.crota2)

print "Source located at (RA,Dec)=(%.2f, %.2f)"%(raSrc,decSrc)

###				   Finished Source finding				 ###
###########################################################################

# Reference maps (using reference observation products)
mapPsw=obs_ref.level2.getProduct("PSW")
mapPmw=obs_ref.level2.getProduct("PMW")
mapPlw=obs_ref.level2.getProduct("PLW")

# Combine reference maps into a MapContext
mapContext = MapContext()
mapContext.setProduct("PSW", mapPsw)
mapContext.setProduct("PMW", mapPmw)
mapContext.setProduct("PLW", mapPlw)

scanLines = 0
differenceMapsBefore=[]
differenceMapsAfter=[]
# Start cycle on observations
for i in range(len(obsids)):

	if i == 0:
		##copy level 1 from obs_ref
		scanLines = obs_ref.level1
		print obs_ref.level1.count, "Scan lines from  observation ", obsids[i]
	else:
		##Don't need to coalign first obsid with itself!
		print "Co-aligning ObsID "+hex(obsids[i])+" with ObsId "+hex(obsids[0])
		#
		# Load observation into variable obs
		query=MetaQuery(ObservationContext,"p","p.meta['obsid'].value==%iL"%obsids[0])
		store = ProductStorage (pool)    # Just the Data
		store.authenticate()
		refs=store.select(query)
		obs=refs[0].product
		#obs = getObservation(obsids[i], poolName=pools[i],instrument="SPIRE")
		##remove Rg maps (crashed astrometryFix)
		obs.level2.getRefs().remove("PSWRg")
		obs.level2.getRefs().remove("PMWRg")
		obs.level2.getRefs().remove("PLWRg")
		# Create a difference map (PSW) before correcting astrometry
		differenceMapsBefore.append(imageSubtract( \
			image1=mapPsw, image2=obs.level2.getProduct("PSW"), ref=1))
		#
		# Change the astrometry in obs to be consistent with reference maps
		obs = astrometryFix(data=obs, reference=mapContext)
		raPsw0=obs.level2.refs["PSW"].meta["astrometryFixOldCrval1"].getValue()
		raPsw1=obs.level2.refs["PSW"].meta["astrometryFixNewCrval1"].getValue()
		decPsw0=obs.level2.refs["PSW"].meta["astrometryFixOldCrval2"].getValue()
		decPsw1=obs.level2.refs["PSW"].meta["astrometryFixNewCrval2"].getValue()
		print 'PSW Offset: (%.2f , %.2f) arcsec'%((raPsw1-raPsw0)*3600.,(decPsw1-decPsw0)*3600.)

		raPmw0=obs.level2.refs["PMW"].meta["astrometryFixOldCrval1"].getValue()
		raPmw1=obs.level2.refs["PMW"].meta["astrometryFixNewCrval1"].getValue()
		decPmw0=obs.level2.refs["PMW"].meta["astrometryFixOldCrval2"].getValue()
		decPmw1=obs.level2.refs["PMW"].meta["astrometryFixNewCrval2"].getValue()
		print 'PMW Offset: (%.2f , %.2f) arcsec'%((raPsw1-raPsw0)*3600.,(decPsw1-decPsw0)*3600.)

		raPlw0=obs.level2.refs["PLW"].meta["astrometryFixOldCrval1"].getValue()
		raPlw1=obs.level2.refs["PLW"].meta["astrometryFixNewCrval1"].getValue()
		decPlw0=obs.level2.refs["PLW"].meta["astrometryFixOldCrval2"].getValue()
		decPlw1=obs.level2.refs["PLW"].meta["astrometryFixNewCrval2"].getValue()
		print 'PLW Offset: (%.2f , %.2f) arcsec'%((raPlw1-raPlw0)*3600.,(decPlw1-decPlw0)*3600.)
		#
		# Create a difference map (PSW) after correcting astrometry
		differenceMapsAfter.append(imageSubtract( \
		image1=mapPsw, image2=obs.level2.getProduct("PSW"), ref=1))

		## obs now contains data shifted to same as obsid[0]
		print "ObsId ",hex(obsids[i])," co-aligned."
	
		#
		# Display PSW map of the input obseravtions
		# To display a different band, change 2 instances of "PLW" in the line below
		#Display(obs.refs["level2"].product.refs["PLW"].product, title= "%i PLW Map"%obs.obsid)
		#
			
		for ref in obs.level1.refs:
			# Attach level 1 product from next observation.
			scanLines.addRef(ref)
                print obs.level1.count, "Scan lines from  observation ", obsids[i]

##################################################
###      Baseline removal                   ###
##################################################

##set mask
roi=SkyMaskCircle(raSrc,decSrc,4.).not()
#Display(roi.masks(mapPsw1Arcsec))

# Using new Level 1 context. Run baseline removal  as an input to the map making
scans=baselineRemovalMedian(input=scanLines,roi=roi)

# Create combined map
psw_combined = naiveScanMapper(scans, array="PSW",resolution=1.,wcs=wcsOut)
pmw_combined = naiveScanMapper(scans, array="PMW",resolution=1.,wcs=wcsOut)
plw_combined = naiveScanMapper(scans, array="PLW",resolution=1.,wcs=wcsOut)


# Display combined map
Display(psw_combined, title='Combined PSW Map')
Display(pmw_combined, title='Combined PMW Map')
Display(plw_combined, title='Combined PLW Map')

# Save Maps to output directory
simpleFitsWriter(psw_combined, outDir+"MarBeam_psw_1arcsec_combined.fits")
simpleFitsWriter(pmw_combined, outDir+"MarBeam_pmw_1arcsec_combined.fits")
simpleFitsWriter(plw_combined, outDir+"MarBeam_plw_1arcsec_combined.fits")
print "Map saved as FITS files to %s"%(outDir)

# End of Script
###########################################################################



