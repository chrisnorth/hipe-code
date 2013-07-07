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
##			mSPIRE Large Map Mode User Reprocessing Script			###
###########################################################################
#  Purpose:  A simplified version of SPIRE Large Map Mode POF5 pipeline   
#			script distributed with HIPE 6.0.  This is for data reprocessing
#			 by a user using the latest SPIRE calibration products.
# 
#			The results are three FITS files with extensions for the PSW,  
#			PMW, PLW arrays containing the final,
#			- image map  
#			- error map  
#			- coverage map  
#	
#  Usage:	The user needs to specify the options in the simple user input 
#			section at the beginning of the script;
#			- Observation ID (obsid)  
#			- Data pool name  
#			- Out put directory for final fits files  
#
# Note that it is possible to save entire observation back to a pool by
#	 uncommenting the saveObservation command at the end of the script
#
#  Assumptions: A data pool has already been created on disk. 
#			   The data has already been processed to Level 0.5
#
#  Updated: 12/10/2011
#
###########################################################################



###########################################################################
###					 User Selectable Options						 ###
###########################################################################
# (A) Specific OBSID in the form of an integer or hexadecimal (0x) number: 
# (B) the name of the data Pool in your Local Store:
# (C) Specify the output directory for writing the maps FITS files:
#
myObsid	= 0x5000241d
myDataPool = "default"
outDir	 = "/home/astrog82/spxcen/Herschel/Calibration/RelGains"
# e.g.
#myObsid	=  0x50001833
#myDataPool = "OD117-ScanNGC5315-0x50001833"
#outDir	 = "/Users/cpearson/jython/localstore/plots/"
#
# Additional Options
# (D) The mapping Algorithm to use naive or madmap
mapping='naive'
###########################################################################



#*************************************************************************
##  Load in an observation context from your data pool into HIPE:
obs=getObservation(myObsid,poolName=myDataPool)
#obs=obsid_1342213459
print
print "Processing observation %i (0x%X) from data pool %s."%(myObsid, myObsid, myDataPool)


#*************************************************************************
# Calibration Context and Calibration Files 
# Read the latest calibration tree relevant to HCSS v6 from the local disc:
# cal = spireCal(pool="spire_cal_7_0")
cal = spireCal(jarFile='/home/user/spxcen/dp-spire/cal/spire_cal_8_phot_27Sep2011.jar')
# TO CORRECT AN ERROR ON THE ABOVE LINE, run the following two lines 
# ONCE only to load the calibration tree from the Archive (may take some time):
# (for more details, see the "Calibration" chapter in the SPIRE Data Reduction Guide)


#cal = spireCal(calTree="spire_cal_7_0")
#localStoreWriter(cal, "spire_cal_7_0")

# Attach the updated calibration tree to the observation context
##obs.calibration.update(cal)


# Extract necessary Calibration Products from the Observation Context
bsmPos			 = cal.phot.bsmPos
lpfPar			 = cal.phot.lpfPar
detAngOff		  = cal.phot.detAngOff
chanTimeConst	  = cal.phot.chanTimeConst
chanNum			= cal.phot.chanNum
fluxConvList	   = cal.phot.fluxConvList
tempDriftCorrList  = cal.phot.tempDriftCorrList
relGains	   = cal.phot.chanRelGain
# Note to read in a single calibration fits file from some location use, e.g., 
# bsmPos = simpleFitsReader("/enter/path/here/"+"YourCalibrationFilename.fits")

# Extract the necessary Auxiliary Products from the Observation Context
hpp	  = obs.auxiliary.pointing
siam	 = obs.auxiliary.siam
timeCorr = obs.auxiliary.timeCorrelation
#*************************************************************************


#*************************************************************************
## Reports how many scan lines there are to process
count=1
bbids=obs.level0_5.getBbids(0xa103)
nlines=len(bbids)
print "Total number of scan lines: ",nlines
print
# Create Level1 context to collect Level one products
level1=Level1Context(obs.meta['obsid'].value)
level1_rg=Level1Context(obs.meta['obsid'].value)

#
###########################################################################
###   Pipeline Level 0.5 to Level 1									 ###
###   Process all Building Blocks for this observation				  ###
###########################################################################
# Loop over scan lines
print 'bbids:',bbids
for bbid in bbids:
	print "Starting BBID=0x%x: scan %i / %i"%(bbid,count,nlines)
	# Get basic level 0.5 data products (detector data and housekeeping data) 
	pdt  = obs.level0_5.get(bbid).pdt
	nhkt = obs.level0_5.get(bbid).nhkt
	# record the calibration tree version used by the pipeline
	pdt.calVersion = cal.version
	#
	# -----------------------------------------------------------
	# (1) join all scan legs and turnarounds together
	bbCount=bbid & 0xFFFF
	pdtLead=None
	nhktLead=None
	pdtTrail=None
	nhktTrail=None
	if bbCount >1:
		blockLead=obs.level0_5.get(0xaf000000L+bbCount-1)
		pdtLead=blockLead.pdt
		nhktLead=blockLead.nhkt
		if pdtLead != None and pdtLead.sampleTime[-1] < pdt.sampleTime[0]-3.0:
			pdtLead=None
			nhktLead=None
	if bbid < MAX(Long1d(bbids)):
		blockTrail=obs.level0_5.get(0xaf000000L+bbCount)
		pdtTrail=blockTrail.pdt
		nhktTrail=blockTrail.nhkt
		if pdtTrail != None and pdtTrail.sampleTime[0] > pdt.sampleTime[-1]+3.0:
			pdtTrail=None
			nhktTrail=None
	pdt=joinPhotDetTimelines(pdt,pdtLead,pdtTrail)
	nhkt=joinNhkTimelines(nhkt,nhktLead,nhktTrail)
	#
	# -----------------------------------------------------------
	# (2) Convert BSM timeline to angles on sky (constant for scan map)
	bat=calcBsmAngles(nhkt,bsmPos=bsmPos)
	#
	# -----------------------------------------------------------
	# (3) Create the Spire Pointing Product for this observation
	spp=createSpirePointing(detAngOff=detAngOff,bat=bat,hpp=hpp,siam=siam)
	#
	# -----------------------------------------------------------
	# (4) Detect jumps in the Thermistor timelines that occasionally occur,
	# leading to map artefacts introduced in the temperature drift correction
	# Also requires the Temperature Drift Correct Calibration File.
	tempDriftCorr=tempDriftCorrList.getProduct(pdt.meta["biasMode"].value,pdt.startDate)
	if pdt.meta["biasMode"].value == "nominal":
		pdt=signalJumpDetector(pdt,tempDriftCorr=tempDriftCorr, kappa=2.0,gamma=6.0,\
			gapWidth=1.0,windowWidth=40.0, filterType="DISCRETEDERIVATIVE",glitchinfo="NULL")
	#
	# -----------------------------------------------------------
	# (5) Run the concurrent deglitcher on the timeline data
	pdt=concurrentGlitchDeglitcher(pdt,chanNum=chanNum,kappa=2.0, size = 15, correctGlitches = True)
	#
	# -----------------------------------------------------------
	# (6) Run the wavelet deglitcher on the timeline data
	pdt=waveletDeglitcher(pdt, scaleMin=1.0, scaleMax=8.0, scaleInterval=7, holderMin=-3.0,\
		holderMax=-0.3, correlationThreshold=0.3, optionReconstruction='linear',\
		reconstructionPointsBefore=1, reconstructionPointsAfter=6)
	#
	# Alternatively, run the sigma-kappa deglitcher.
	# The following task can be uncommented to try the sigma-kappa deglitcher.
	# In this case the wavelet deglitcher should be commented out.
	# This module is still a prototype and should be used with caution.
	#pdt = sigmaKappaDeglitcher(pdt,
	#			filterType="BOXCAR", boxFilterWidth = 3, \
	#			boxFilterCascade = 1, kappa = 4.0, \
	#			disableSigmaKappaDetection = 'NULL', \
	#			largeGlitchMode = 'ADDITIVE', \
	#			largeGlitchDiscriminatorTimeConstant = 4, \
	#			largeGlitchRemovalTimeConstant = 6, \
	#			disableLargeGlitchDetection = 'NULL', \
	#			correctionMode = 'DIRECT', gamma = 1.0, \
	#			randomSeed = 1984574303L, \
	#			disableGlitchReconstruction = 'NULL', \
	#			iterationNumber = 1)
	#
	# -----------------------------------------------------------
	# (7) Apply the Low Pass Filter response correction
	pdt=lpfResponseCorrection(pdt,lpfPar=lpfPar)
	#
	# -----------------------------------------------------------
	# (8) Apply the flux conversion 
	fluxConv=fluxConvList.getProduct(pdt.meta["biasMode"].value,pdt.startDate)
	pdt=photFluxConversion(pdt,fluxConv=fluxConv)
	#
	# -----------------------------------------------------------
	# (9) Make the temperature drift correction
	pdt=temperatureDriftCorrection(pdt,tempDriftCorr=tempDriftCorr)
	#
	# -----------------------------------------------------------
	# (10) Apply the bolometer time response correction
	pdt=bolometerResponseCorrection(pdt,chanTimeConst=chanTimeConst)
	#
	# -----------------------------------------------------------
	# (11) Add pointing timelines to the data
	psp=associateSkyPosition(pdt,spp=spp)
	#
	# -----------------------------------------------------------
	# (12) Cut the timeline back into individual scan lines .
	# NOTE: If you want include turnaround data in map making, call the following
	# task with the option "extend=True"
	psp=cutPhotDetTimelines(psp,extend=False)
	#
	# -----------------------------------------------------------
	# (13) Apply the time correlation 
	psp=timeCorrelation(psp,timeCorr)
	#
	#
	#
	#####Apply relative gains
	psp_rg = applyRelativeGains(psp, relGains)
	#
	# -----------------------------------------------------------
	# Add Photometer Scan Product to Level 1 context
	# Both with AND without Relative Gains applied
	level1.addProduct(psp)
	level1_rg.addProduct(psp_rg)
	#
	# -----------------------------------------------------------
	# Save scans to disk (with & without relative gains)
	dirPath = '/home/astrog82/spxcen/Herschel/Calibration/RelGains/Obs/'
	#
	#without relative gains
	#filename = 'psp_'+hex(obs.meta['obsid'].value)[2:10]+"_scan"+str(bbCount)+'.fits'
	#print dirPath+filename
	#FitsArchive().save(dirPath+filename, psp)
	#
	# with relative gains
	#filename_rg = 'psp_'+hex(obs.meta['obsid'].value)[2:10]+"_scan"+str(bbCount)+'_RelGains.fits'
	#print dirPath+filename
	#FitsArchive().save(dirPath+filename_rg, psp_rg)
	
	print "Completed BBID=0x%x (scan %i/%i)"%(bbid,count+1,nlines)
	# set the progress
	count=count+1
print
print "Finished the Level 0.5 to Level 1 processing for OBSID= %i, (0x%x)"%(myObsid,myObsid)
# Update the Level 1 Context in the Observation Context
obs.level1 = level1
#
###			Finished the Level 0.5 to Level 1 processing			 ###
###########################################################################

######################################################################
###        find source position                     ###
##################################################################
raEst=326.026 #RA from JPL Horizons
decEst=-14.064 #Dec from JPL Horizons

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

###				   Finished Source finding				 ###
###########################################################################

###########################################################################
###					  Baseline Subtraction						   ###
###########################################################################
print
print "Begin the Baseline Subtraction for OBSID= %i, (0x%x)"%(myObsid,myObsid)

##set mask
roi=SkyMaskCircle(raSrc,decSrc,4.).not()
#Display(roi.masks(mapPsw1Arcsec))

#
# Using Level 1 context. Run baseline removal  as an input to the map making
scans=baselineRemovalMedian(level1,roi=roi)
scans_rg=baselineRemovalMedian(level1_rg,roi=roi)
print "Finished the Baseline Subtraction for OBSID= %i, (0x%x)"%(myObsid,myObsid)
print
###				   Finished the Baseline Subtraction				 ###
###########################################################################

# --------------------------------------
# save baseline subtracted scans to disk
dirPath = '/home/astrog82/spxcen/Herschel/Calibration/RelGains/Obs/'
#counter=0
#for i in range (scans.count):
#	counter+=1
#	filename = 'psp_'+hex(obs.meta['obsid'].value)[2:10]+"_scan"+str(counter)+'_BLsub.fits'
#	print dirPath+filename
#	FitsArchive().save(dirPath+filename, scans.getProduct(i))
#	#
#	filename_rg = 'psp_'+hex(obs.meta['obsid'].value)[2:10]+"_scan"+str(counter)+'_RelGains_BLsub.fits'
#	FitsArchive().save(dirPath+filename_rg, scans_rg.getProduct(i))

#pspList=[]
#loop
#	pspList.append(FitsArchive().load(path))
	



##scanCon = ScanContext(pspList)
###mapPsw=naiveScanMapper(scanCon, array="PSW")


###########################################################################
###						  Mapmaking								  ###
###########################################################################
# 
# Either the Naive Map maker
print 'beginning Map Maker'
mapping = 'naive'
if mapping == 'naive':
	print 'Starting Naive Map maker'
	mapPlw=naiveScanMapper(scans_rg, array="PLW",resolution=1.,wcs=wcsOut)
	mapPmw=naiveScanMapper(scans_rg, array="PMW",resolution=1.,wcs=wcsOut)
	mapPsw=naiveScanMapper(scans_rg, array="PSW",resolution=1.,wcs=wcsOut)
	#mapPlw_rg=naiveScanMapper(scans_rg, array="PLW")
	#mapPmw_rg=naiveScanMapper(scans_rg, array="PMW")
	#mapPsw_rg=naiveScanMapper(scans_rg, array="PSW")
else:
	# -----------------------------------------------------------
	# or the Mad Map Map maker (requires Channel Noise Table Calibration Product)
	print 'Starting Mad Mapper'
	chanNoise=cal.phot.chanNoiseList.getProduct(level1.getProduct(0).meta["biasMode"].value,level1.getProduct(0).startDate)
	mapPlw=madScanMapper(scans, array="PLW",chanNoise=chanNoise)
	mapPmw=madScanMapper(scans, array="PMW",chanNoise=chanNoise)
	mapPsw=madScanMapper(scans, array="PSW",chanNoise=chanNoise)
pass
# Update the Level 2 (map) Context in the Observation Context
obs.level2.setProduct("PSW", mapPsw)
obs.level2.setProduct("PMW", mapPmw)
obs.level2.setProduct("PLW", mapPlw)

#obs.level2.setProduct("PSW", mapPsw_rg)
#obs.level2.setProduct("PMW", mapPmw_rg)
#obs.level2.setProduct("PLW", mapPlw_rg)

#
print "Finished the map making for OBSID= %i, (0x%x)"%(myObsid,myObsid)

print
#
#
# -----------------------------------------------------------
# Save Maps to output directory
simpleFitsWriter(mapPsw, "%smapPSW_%s.fits"%(dirPath, str(hex(obs.meta['obsid'].value)[2:10])))
simpleFitsWriter(mapPmw, "%smapPMW_%s.fits"%(dirPath, str(hex(obs.meta['obsid'].value)[2:10])))
simpleFitsWriter(mapPlw, "%smapPLW_%s.fits"%(dirPath, str(hex(obs.meta['obsid'].value)[2:10])))
#simpleFitsWriter(mapPsw_rg, "%smapPSW_rg_%s.fits"%(dirPath, str(hex(obs.meta['obsid'].value)[2:10])))
#simpleFitsWriter(mapPmw_rg, "%smapPMW_rg_%s.fits"%(dirPath, str(hex(obs.meta['obsid'].value)[2:10])))
#simpleFitsWriter(mapPlw_rg, "%smapPLW_rg_%s.fits"%(dirPath, str(hex(obs.meta['obsid'].value)[2:10])))
print "Map saved as FITS files to %s"%(dirPath)
#
###				Finished the Mapmaking								###
############################################################################

# Finally we can save the new reprocessed observation back to your hard disk
# Uncomment the next line and choose a poolName, either the existing one or a new one
saveObservation(obs,poolName=myDataPool,saveCalTree=True)

#
print
print "Completed the processing of OBSID= %i, (0x%x)"%(myObsid,myObsid)



#### End of the script ####
